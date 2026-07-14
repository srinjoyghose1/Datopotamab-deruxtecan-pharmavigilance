"""FAERS comparator sensitivity analyses for the breast-cancer Dato-DXd cohort.

Mirrors scripts/07_jader_comparator.py with three comparators over the same
2025Q1-2026Q1 window:

1. all deduplicated FAERS reports, excluding every Dato-DXd report;
2. breast-cancer reports with a prespecified active comparator as Primary Suspect;
3. (2), excluding reports with sacituzumab govitecan or trastuzumab deruxtecan
   as a Primary Suspect.

Indications are linked to the relevant drug through PRIMARYID + DRUG_SEQ. Thus a
breast-cancer indication belonging to a concomitant medication cannot qualify a
case. Breast cancer includes all explicitly reported subtypes, including TNBC.
"""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma


REPO_ROOT = Path(__file__).resolve().parents[1]
EXTRACT_ROOT = REPO_ROOT / "data" / "raw" / "faers_quarterly" / "extracted"
RAW_FAERS_DIR = REPO_ROOT / "data" / "raw" / "faers"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
TABLES_DIR = REPO_ROOT / "outputs" / "tables"
QUARTERS = ["2025q1", "2025q2", "2025q3", "2025q4", "2026q1"]

ACTIVE_COMPARATORS = {
    "capecitabine": ["capecitabine"],
    "eribulin": ["eribulin"],
    "gemcitabine": ["gemcitabine"],
    "vinorelbine": ["vinorelbine"],
    "paclitaxel": ["paclitaxel"],
    "nab_paclitaxel": ["nab-paclitaxel", "paclitaxel protein-bound", "abraxane"],
    "carboplatin": ["carboplatin"],
    "sacituzumab_govitecan": ["sacituzumab govitecan", "trodelvy"],
    "trastuzumab_deruxtecan": ["trastuzumab deruxtecan", "enhertu"],
}
EXCLUDED_ADC_COMPARATORS = {"sacituzumab_govitecan", "trastuzumab_deruxtecan"}


def find_table_file(quarter: str, table: str) -> Path:
    matches = list((EXTRACT_ROOT / quarter / "ASCII").glob(f"{table}*.txt"))
    matches += list((EXTRACT_ROOT / quarter / "ASCII").glob(f"{table}*.TXT"))
    return matches[0]


def load_deleted_caseids() -> set[str]:
    deleted: set[str] = set()
    for path in EXTRACT_ROOT.glob("*/Deleted/DELETE*.txt"):
        deleted.update(line.strip() for line in path.read_text().splitlines() if line.strip())
    return deleted


def background_case_universe() -> pd.DataFrame:
    frames = []
    for quarter in QUARTERS:
        frame = pd.read_csv(
            find_table_file(quarter, "DEMO"), delimiter="$",
            usecols=["primaryid", "caseid", "fda_dt"], dtype=str,
            low_memory=False, encoding="latin1",
        )
        frame.columns = [column.strip().lower() for column in frame.columns]
        frames.append(frame)
    demo = pd.concat(frames, ignore_index=True)
    demo["_fda"] = pd.to_numeric(demo["fda_dt"], errors="coerce").fillna(-1)
    demo["_primary"] = pd.to_numeric(demo["primaryid"], errors="coerce").fillna(-1)
    demo = demo.sort_values(["_fda", "_primary"], ascending=False)
    demo = demo.drop_duplicates("caseid", keep="first")
    demo = demo[~demo["caseid"].isin(load_deleted_caseids())]
    return demo[["primaryid", "caseid"]]


def name_mask(frame: pd.DataFrame, terms: list[str]) -> pd.Series:
    combined = frame["drugname"].fillna("") + " " + frame["prod_ai"].fillna("")
    combined = combined.str.casefold()
    mask = pd.Series(False, index=frame.index)
    for term in terms:
        mask |= combined.str.contains(term.casefold(), regex=False)
    return mask


def target_breast_keys() -> tuple[set[str], set[str]]:
    """Return explicitly breast-cancer Dato cases and every Dato case."""
    cases = pd.read_csv(PROCESSED_DIR / "analysis_cases.csv", dtype=str)
    all_dato = set(cases["primaryid"])
    drug = pd.read_csv(RAW_FAERS_DIR / "drug.csv", dtype=str)
    indi = pd.read_csv(RAW_FAERS_DIR / "indi.csv", dtype=str)
    links = drug[["primaryid", "drug_seq"]].drop_duplicates().merge(
        indi, left_on=["primaryid", "drug_seq"],
        right_on=["primaryid", "indi_drug_seq"], how="inner",
    )
    breast = links["indi_pt"].fillna("").str.contains("breast", case=False, regex=False)
    return set(links.loc[breast, "primaryid"]) & all_dato, all_dato


def active_comparator_keys(kept_primaryids: set[str]) -> tuple[set[str], set[str]]:
    active: set[str] = set()
    excluded_adc: set[str] = set()
    for quarter in QUARTERS:
        drug = pd.read_csv(
            find_table_file(quarter, "DRUG"), delimiter="$",
            usecols=["primaryid", "drug_seq", "role_cod", "drugname", "prod_ai"],
            dtype=str, low_memory=False, encoding="latin1",
        )
        drug.columns = [column.strip().lower() for column in drug.columns]
        drug = drug[drug["primaryid"].isin(kept_primaryids) & drug["role_cod"].eq("PS")].copy()
        matched_frames = []
        for label, terms in ACTIVE_COMPARATORS.items():
            matched = drug[name_mask(drug, terms)].copy()
            if matched.empty:
                continue
            matched["comparator"] = label
            matched_frames.append(matched)
        if not matched_frames:
            continue
        candidates = pd.concat(matched_frames, ignore_index=True)

        indi = pd.read_csv(
            find_table_file(quarter, "INDI"), delimiter="$",
            usecols=["primaryid", "indi_drug_seq", "indi_pt"],
            dtype=str, low_memory=False, encoding="latin1",
        )
        indi.columns = [column.strip().lower() for column in indi.columns]
        linked = candidates.merge(
            indi, left_on=["primaryid", "drug_seq"],
            right_on=["primaryid", "indi_drug_seq"], how="inner",
        )
        linked = linked[linked["indi_pt"].fillna("").str.contains("breast", case=False, regex=False)]
        active.update(linked["primaryid"])
        excluded_adc.update(linked.loc[linked["comparator"].isin(EXCLUDED_ADC_COMPARATORS), "primaryid"])
    return active, excluded_adc


def reaction_counts(case_sets: dict[str, set[str]]) -> dict[str, Counter]:
    counters = {name: Counter() for name in case_sets}
    for quarter in QUARTERS:
        frame = pd.read_csv(
            find_table_file(quarter, "REAC"), delimiter="$",
            usecols=["primaryid", "pt"], dtype=str,
            low_memory=False, encoding="latin1",
        )
        frame.columns = [column.strip().lower() for column in frame.columns]
        frame = frame.drop_duplicates(["primaryid", "pt"])
        for name, keys in case_sets.items():
            counters[name].update(frame.loc[frame["primaryid"].isin(keys), "pt"].dropna())
    return counters


def ror_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    if min(a, b, c, d) == 0:
        a, b, c, d = (value + 0.5 for value in (a, b, c, d))
    ror = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return ror, math.exp(math.log(ror) - 1.96 * se), math.exp(math.log(ror) + 1.96 * se)


def prr_chi2(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    raw = [[a, b], [c, d]]
    if min(a, b, c, d) == 0:
        a, b, c, d = (value + 0.5 for value in (a, b, c, d))
    prr = (a / (a + b)) / (c / (c + d))
    table = raw if min(value for row in raw for value in row) > 0 else [[a, b], [c, d]]
    chi2, _, _, _ = stats.chi2_contingency(table, correction=True)
    return prr, chi2


def bcpnn_ic(a: int, n_target: int, n_pt: int, n_total: int) -> tuple[float, float, float]:
    expected = n_target * n_pt / n_total
    ic = np.log2((a + 0.5) / (expected + 0.5))
    variance = (1 / np.log(2) ** 2) * (1 / (a + 0.5) + 1 / (expected + 0.5))
    return ic, ic - 1.96 * np.sqrt(variance), expected


def fit_gamma_poisson(a: np.ndarray, expected: np.ndarray) -> tuple[float, float]:
    ratios = a / expected
    mean, variance = ratios.mean(), ratios.var(ddof=1)
    if variance <= 0:
        return 1.0, 1.0
    beta = mean / variance
    return mean * beta, beta


def calculate_signals(target: Counter, comparator: Counter, n_target: int, n_comparator: int, name: str) -> pd.DataFrame:
    rows = []
    for pt, a in target.items():
        b, c = n_target - a, comparator.get(pt, 0)
        d = n_comparator - c
        ror, lo, hi = ror_ci(a, b, c, d)
        prr, chi2 = prr_chi2(a, b, c, d)
        ic, ic025, expected = bcpnn_ic(a, n_target, a + c, n_target + n_comparator)
        rows.append({
            "pt": pt, "comparator": name, "a_dato_cases": a,
            "b_dato_without_pt": b, "c_comparator_cases": c,
            "d_comparator_without_pt": d, "n_dato": n_target,
            "n_comparator": n_comparator, "ror": ror,
            "ror_ci_lower": lo, "ror_ci_upper": hi, "prr": prr,
            "chi2_yates": chi2, "ic": ic, "ic025": ic025,
            "expected_count": expected,
        })
    out = pd.DataFrame(rows)
    alpha, beta = fit_gamma_poisson(out["a_dato_cases"].to_numpy(), out["expected_count"].to_numpy())
    out["ebgm"] = [float(np.exp(digamma(alpha + a) - np.log(beta + e))) for a, e in zip(out.a_dato_cases, out.expected_count)]
    out["eb05"] = [float(stats.gamma.ppf(0.05, a=alpha + a, scale=1 / (beta + e))) for a, e in zip(out.a_dato_cases, out.expected_count)]
    out["ror_signal"] = (out.a_dato_cases >= 3) & (out.ror_ci_lower > 1)
    out["prr_signal"] = (out.a_dato_cases >= 3) & (out.prr >= 2) & (out.chi2_yates >= 4)
    out["bcpnn_signal"] = (out.a_dato_cases >= 3) & (out.ic025 > 0)
    out["mgps_signal"] = (out.a_dato_cases >= 3) & (out.eb05 >= 2)
    out["consensus_signal"] = out[["ror_signal", "prr_signal", "bcpnn_signal"]].all(axis=1)
    return out.sort_values(["consensus_signal", "ror", "a_dato_cases"], ascending=[False, False, False])


def main() -> None:
    background = background_case_universe()
    kept = set(background["primaryid"])
    target, all_dato = target_breast_keys()
    active, excluded_adc = active_comparator_keys(kept)
    active -= all_dato
    class_excluded = active - excluded_adc
    sets = {
        "dato_breast_all_subtypes": target,
        "full_faers_excluding_all_dato": kept - all_dato,
        "active_breast_comparator": active,
        "active_breast_comparator_class_exclusion": class_excluded,
    }
    counts = reaction_counts(sets)

    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame([{"cohort": name, "n_cases": len(keys)} for name, keys in sets.items()])
    summary.to_csv(TABLES_DIR / "faers_comparator_cohort_sizes.csv", index=False)
    outputs = []
    for name in list(sets)[1:]:
        result = calculate_signals(counts["dato_breast_all_subtypes"], counts[name], len(target), len(sets[name]), name)
        result.to_csv(TABLES_DIR / f"faers_signals_{name}.csv", index=False)
        outputs.append(result)
    pd.concat(outputs, ignore_index=True).to_csv(TABLES_DIR / "faers_signals_comparator_all.csv", index=False)

    print(summary.to_string(index=False))
    for result in outputs:
        name = result["comparator"].iloc[0]
        print(f"\n{name}: consensus signals")
        columns = ["pt", "a_dato_cases", "ror", "ror_ci_lower", "prr", "ic025", "eb05"]
        print(result.loc[result.consensus_signal, columns].head(25).to_string(index=False))


if __name__ == "__main__":
    main()
