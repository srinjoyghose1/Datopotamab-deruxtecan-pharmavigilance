"""FAERS comparator hierarchy for the breast-cancer Dato-DXd cohort.

Primary tiers over 2025Q1-2026Q1 are the full FAERS background, the four
TROPION-Breast01 trial-aligned chemotherapies, trastuzumab deruxtecan alone, and
sacituzumab govitecan alone. Supplemental tiers retain the expanded chemotherapy
and legacy pooled active-comparator designs. Leave-one-drug-out analyses test the
stability of primary trial-aligned signals.

Indications are linked to the relevant drug through PRIMARYID + DRUG_SEQ. Thus a
breast-cancer indication belonging to a concomitant medication cannot qualify a
case. To approximate the HR-positive/HER2-negative population without discarding
reports whose receptor status is simply missing, explicitly discordant TNBC,
HER2-positive, and hormone-receptor-negative indications are excluded at case
level; nonspecific breast-cancer indications remain eligible.
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
TRIAL_ALIGNED_DRUGS = ("capecitabine", "eribulin", "gemcitabine", "vinorelbine")
EXPANDED_CHEMO_DRUGS = TRIAL_ALIGNED_DRUGS + ("paclitaxel", "nab_paclitaxel", "carboplatin")
EXCLUDED_ADC_COMPARATORS = {"sacituzumab_govitecan", "trastuzumab_deruxtecan"}
DISCORDANT_INDICATION_TERMS = (
    "triple negative",
    "triple-negative",
    "her2 positive",
    "her2-positive",
    "hormone receptor negative",
    "hormone receptor-negative",
    "hr negative",
    "hr-negative",
)


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


def eligible_breast_ids(frame: pd.DataFrame) -> set[str]:
    """Return breast-cancer report IDs after case-level subtype exclusions."""
    indication = frame["indi_pt"].fillna("").str.casefold()
    breast = indication.str.contains("breast", regex=False)
    breast_ids = set(frame.loc[breast, "primaryid"])
    discordant = pd.Series(False, index=frame.index)
    for term in DISCORDANT_INDICATION_TERMS:
        discordant |= indication.str.contains(term, regex=False)
    discordant_ids = set(frame.loc[breast & discordant, "primaryid"])
    return breast_ids - discordant_ids


def target_breast_keys() -> tuple[set[str], set[str]]:
    """Return HR+/HER2- compatible breast-cancer Dato cases and all Dato cases."""
    cases = pd.read_csv(PROCESSED_DIR / "analysis_cases.csv", dtype=str)
    all_dato = set(cases["primaryid"])
    drug = pd.read_csv(RAW_FAERS_DIR / "drug.csv", dtype=str)
    indi = pd.read_csv(RAW_FAERS_DIR / "indi.csv", dtype=str)
    links = drug[["primaryid", "drug_seq"]].drop_duplicates().merge(
        indi, left_on=["primaryid", "drug_seq"],
        right_on=["primaryid", "indi_drug_seq"], how="inner",
    )
    return eligible_breast_ids(links) & all_dato, all_dato


def comparator_keys_by_drug(kept_primaryids: set[str]) -> dict[str, set[str]]:
    """Build indication-compatible PS report sets for every comparator drug."""
    keys_by_drug = {label: set() for label in ACTIVE_COMPARATORS}
    for quarter in QUARTERS:
        drug = pd.read_csv(
            find_table_file(quarter, "DRUG"), delimiter="$",
            usecols=["primaryid", "drug_seq", "role_cod", "drugname", "prod_ai"],
            dtype=str, low_memory=False, encoding="latin1",
        )
        drug.columns = [column.strip().lower() for column in drug.columns]
        drug = drug[drug["primaryid"].isin(kept_primaryids) & drug["role_cod"].eq("PS")].copy()
        indi = pd.read_csv(
            find_table_file(quarter, "INDI"), delimiter="$",
            usecols=["primaryid", "indi_drug_seq", "indi_pt"],
            dtype=str, low_memory=False, encoding="latin1",
        )
        indi.columns = [column.strip().lower() for column in indi.columns]
        for label, terms in ACTIVE_COMPARATORS.items():
            matched = drug[name_mask(drug, terms)]
            if matched.empty:
                continue
            linked = matched.merge(
                indi, left_on=["primaryid", "drug_seq"],
                right_on=["primaryid", "indi_drug_seq"], how="inner",
            )
            keys_by_drug[label].update(eligible_breast_ids(linked))
    return keys_by_drug


def union_sets(keys_by_drug: dict[str, set[str]], labels) -> set[str]:
    result: set[str] = set()
    for label in labels:
        result.update(keys_by_drug[label])
    return result


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
    out["consensus_signal"] = out[["ror_signal", "bcpnn_signal"]].all(axis=1)
    out["four_algorithm_signal"] = out[
        ["ror_signal", "prr_signal", "bcpnn_signal", "mgps_signal"]
    ].all(axis=1)
    return out.sort_values(["consensus_signal", "ror", "a_dato_cases"], ascending=[False, False, False])


def main() -> None:
    background = background_case_universe()
    kept = set(background["primaryid"])
    target, all_dato = target_breast_keys()
    keys_by_drug = comparator_keys_by_drug(kept)
    keys_by_drug = {label: keys - all_dato for label, keys in keys_by_drug.items()}
    trial_aligned = union_sets(keys_by_drug, TRIAL_ALIGNED_DRUGS)
    expanded_chemo = union_sets(keys_by_drug, EXPANDED_CHEMO_DRUGS)
    tdxd_alone = keys_by_drug["trastuzumab_deruxtecan"]
    sg_alone = keys_by_drug["sacituzumab_govitecan"]
    active = union_sets(keys_by_drug, ACTIVE_COMPARATORS)
    class_excluded = union_sets(
        keys_by_drug, [label for label in ACTIVE_COMPARATORS if label not in EXCLUDED_ADC_COMPARATORS]
    )
    sets = {
        "dato_breast_hr_compatible": target,
        "full_faers_excluding_all_dato": kept - all_dato,
        "trial_aligned_chemo": trial_aligned,
        "trastuzumab_deruxtecan_alone": tdxd_alone,
        "sacituzumab_govitecan_alone": sg_alone,
        "expanded_chemo": expanded_chemo,
        "active_breast_comparator": active,
        "active_breast_comparator_class_exclusion": class_excluded,
    }
    counts = reaction_counts(sets)

    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame([{"cohort": name, "n_cases": len(keys)} for name, keys in sets.items()])
    summary.to_csv(TABLES_DIR / "faers_comparator_cohort_sizes.csv", index=False)
    outputs = []
    results_by_name = {}
    for name in list(sets)[1:]:
        result = calculate_signals(counts["dato_breast_hr_compatible"], counts[name], len(target), len(sets[name]), name)
        result.to_csv(TABLES_DIR / f"faers_signals_{name}.csv", index=False)
        outputs.append(result)
        results_by_name[name] = result
    pd.concat(outputs, ignore_index=True).to_csv(TABLES_DIR / "faers_signals_comparator_all.csv", index=False)

    trial_signal_pts = set(
        results_by_name["trial_aligned_chemo"].loc[
            results_by_name["trial_aligned_chemo"]["consensus_signal"], "pt"
        ]
    )
    loo_id_sets = {
        drug_dropped: union_sets(
            keys_by_drug, [drug for drug in TRIAL_ALIGNED_DRUGS if drug != drug_dropped]
        )
        for drug_dropped in TRIAL_ALIGNED_DRUGS
    }
    loo_counts = reaction_counts(loo_id_sets)
    loo_rows = []
    for drug_dropped in TRIAL_ALIGNED_DRUGS:
        comparator_ids = loo_id_sets[drug_dropped]
        loo = calculate_signals(
            counts["dato_breast_hr_compatible"], loo_counts[drug_dropped],
            len(target), len(comparator_ids), f"trial_aligned_drop_{drug_dropped}",
        )
        loo = loo[loo["pt"].isin(trial_signal_pts)].copy()
        loo["drug_dropped"] = drug_dropped
        loo = loo.rename(columns={"n_comparator": "comparator_n"})
        loo_rows.append(loo[[
            "pt", "drug_dropped", "a_dato_cases", "ror", "ror_ci_lower",
            "ror_ci_upper", "ic025", "comparator_n", "consensus_signal",
        ]])
    leave_one_out = pd.concat(loo_rows, ignore_index=True) if loo_rows else pd.DataFrame(columns=[
        "pt", "drug_dropped", "a_dato_cases", "ror", "ror_ci_lower",
        "ror_ci_upper", "ic025", "comparator_n", "consensus_signal",
    ])
    leave_one_out.to_csv(TABLES_DIR / "faers_signals_leave_one_out.csv", index=False)

    print(summary.to_string(index=False))
    for result in outputs:
        name = result["comparator"].iloc[0]
        print(f"\n{name}: consensus signals")
        columns = ["pt", "a_dato_cases", "ror", "ror_ci_lower", "prr", "ic025", "eb05"]
        print(result.loc[result.consensus_signal, columns].head(25).to_string(index=False))


if __name__ == "__main__":
    main()
