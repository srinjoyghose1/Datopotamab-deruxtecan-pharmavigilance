"""JADER comparator sensitivity analyses for datopotamab deruxtecan.

Uses PMDA's full JADER CSV release (CP932 encoding) and treats one
``identifier + report_count`` pair as one report version. JADER does not expose
a cross-reporter master case identifier, so duplicate reports can remain.

The target cohort is datopotamab deruxtecan as a suspected drug with an explicit
breast-cancer reason for use. Three case-level comparators are evaluated:

1. all JADER report versions excluding target-cohort reports;
2. breast-cancer reports for the prespecified active comparators;
3. (2), excluding any report that also names sacituzumab govitecan or
   trastuzumab deruxtecan as a suspected drug.

Run, for example:
  python scripts/07_jader_comparator.py \
      --jader-dir /Users/srinjoy/Downloads/pmdacasereport202606
"""

from __future__ import annotations

import argparse
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma


REPO_ROOT = Path(__file__).resolve().parents[1]
TABLES_DIR = REPO_ROOT / "outputs" / "tables"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
ENCODING = "cp932"
CHUNK_SIZE = 200_000

# Curated display translations for PTs observed in the Dato-DXd cohort. The
# source Japanese term is always retained because this is not a redistributable
# MedDRA bilingual dictionary.
PT_ENGLISH = {
    "間質性肺疾患": "Interstitial lung disease",
    "感染性虹彩毛様体炎": "Infectious iridocyclitis",
    "角膜上皮症": "Corneal epitheliopathy",
    "アフタ性潰瘍": "Aphthous ulcer",
    "咽喉頭炎": "Pharyngolaryngitis",
    "角膜炎": "Keratitis",
    "疾患進行": "Disease progression",
    "眼圧上昇": "Intraocular pressure increased",
    "口内炎": "Stomatitis",
    "炎症": "Inflammation",
    "注入に伴う反応": "Infusion related reaction",
    "出血性膀胱炎": "Haemorrhagic cystitis",
    "骨髄炎": "Osteomyelitis",
    "背部痛": "Back pain",
    "圧迫骨折": "Compression fracture",
    "肺毒性": "Pulmonary toxicity",
    "うっ血性心不全": "Cardiac failure congestive",
    "脳出血": "Cerebral haemorrhage",
    "嘔吐": "Vomiting",
    "悪心": "Nausea",
    "顎骨壊死": "Osteonecrosis of jaw",
    "多形紅斑": "Erythema multiforme",
    "骨髄抑制": "Bone marrow suppression",
    "呼吸困難": "Dyspnoea",
    "倦怠感": "Malaise",
    "好中球数減少": "Neutrophil count decreased",
    "急性腎障害": "Acute kidney injury",
    "下痢": "Diarrhoea",
    "白血球数減少": "White blood cell count decreased",
    "アナフィラキシー反応": "Anaphylactic reaction",
    "貧血": "Anaemia",
    "肝機能異常": "Hepatic function abnormal",
    "発熱": "Pyrexia",
}
OUTCOME_ENGLISH = {
    "回復": "Recovered/resolved",
    "軽快": "Recovering/resolving",
    "未回復": "Not recovered/not resolved",
    "後遺症あり": "Recovered/resolved with sequelae",
    "死亡": "Fatal",
    "不明": "Unknown",
}

# Generic-name substrings used by the PMDA June 2026 export. Paclitaxel captures
# both conventional and albumin-bound paclitaxel records, as intended.
ACTIVE_COMPARATORS = {
    "capecitabine": "カペシタビン",
    "eribulin": "エリブリン",
    "gemcitabine": "ゲムシタビン",
    "vinorelbine": "ビノレルビン",
    "paclitaxel_or_nab_paclitaxel": "パクリタキセル",
    "carboplatin": "カルボプラチン",
    "sacituzumab_govitecan": "サシツズマブ",
    "trastuzumab_deruxtecan": "トラスツズマブ　デルクステカン",
}
EXCLUDED_ADC_COMPARATORS = {
    "sacituzumab_govitecan",
    "trastuzumab_deruxtecan",
}


def key_frame(df: pd.DataFrame) -> pd.Series:
    return df["識別番号"].astype(str) + "|" + df["報告回数"].astype(str)


def ror_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    """ROR and log-method 95% CI; apply Haldane correction only if needed."""
    if min(a, b, c, d) == 0:
        a, b, c, d = (value + 0.5 for value in (a, b, c, d))
    ror = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return ror, math.exp(math.log(ror) - 1.96 * se), math.exp(math.log(ror) + 1.96 * se)


def prr_chi2(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    raw_table = [[a, b], [c, d]]
    if min(a, b, c, d) == 0:
        a, b, c, d = (value + 0.5 for value in (a, b, c, d))
    prr = (a / (a + b)) / (c / (c + d))
    # The corrected table permits a finite PRR. Chi-square is undefined for a
    # zero-margin raw table, so use the same Haldane-corrected table there.
    chi_table = raw_table if min(value for row in raw_table for value in row) > 0 else [[a, b], [c, d]]
    chi2, _, _, _ = stats.chi2_contingency(chi_table, correction=True)
    return prr, chi2


def bcpnn_ic(a: int, n_drug: int, n_pt: int, n_total: int) -> tuple[float, float, float]:
    expected = (n_drug * n_pt) / n_total
    ic = np.log2((a + 0.5) / (expected + 0.5))
    variance = (1 / np.log(2) ** 2) * (1 / (a + 0.5) + 1 / (expected + 0.5))
    return ic, ic - 1.96 * np.sqrt(variance), expected


def fit_gamma_poisson(a: np.ndarray, expected: np.ndarray) -> tuple[float, float]:
    ratios = a / expected
    mean = ratios.mean()
    variance = ratios.var(ddof=1)
    if variance <= 0:
        return 1.0, 1.0
    beta = mean / variance
    return mean * beta, beta


def ebgm_eb05(a: int, expected: float, alpha: float, beta: float) -> tuple[float, float]:
    post_alpha, post_beta = alpha + a, beta + expected
    return (
        float(np.exp(digamma(post_alpha) - np.log(post_beta))),
        float(stats.gamma.ppf(0.05, a=post_alpha, scale=1 / post_beta)),
    )


def load_drug_sets(drug_path: Path) -> tuple[set[str], set[str], set[str], pd.DataFrame]:
    """Build target and active-comparator report-version sets from DRUG."""
    target_keys: set[str] = set()
    active_keys: set[str] = set()
    excluded_adc_keys: set[str] = set()
    target_rows: list[pd.DataFrame] = []
    columns = ["識別番号", "報告回数", "医薬品の関与", "医薬品（一般名）", "使用理由"]

    for chunk in pd.read_csv(drug_path, encoding=ENCODING, dtype=str, usecols=columns, chunksize=CHUNK_SIZE):
        suspect = chunk[chunk["医薬品の関与"].eq("被疑薬")].copy()
        suspect["key"] = key_frame(suspect)
        name = suspect["医薬品（一般名）"].fillna("")
        reason = suspect["使用理由"].fillna("")
        breast = reason.str.contains("乳癌", regex=False)

        dato = name.str.contains("ダトポタマブ", regex=False)
        target = suspect[dato & breast]
        target_keys.update(target["key"])
        target_rows.append(target)

        for label, term in ACTIVE_COMPARATORS.items():
            matched = name.str.contains(term, regex=False)
            matched_keys = set(suspect.loc[matched & breast, "key"])
            active_keys.update(matched_keys)
            if label in EXCLUDED_ADC_COMPARATORS:
                excluded_adc_keys.update(set(suspect.loc[matched, "key"]))

    target_df = pd.concat(target_rows, ignore_index=True)
    if len(target_keys) != len(target_df):
        raise ValueError("Target drug has multiple suspected-drug rows per report version; inspect before proceeding.")
    if target_df.empty:
        raise ValueError("No Dato-DXd breast-cancer suspected-drug records found.")
    return target_keys, active_keys, excluded_adc_keys, target_df


def count_reactions(reac_path: Path, sets: dict[str, set[str]]) -> dict[str, Counter]:
    counters = {name: Counter() for name in sets}
    for chunk in pd.read_csv(reac_path, encoding=ENCODING, dtype=str, chunksize=CHUNK_SIZE):
        chunk["key"] = key_frame(chunk)
        unique_pt = chunk.drop_duplicates(["key", "有害事象"])
        for name, case_keys in sets.items():
            matched = unique_pt[unique_pt["key"].isin(case_keys)]
            counters[name].update(matched["有害事象"].dropna())
    return counters


def write_normalized_target_files(
    demo_path: Path, reac_path: Path, target_rows: pd.DataFrame, target_keys: set[str], release: str
) -> None:
    """Write small UTF-8 files with English headers for the Dato-DXd cohort."""
    demo = pd.read_csv(demo_path, encoding=ENCODING, dtype=str)
    demo["case_key"] = key_frame(demo)
    demo = demo[demo["case_key"].isin(target_keys)].copy()
    demo = demo.rename(columns={
        "識別番号": "identifier", "報告回数": "report_count", "性別": "sex_japanese",
        "年齢": "age_japanese", "体重": "weight_japanese", "身長": "height_japanese",
        "報告年度・四半期": "report_quarter_japanese", "状況": "status_japanese",
        "報告の種類": "report_type_japanese", "報告者の資格": "reporter_qualification_japanese",
        "E2B": "e2b_version",
    })
    quarter_text = demo["report_quarter_japanese"].fillna("")
    demo["report_year"] = quarter_text.str.extract(r"(\d{4})", expand=False)
    quarter_number = quarter_text.str.extract(r"([1-4])", expand=False)
    for japanese, number in {"第一": "1", "第二": "2", "第三": "3", "第四": "4"}.items():
        quarter_number = quarter_number.mask(quarter_text.str.contains(japanese, regex=False), number)
    demo["report_quarter"] = quarter_number.map(lambda value: f"Q{value}" if pd.notna(value) else pd.NA)
    demo["jader_release"] = release

    target = target_rows.rename(columns={
        "識別番号": "identifier", "報告回数": "report_count",
        "医薬品の関与": "drug_role_japanese", "医薬品（一般名）": "generic_name_japanese",
        "使用理由": "indication_japanese",
    }).drop(columns="key")
    cases = demo.merge(target, on=["identifier", "report_count"], how="left")
    ordered = [
        "case_key", "identifier", "report_count", "report_year", "report_quarter",
        "sex_japanese", "age_japanese", "weight_japanese", "height_japanese",
        "status_japanese", "report_type_japanese", "reporter_qualification_japanese",
        "e2b_version", "drug_role_japanese", "generic_name_japanese",
        "indication_japanese", "report_quarter_japanese", "jader_release",
    ]

    reactions = pd.read_csv(reac_path, encoding=ENCODING, dtype=str)
    reactions["case_key"] = key_frame(reactions)
    reactions = reactions[reactions["case_key"].isin(target_keys)].rename(columns={
        "識別番号": "identifier", "報告回数": "report_count",
        "有害事象連番": "reaction_sequence", "有害事象": "pt_japanese",
        "転帰": "outcome_japanese", "有害事象の発現日": "event_onset_date",
    })
    reactions["pt_english"] = reactions["pt_japanese"].map(PT_ENGLISH)
    reactions["outcome_english"] = reactions["outcome_japanese"].map(OUTCOME_ENGLISH)
    reactions["jader_release"] = release
    reaction_order = [
        "case_key", "identifier", "report_count", "reaction_sequence", "pt_japanese",
        "pt_english", "outcome_japanese", "outcome_english", "event_onset_date", "jader_release",
    ]

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    cases[ordered].sort_values("case_key").to_csv(PROCESSED_DIR / "jader_dato_cases.csv", index=False)
    reactions[reaction_order].sort_values(["case_key", "reaction_sequence"]).to_csv(
        PROCESSED_DIR / "jader_dato_case_pt.csv", index=False
    )


def calculate_signals(
    target_counts: Counter, comparator_counts: Counter, n_target: int, n_comparator: int, comparator: str
) -> pd.DataFrame:
    rows = []
    for pt, a in target_counts.items():
        b = n_target - a
        c = comparator_counts.get(pt, 0)
        d = n_comparator - c
        # c/d describe the comparator only; target records are already excluded.
        ror, lo, hi = ror_ci(a, b, c, d)
        prr, chi2 = prr_chi2(a, b, c, d)
        ic, ic025, expected = bcpnn_ic(a, n_target, a + c, n_target + n_comparator)
        rows.append({
            "pt_japanese": pt,
            "comparator": comparator,
            "a_dato_cases": a,
            "b_dato_without_pt": b,
            "c_comparator_cases": c,
            "d_comparator_without_pt": d,
            "n_dato": n_target,
            "n_comparator": n_comparator,
            "ror": ror,
            "ror_ci_lower": lo,
            "ror_ci_upper": hi,
            "prr": prr,
            "chi2_yates": chi2,
            "ic": ic,
            "ic025": ic025,
            "expected_count": expected,
        })

    out = pd.DataFrame(rows)
    out.insert(1, "pt_english", out["pt_japanese"].map(PT_ENGLISH))
    alpha, beta = fit_gamma_poisson(out["a_dato_cases"].to_numpy(), out["expected_count"].to_numpy())
    eb = [ebgm_eb05(a, e, alpha, beta) for a, e in zip(out["a_dato_cases"], out["expected_count"])]
    out["ebgm"], out["eb05"] = zip(*eb)
    out["ror_signal"] = (out["a_dato_cases"] >= 3) & (out["ror_ci_lower"] > 1)
    out["prr_signal"] = (out["a_dato_cases"] >= 3) & (out["prr"] >= 2) & (out["chi2_yates"] >= 4)
    out["bcpnn_signal"] = (out["a_dato_cases"] >= 3) & (out["ic025"] > 0)
    out["mgps_signal"] = (out["a_dato_cases"] >= 3) & (out["eb05"] >= 2)
    out["consensus_signal"] = out[["ror_signal", "prr_signal", "bcpnn_signal"]].all(axis=1)
    return out.sort_values(["consensus_signal", "ror", "a_dato_cases"], ascending=[False, False, False])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jader-dir", type=Path, required=True, help="Directory holding demo/drug/reac YYYYMM.csv files")
    args = parser.parse_args()
    jader_dir = args.jader_dir
    drug_path = next(jader_dir.glob("drug*.csv"))
    demo_path = next(jader_dir.glob("demo*.csv"))
    reac_path = next(jader_dir.glob("reac*.csv"))

    target_keys, active_keys, excluded_adc_keys, target_rows = load_drug_sets(drug_path)
    demo = pd.read_csv(demo_path, encoding=ENCODING, dtype=str, usecols=["識別番号", "報告回数", "報告年度・四半期"])
    demo["key"] = key_frame(demo)
    all_keys = set(demo["key"])
    write_normalized_target_files(demo_path, reac_path, target_rows, target_keys, jader_dir.name)

    active_keys -= target_keys
    active_class_excluded_keys = active_keys - excluded_adc_keys
    comparators = {
        "full_jader_excluding_dato": all_keys - target_keys,
        "active_breast_comparator": active_keys,
        "active_breast_comparator_class_exclusion": active_class_excluded_keys,
    }
    counters = count_reactions(reac_path, {"dato": target_keys, **comparators})

    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame([
        {"cohort": "dato_breast", "n_report_versions": len(target_keys)},
        *[{"cohort": name, "n_report_versions": len(keys)} for name, keys in comparators.items()],
    ])
    summary["jader_release"] = jader_dir.name
    summary.to_csv(TABLES_DIR / "jader_comparator_cohort_sizes.csv", index=False)

    signals = []
    for name, keys in comparators.items():
        frame = calculate_signals(counters["dato"], counters[name], len(target_keys), len(keys), name)
        frame.to_csv(TABLES_DIR / f"jader_signals_{name}.csv", index=False)
        signals.append(frame)
    pd.concat(signals, ignore_index=True).to_csv(TABLES_DIR / "jader_signals_comparator_all.csv", index=False)

    print(summary.to_string(index=False))
    for name in comparators:
        top = pd.read_csv(TABLES_DIR / f"jader_signals_{name}.csv")
        print(f"\n{name}: consensus signals")
        print(top.loc[top["consensus_signal"], ["pt_japanese", "a_dato_cases", "ror", "ror_ci_lower", "prr", "ic025", "eb05"]].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
