"""
Demographic/context tabulation, serious-outcome analysis, subgroup disproportionality,
SOC rollup, and Weibull time-to-onset (TTO) analysis for the datopotamab deruxtecan
FAERS analytic cohort (data/processed/{analysis_cases,case_pt}.csv, n=416 cases).

Five parts, each writing to outputs/tables/:

1. Demographics & context: sex, age band, reporting country, reporter type, report year.
2. Serious-outcome tabulation: FAERS outcome codes, proportion serious, top PTs among
   serious/fatal reports.
3. Subgroup disproportionality: ROR with CI for the 8 consensus signals from
   scripts/05_disproportionality.py, stratified by sex and by age band, against the
   same whole-database background rebuilt here with sex/age attached.
4. SOC aggregation: rolls up signals by MedDRA System Organ Class. No licensed
   MedDRA PT->SOC crosswalk is available (documented gap since scripts/04), so this
   currently produces a single "Unmapped" bucket -- the code is correct and will
   populate properly the moment a crosswalk lands at refs/meddra_pt_soc.csv.
5. Time-to-onset: onset = EVENT_DT (case-level) - START_DT (specific to this drug's
   own drug_seq in THER, not any concomitant drug's), with implausible/negative
   values dropped and logged, medians/IQR overall and per label AESI group, and a
   Weibull fit (lifelines) with shape parameter (rho_) 95% CI and hazard
   classification (rho_<1 early-failure, ~1 random, >1 wear-out).
"""

from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import WeibullFitter

REPO_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
RAW_FAERS_DIR = REPO_ROOT / "data" / "raw" / "faers"
EXTRACT_ROOT = REPO_ROOT / "data" / "raw" / "faers_quarterly" / "extracted"
TABLES_DIR = REPO_ROOT / "outputs" / "tables"

QUARTERS = ["2025q1", "2025q2", "2025q3", "2025q4", "2026q1"]

AGE_BAND_ORDER = ["<18", "18-64", "65-74", ">=75", "Unknown"]

# FAERS age_cod -> multiplier to convert AGE to years.
AGE_COD_TO_YEARS = {"YR": 1, "DEC": 10, "MON": 1 / 12, "WK": 1 / 52.1775, "DY": 1 / 365.25, "HR": 1 / (365.25 * 24)}

OCCP_COD_LABELS = {"MD": "Physician", "PH": "Pharmacist", "HP": "Other health professional",
                    "NR": "Nurse", "LW": "Lawyer", "CN": "Consumer"}

SERIOUS_OUTCOME_CODES = {
    "DE": "Death", "LT": "Life-threatening", "HO": "Hospitalization",
    "DS": "Disability", "CA": "Congenital anomaly", "RI": "Required intervention",
    "OT": "Other serious",
}

# Label-anchored AESI PT groupings. Curated by hand against the actual PT list in
# case_pt.csv, not a blind keyword match -- e.g. infectious pneumonia variants
# ("Pneumonia", "Pneumonia aspiration", "Pneumocystis jirovecii pneumonia") and
# lung-cancer/indication terms are deliberately EXCLUDED from ILD/pneumonitis, since
# they are clinically distinct from the drug-induced ILD/pneumonitis the label
# describes (Section 5.1). Similarly "Dry mouth" (xerostomia, a different mechanism)
# is excluded from Stomatitis, and "Intraocular pressure increased" (a distinct
# glaucoma-type finding, not part of the described corneal/surface toxicity) is
# excluded from Ocular. See docs/candidate_justification.md for the label AESI list.
AESI_GROUPS = {
    "ILD/pneumonitis": ["Interstitial lung disease", "Pneumonitis"],
    "Ocular": [
        "Dry eye", "Keratitis", "Ulcerative keratitis", "Conjunctivitis",
        "Corneal epitheliopathy", "Corneal erosion", "Eye irritation", "Eye pain",
        "Eye pruritus", "Ocular discomfort", "Ocular hyperaemia", "Ocular toxicity",
        "Vision blurred", "Abnormal sensation in eye",
    ],
    "Stomatitis": ["Stomatitis", "Mucosal inflammation"],
    "Infusion-related reaction": [
        "Infusion related reaction", "Infusion related hypersensitivity reaction",
        "Infusion site extravasation", "Infusion site pain", "Infusion site rash",
    ],
}

MAX_PLAUSIBLE_TTO_DAYS = 3650  # 10 years, per protocol


# ============================================================================
# Shared helpers
# ============================================================================

def age_to_years(age, age_cod):
    if pd.isna(age) or pd.isna(age_cod) or age_cod not in AGE_COD_TO_YEARS:
        return np.nan
    return float(age) * AGE_COD_TO_YEARS[age_cod]


def age_band(age_years):
    if pd.isna(age_years):
        return "Unknown"
    if age_years < 18:
        return "<18"
    if age_years < 65:
        return "18-64"
    if age_years < 75:
        return "65-74"
    return ">=75"


def parse_faers_date(s):
    """FAERS dates are YYYYMMDD but can be truncated to YYYYMM or YYYY. Partial
    dates can't give a day-level onset, so they're parsed as NaT (not guessed)."""
    if pd.isna(s):
        return pd.NaT
    s = str(s).strip()
    if len(s) == 8:
        return pd.to_datetime(s, format="%Y%m%d", errors="coerce")
    return pd.NaT


def find_table_file(quarter: str, table: str) -> Path:
    matches = list((EXTRACT_ROOT / quarter / "ASCII").glob(f"{table}*.txt")) + list(
        (EXTRACT_ROOT / quarter / "ASCII").glob(f"{table}*.TXT")
    )
    return matches[0]


def load_deleted_caseids() -> set:
    deleted = set()
    for delete_file in EXTRACT_ROOT.glob("*/Deleted/DELETE*.txt"):
        with open(delete_file) as f:
            for line in f:
                line = line.strip()
                if line:
                    deleted.add(line)
    return deleted


# ============================================================================
# Part 1: demographics & context
# ============================================================================

def part1_demographics(cases: pd.DataFrame):
    print("\n=== PART 1: demographics & context ===")

    sex_tab = cases["sex"].fillna("Unknown").value_counts().rename_axis("sex").reset_index(name="n")
    sex_tab["pct"] = (100 * sex_tab["n"] / len(cases)).round(1)
    sex_tab.to_csv(TABLES_DIR / "demo_sex.csv", index=False)
    print(f"Sex: {dict(zip(sex_tab['sex'], sex_tab['n']))}")

    age_years = cases.apply(lambda r: age_to_years(r["age"], r["age_cod"]), axis=1)
    bands = age_years.apply(age_band)
    age_tab = bands.value_counts().reindex(AGE_BAND_ORDER, fill_value=0).rename_axis("age_band").reset_index(name="n")
    age_tab["pct"] = (100 * age_tab["n"] / len(cases)).round(1)
    age_tab.to_csv(TABLES_DIR / "demo_age_band.csv", index=False)
    print(f"Age bands: {dict(zip(age_tab['age_band'], age_tab['n']))}")

    country_rows = []
    for col, label in [("reporter_country", "reporter"), ("occr_country", "occurrence")]:
        tab = cases[col].fillna("Unknown/missing").value_counts().rename_axis("country").reset_index(name="n")
        tab["country_type"] = label
        tab["pct"] = (100 * tab["n"] / len(cases)).round(1)
        country_rows.append(tab)
    country_tab = pd.concat(country_rows, ignore_index=True)[["country_type", "country", "n", "pct"]]
    country_tab.to_csv(TABLES_DIR / "demo_country.csv", index=False)
    n_occr_missing = (cases["occr_country"].isna()).sum()
    print(f"Country: reporter_country complete; occr_country missing for {n_occr_missing}/{len(cases)} cases")

    occp = cases["occp_cod"].map(OCCP_COD_LABELS).fillna("Unknown/not reported")
    reporter_tab = occp.value_counts().rename_axis("reporter_type").reset_index(name="n")
    reporter_tab["pct"] = (100 * reporter_tab["n"] / len(cases)).round(1)
    reporter_tab.to_csv(TABLES_DIR / "demo_reporter_type.csv", index=False)
    print(f"Reporter type: {dict(zip(reporter_tab['reporter_type'], reporter_tab['n']))}")

    report_year = cases["init_fda_dt"].astype(str).str[:4]
    year_tab = report_year.value_counts().sort_index().rename_axis("report_year").reset_index(name="n")
    year_tab.to_csv(TABLES_DIR / "demo_report_year.csv", index=False)
    print(f"Report year (by init_fda_dt, i.e. first FDA receipt): {dict(zip(year_tab['report_year'], year_tab['n']))}")

    return age_years, bands


# ============================================================================
# Part 2: serious outcomes
# ============================================================================

def part2_serious_outcomes(cases: pd.DataFrame, case_pt: pd.DataFrame):
    print("\n=== PART 2: serious-outcome analysis ===")

    outcome_counts = {code: 0 for code in SERIOUS_OUTCOME_CODES}
    for outcomes_str in cases["outcomes"].fillna(""):
        for code in outcomes_str.split(";"):
            if code in outcome_counts:
                outcome_counts[code] += 1

    n_serious = cases["serious"].sum()
    pct_serious = 100 * n_serious / len(cases)
    print(f"Serious cases: {n_serious}/{len(cases)} ({pct_serious:.1f}%)")

    rows = [{"outc_cod": code, "label": label, "n_cases": outcome_counts[code],
             "pct_of_all_cases": round(100 * outcome_counts[code] / len(cases), 1)}
            for code, label in SERIOUS_OUTCOME_CODES.items()]
    rows.append({"outc_cod": "ANY_SERIOUS", "label": "Any serious outcome",
                 "n_cases": int(n_serious), "pct_of_all_cases": round(pct_serious, 1)})
    outcome_tab = pd.DataFrame(rows)
    outcome_tab.to_csv(TABLES_DIR / "serious_outcomes.csv", index=False)
    for r in rows:
        print(f"  {r['label']} ({r['outc_cod']}): {r['n_cases']} ({r['pct_of_all_cases']}%)")

    fatal_caseids = set(cases.loc[cases["outcomes"].fillna("").str.contains("DE"), "caseid"])
    serious_caseids = set(cases.loc[cases["serious"], "caseid"])

    top_pts_fatal = (
        case_pt[case_pt["caseid"].astype(str).isin(fatal_caseids)]
        .drop_duplicates(subset=["caseid", "pt"])["pt"].value_counts()
        .rename_axis("pt").reset_index(name="n_fatal_cases")
    )
    top_pts_fatal.to_csv(TABLES_DIR / "top_pts_fatal.csv", index=False)
    print(f"\nTop PTs among {len(fatal_caseids)} fatal cases:")
    print(top_pts_fatal.head(10).to_string(index=False))

    top_pts_serious = (
        case_pt[case_pt["caseid"].astype(str).isin(serious_caseids)]
        .drop_duplicates(subset=["caseid", "pt"])["pt"].value_counts()
        .rename_axis("pt").reset_index(name="n_serious_cases")
    )
    top_pts_serious.to_csv(TABLES_DIR / "top_pts_serious.csv", index=False)
    print(f"\nTop PTs among {len(serious_caseids)} serious cases:")
    print(top_pts_serious.head(10).to_string(index=False))


# ============================================================================
# Part 3: subgroup disproportionality
# ============================================================================

def build_background_with_demo():
    """Whole-database background (same FDA-standard dedup as scripts/04 and /05),
    but keeping sex/age so it can be stratified for subgroup ROR."""
    frames = []
    for quarter in QUARTERS:
        demo_file = find_table_file(quarter, "DEMO")
        df = pd.read_csv(
            demo_file, delimiter="$", usecols=["primaryid", "caseid", "fda_dt", "sex", "age", "age_cod"],
            dtype=str, low_memory=False, encoding="latin1",
        )
        df.columns = [c.strip().lower() for c in df.columns]
        frames.append(df)
    demo = pd.concat(frames, ignore_index=True)

    demo["_fda_dt_sort"] = pd.to_numeric(demo["fda_dt"], errors="coerce").fillna(-1)
    demo["_primaryid_sort"] = pd.to_numeric(demo["primaryid"], errors="coerce").fillna(-1)
    demo = demo.sort_values(["_fda_dt_sort", "_primaryid_sort"], ascending=False)
    deduped = demo.drop_duplicates(subset="caseid", keep="first").drop(columns=["_fda_dt_sort", "_primaryid_sort"])

    deleted_caseids = load_deleted_caseids()
    deduped = deduped[~deduped["caseid"].isin(deleted_caseids)].copy()

    deduped["age_years"] = deduped.apply(lambda r: age_to_years(r["age"], r["age_cod"]), axis=1)
    deduped["age_band"] = deduped["age_years"].apply(age_band)
    deduped["sex_group"] = deduped["sex"].fillna("Unknown")
    return deduped


def background_pt_case_ids(kept_primaryids: set, pts: list) -> pd.DataFrame:
    """Returns (primaryid, pt) rows for the given PTs, restricted to kept_primaryids,
    across all quarters -- used to attach sex/age_band to background PT reporters."""
    frames = []
    for quarter in QUARTERS:
        reac_file = find_table_file(quarter, "REAC")
        df = pd.read_csv(reac_file, delimiter="$", usecols=["primaryid", "pt"], dtype=str,
                          low_memory=False, encoding="latin1")
        df.columns = [c.strip().lower() for c in df.columns]
        df = df[df["primaryid"].isin(kept_primaryids) & df["pt"].isin(pts)]
        df = df.drop_duplicates(subset=["primaryid", "pt"])
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def ror_with_ci(a, b, c, d):
    # Haldane-Anscombe continuity correction when any cell is zero -- routine in
    # subgroup analyses since small strata make zero cells common.
    corrected = False
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        corrected = True
    ror = (a * d) / (b * c)
    log_se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lo = np.exp(np.log(ror) - 1.96 * log_se)
    hi = np.exp(np.log(ror) + 1.96 * log_se)
    return ror, lo, hi, corrected


def part3_subgroup_disproportionality(case_pt: pd.DataFrame, n_drug_total: int):
    print("\n=== PART 3: subgroup disproportionality ===")

    signals_path = TABLES_DIR / "signals_significant.csv"
    if not signals_path.exists():
        print(f"  {signals_path} not found -- run scripts/05_disproportionality.py first. Skipping part 3.")
        return
    consensus_pts = pd.read_csv(signals_path)["pt"].tolist()
    print(f"Stratifying the {len(consensus_pts)} consensus-signal PTs by sex and age band.")

    print("Building whole-database background with sex/age (this re-derives the")
    print("scripts/05 background, extended with demographics -- takes ~30s)...")
    background = build_background_with_demo()
    kept_primaryids = set(background["primaryid"])
    bg_pt_rows = background_pt_case_ids(kept_primaryids, consensus_pts)
    bg_pt_rows = bg_pt_rows.merge(background[["primaryid", "sex_group", "age_band"]], on="primaryid", how="left")

    drug_cases = pd.read_csv(PROCESSED_DIR / "analysis_cases.csv", dtype=str)
    drug_cases["age_years"] = drug_cases.apply(lambda r: age_to_years(r["age"], r["age_cod"]), axis=1)
    drug_cases["age_band"] = drug_cases["age_years"].apply(age_band)
    drug_cases["sex_group"] = drug_cases["sex"].fillna("Unknown")
    caseid_to_subgroup = drug_cases.set_index("caseid")[["sex_group", "age_band"]]

    drug_pt = case_pt.drop_duplicates(subset=["caseid", "pt"]).merge(
        caseid_to_subgroup, left_on="caseid", right_index=True, how="left"
    )

    def stratify(dimension: str, groups: list, out_name: str):
        rows = []
        bg_totals = background[dimension].value_counts()
        drug_totals = drug_cases[dimension].value_counts()
        for pt in consensus_pts:
            drug_a_by_group = drug_pt.loc[drug_pt["pt"] == pt, dimension].value_counts()
            bg_pt_by_group = bg_pt_rows.loc[bg_pt_rows["pt"] == pt, dimension].value_counts()
            for group in groups:
                a = int(drug_a_by_group.get(group, 0))
                n_drug_group = int(drug_totals.get(group, 0))
                n_pt_group = int(bg_pt_by_group.get(group, 0))
                n_total_group = int(bg_totals.get(group, 0))
                if n_drug_group == 0 or n_total_group == 0:
                    continue
                b = n_drug_group - a
                c = n_pt_group - a
                d = n_total_group - n_drug_group - n_pt_group + a
                if min(b, c, d) < 0:
                    c = max(c, 0)
                    d = n_total_group - n_drug_group - n_pt_group + a
                    if d < 0:
                        continue
                ror, lo, hi, corrected = ror_with_ci(a, b, c, d)
                rows.append({
                    "pt": pt, dimension: group, "a": a, "n_drug_group": n_drug_group,
                    "n_pt_group": n_pt_group, "n_total_group": n_total_group,
                    "ror": ror, "ror_ci_lower": lo, "ror_ci_upper": hi,
                    "continuity_corrected": corrected,
                })
        out = pd.DataFrame(rows)
        if out.empty:
            print(f"  No subgroup rows computed for {dimension} -- skipping {out_name}.")
            return out

        # Flag PTs whose drug-side cases are concentrated (>=80%) in one subgroup.
        concentration = out.groupby("pt").apply(
            lambda g: g.loc[g["a"].idxmax(), dimension] if g["a"].sum() > 0 and
            g["a"].max() / g["a"].sum() >= 0.8 else None,
            include_groups=False,
        )
        out["concentrated_in"] = out["pt"].map(concentration)
        out.to_csv(TABLES_DIR / out_name, index=False)
        print(f"  Wrote {TABLES_DIR / out_name} ({len(out)} rows)")
        concentrated = out.drop_duplicates("pt").dropna(subset=["concentrated_in"])
        for _, r in concentrated.iterrows():
            note = ""
            if r["concentrated_in"] == "Unknown":
                note = (" -- CAVEAT: this likely reflects missingness (53% of cases have "
                        "unknown age; some have unknown sex), not a true subgroup effect. "
                        "Do not interpret as a real demographic concentration.")
            print(f"    CONCENTRATED: {r['pt']!r} -- >=80% of drug cases in {dimension}={r['concentrated_in']!r}{note}")
        return out

    stratify("sex_group", ["M", "F", "Unknown"], "subgroup_ror_sex.csv")
    stratify("age_band", AGE_BAND_ORDER, "subgroup_ror_age_band.csv")


# ============================================================================
# Part 4: SOC rollup
# ============================================================================

def part4_soc_rollup(case_pt: pd.DataFrame):
    print("\n=== PART 4: SOC aggregation ===")
    signals_path = TABLES_DIR / "signals_all.csv"
    if not signals_path.exists():
        print(f"  {signals_path} not found -- run scripts/05_disproportionality.py first. Skipping part 4.")
        return

    signals = pd.read_csv(signals_path)
    pt_to_soc = case_pt.drop_duplicates("pt").set_index("pt")["soc"]
    signals["soc"] = signals["pt"].map(pt_to_soc)
    n_unmapped = signals["soc"].isna().sum()
    signals["soc"] = signals["soc"].fillna("Unmapped (no licensed MedDRA PT->SOC crosswalk available)")

    if n_unmapped == len(signals):
        print(
            "  GAP: no MedDRA PT->SOC crosswalk available (same gap documented in "
            "scripts/04_clean_dedup.py -- MedDRA's hierarchy is MSSO-licensed). All "
            f"{len(signals)} PTs fall into a single 'Unmapped' bucket below. This code "
            "is correct and will roll up properly the moment a crosswalk is placed at "
            "refs/meddra_pt_soc.csv (scripts/04 will pick it up on re-run)."
        )

    soc_rollup = signals.groupby("soc").agg(
        n_pts=("pt", "nunique"),
        n_consensus_signals=("consensus_signal", "sum"),
        total_cases_a=("a", "sum"),
    ).reset_index().sort_values("n_consensus_signals", ascending=False)
    soc_rollup.to_csv(TABLES_DIR / "soc_rollup.csv", index=False)
    print(soc_rollup.to_string(index=False))


# ============================================================================
# Part 5: time-to-onset
# ============================================================================

def part5_time_to_onset(cases: pd.DataFrame, case_pt: pd.DataFrame):
    print("\n=== PART 5: time-to-onset ===")

    drug = pd.read_csv(RAW_FAERS_DIR / "drug.csv", dtype=str)
    ther = pd.read_csv(RAW_FAERS_DIR / "ther.csv", dtype=str)

    # Isolate THIS drug's own drug_seq per case (drug.csv is already PS-role-filtered
    # to our synonyms from scripts/03), then join THER on primaryid + drug_seq so we
    # get the therapy start date for datopotamab deruxtecan specifically, not a
    # concomitant medication's start date.
    our_drug_seq = drug[["primaryid", "drug_seq"]].drop_duplicates()
    our_ther = ther.merge(
        our_drug_seq, left_on=["primaryid", "dsg_drug_seq"], right_on=["primaryid", "drug_seq"], how="inner"
    )[["primaryid", "start_dt"]].drop_duplicates(subset="primaryid")

    tto = cases[["primaryid", "caseid", "event_dt"]].merge(our_ther, on="primaryid", how="left")
    n_no_ther_row = tto["start_dt"].isna().sum()
    n_start_dt_unparseable = (tto["start_dt"].notna() & tto["start_dt"].astype(str).str.len().ne(8)).sum()
    n_no_event_dt = tto["event_dt"].isna().sum()

    tto["event_date"] = tto["event_dt"].apply(parse_faers_date)
    tto["start_date"] = tto["start_dt"].apply(parse_faers_date)
    tto["onset_days"] = (tto["event_date"] - tto["start_date"]).dt.days

    n_negative = (tto["onset_days"] < 0).sum()
    n_implausible = (tto["onset_days"] > MAX_PLAUSIBLE_TTO_DAYS).sum()
    valid = tto[(tto["onset_days"] >= 0) & (tto["onset_days"] <= MAX_PLAUSIBLE_TTO_DAYS)].copy()

    print(f"Cases with no THER row at all for this drug's drug_seq: {n_no_ther_row}/{len(tto)}")
    print(f"Cases with a THER row but a partial/unparseable START_DT (not full YYYYMMDD): {n_start_dt_unparseable}/{len(tto)}")
    print(f"Cases with no EVENT_DT: {n_no_event_dt}/{len(tto)}")
    print(f"Dropped as negative (event before therapy start): {n_negative}")
    print(f"Dropped as implausible (> {MAX_PLAUSIBLE_TTO_DAYS} days / 10 years): {n_implausible}")
    print(f"Final cases with valid onset_days: {len(valid)}/{len(tto)} "
          f"({100 * len(valid) / len(tto):.1f}%) -- state this completeness rate "
          "explicitly wherever TTO results are reported; it is well under half the cohort.")

    valid.to_csv(TABLES_DIR / "tto_case_level.csv", index=False)

    MIN_N_FOR_CONFIDENT_WEIBULL = 10

    def summarize_group(name: str, days: pd.Series) -> dict:
        n = len(days)
        row = {"group": name, "n": n, "low_n_caveat": n < MIN_N_FOR_CONFIDENT_WEIBULL}
        if n < 2:
            row.update({"median_days": np.nan, "q1_days": np.nan, "q3_days": np.nan,
                        "weibull_shape": np.nan, "shape_ci_lower": np.nan,
                        "shape_ci_upper": np.nan, "hazard_pattern": "insufficient data (n<2)"})
            return row
        row["median_days"] = days.median()
        row["q1_days"] = days.quantile(0.25)
        row["q3_days"] = days.quantile(0.75)

        durations = days.clip(lower=0.5)  # Weibull requires positive support
        wf = WeibullFitter()
        try:
            wf.fit(durations)
            shape = wf.rho_
            ci_lo = wf.summary.loc["rho_", "coef lower 95%"]
            ci_hi = wf.summary.loc["rho_", "coef upper 95%"]
            if ci_hi < 1:
                pattern = "early-failure (decreasing hazard)"
            elif ci_lo > 1:
                pattern = "wear-out (increasing hazard)"
            else:
                pattern = "consistent with random/constant hazard (CI spans 1)"
            row.update({"weibull_shape": shape, "shape_ci_lower": ci_lo,
                        "shape_ci_upper": ci_hi, "hazard_pattern": pattern})
        except Exception as e:
            row.update({"weibull_shape": np.nan, "shape_ci_lower": np.nan,
                        "shape_ci_upper": np.nan, "hazard_pattern": f"fit failed: {e}"})
        return row

    summary_rows = [summarize_group("Overall", valid["onset_days"])]

    for aesi_name, pts in AESI_GROUPS.items():
        aesi_caseids = set(case_pt.loc[case_pt["pt"].isin(pts), "caseid"].astype(str))
        group_days = valid.loc[valid["caseid"].astype(str).isin(aesi_caseids), "onset_days"]
        summary_rows.append(summarize_group(aesi_name, group_days))

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TABLES_DIR / "tto_summary.csv", index=False)
    print("\nTTO summary:")
    print(summary.to_string(index=False))
    low_n_groups = summary.loc[summary["low_n_caveat"] & (summary["n"] >= 2), "group"].tolist()
    if low_n_groups:
        print(
            f"\nCAVEAT: {low_n_groups} have n < {MIN_N_FOR_CONFIDENT_WEIBULL} valid-onset "
            "cases. Their Weibull shape/CI and hazard classification are reported for "
            "completeness but are not reliable at this sample size -- state this plainly "
            "in the manuscript rather than presenting these hazard classifications with "
            "the same confidence as the overall estimate."
        )


def main():
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    cases = pd.read_csv(PROCESSED_DIR / "analysis_cases.csv", dtype=str)
    cases["serious"] = cases["serious"].map({"True": True, "False": False}).fillna(False)
    case_pt = pd.read_csv(PROCESSED_DIR / "case_pt.csv", dtype=str)

    part1_demographics(cases)
    part2_serious_outcomes(cases, case_pt)
    part3_subgroup_disproportionality(case_pt, n_drug_total=cases["caseid"].nunique())
    part4_soc_rollup(case_pt)
    part5_time_to_onset(cases, case_pt)

    print("\nAll outputs written to outputs/tables/.")


if __name__ == "__main__":
    main()
