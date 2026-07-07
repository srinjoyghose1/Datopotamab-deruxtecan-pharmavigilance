"""
Build the analysis dataset from the raw FAERS pull in data/raw/faers/.

Pipeline (see docs/protocol.md for the governing protocol; as of this script's
writing that file does not yet exist in the repo, so the steps below follow the
parameters given directly in the task that requested this script):

1. FDA-standard deduplication: within duplicate CASEIDs, keep the record with the
   latest FDA_DT, breaking remaining ties by the highest PRIMARYID. Cases present
   on any FDA quarterly deleted-cases list (data/raw/faers_quarterly/extracted/*/
   Deleted/DELETE*.txt) are removed entirely.

2. Case selection (Primary Suspect): the raw pull (scripts/03_pull_faers.py) already
   restricts to drug records where datopotamab deruxtecan is role_cod == "PS" before
   building the case universe -- every primaryid in data/raw/faers/demo.csv is
   therefore guaranteed to have a PS-role match. This step is consequently a
   verification, not a re-filter: it asserts every drug record backing the final
   case set is PS and logs the count. This ordering (PS restriction happens at
   acquisition, before dedup, rather than after) is a deliberate simplification and
   is stated here explicitly rather than silently assumed -- see the note in
   data/processed/README.md.

3. Event coding: reaction terms come from REAC.pt, which FAERS already codes as
   MedDRA Preferred Terms. This script does not itself determine the MedDRA
   version -- the ASCII quarterly extract does not carry a per-record version tag
   (only the XML extract's <reactionmeddraversionpt> tag does; see
   docs/environment.md for the version recovered from a sample XML pull).
   PT -> SOC mapping is attempted only if a local crosswalk file is present (see
   MEDDRA_PT_SOC_PATH below); MedDRA's hierarchy is MSSO-licensed and not freely
   redistributable, so if that file is absent, the gap is logged explicitly and
   the soc column is left null rather than guessed.

Outputs:
  data/processed/analysis_cases.csv  -- one row per unique (deduped, PS, non-deleted) case
  data/processed/case_pt.csv         -- one row per case-PT pair
"""

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = REPO_ROOT / "data" / "raw" / "faers"
EXTRACT_ROOT = REPO_ROOT / "data" / "raw" / "faers_quarterly" / "extracted"
OUT_DIR = REPO_ROOT / "data" / "processed"

# Optional local MedDRA PT->SOC crosswalk. Not shipped with this repo (MedDRA's
# hierarchy is MSSO-licensed); if you have access to one, drop it here as a CSV
# with at minimum "pt" and "soc" columns.
MEDDRA_PT_SOC_PATH = REPO_ROOT / "refs" / "meddra_pt_soc.csv"


def load_deleted_caseids() -> set:
    deleted = set()
    for delete_file in EXTRACT_ROOT.glob("*/Deleted/DELETE*.txt"):
        with open(delete_file) as f:
            for line in f:
                line = line.strip()
                if line:
                    deleted.add(line)
    return deleted


def dedup_cases(demo: pd.DataFrame, drug: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    log = {"raw_report_versions": len(demo)}

    # Step 2 verification: every primaryid in demo must have a PS-role match in
    # drug (guaranteed by construction in scripts/03_pull_faers.py, checked here
    # rather than assumed).
    ps_primaryids = set(drug.loc[drug["role_cod"] == "PS", "primaryid"])
    non_ps = set(demo["primaryid"]) - ps_primaryids
    assert not non_ps, (
        f"{len(non_ps)} primaryid(s) in demo.csv lack a PS-role drug record -- "
        "acquisition-time PS filtering assumption violated, investigate before proceeding."
    )
    log["ps_verified_report_versions"] = len(demo)

    # Step 1: FDA-standard dedup. fda_dt is YYYYMMDD; higher = more recent.
    # primaryid is numeric; higher = later version, used as tiebreak.
    demo = demo.copy()
    demo["_fda_dt_sort"] = pd.to_numeric(demo["fda_dt"], errors="coerce").fillna(-1)
    demo["_primaryid_sort"] = pd.to_numeric(demo["primaryid"], errors="coerce").fillna(-1)
    demo = demo.sort_values(["_fda_dt_sort", "_primaryid_sort"], ascending=False)
    deduped = demo.drop_duplicates(subset="caseid", keep="first").drop(
        columns=["_fda_dt_sort", "_primaryid_sort"]
    )
    log["deduped_unique_cases"] = len(deduped)

    # Remove FDA-deleted cases.
    deleted_caseids = load_deleted_caseids()
    before = len(deduped)
    deduped = deduped[~deduped["caseid"].isin(deleted_caseids)]
    log["removed_as_fda_deleted"] = before - len(deduped)
    log["final_cases"] = len(deduped)

    return deduped, log


def build_analysis_cases(deduped_demo: pd.DataFrame, outc: pd.DataFrame, reac: pd.DataFrame) -> pd.DataFrame:
    kept_primaryids = set(deduped_demo["primaryid"])

    outc_subset = outc[outc["primaryid"].isin(kept_primaryids)]
    outcomes_by_case = (
        outc_subset.groupby("primaryid")["outc_cod"]
        .apply(lambda s: ";".join(sorted(s.dropna().unique())))
        .rename("outcomes")
    )
    serious_by_case = outc_subset.groupby("primaryid").size().rename("n_outcomes") > 0

    reac_subset = reac[reac["primaryid"].isin(kept_primaryids)]
    n_reactions_by_case = reac_subset.groupby("primaryid").size().rename("n_reactions")

    cases = deduped_demo[[
        "primaryid", "caseid", "quarter", "sex", "age", "age_cod", "age_grp",
        "wt", "wt_cod", "event_dt", "init_fda_dt", "fda_dt", "rept_dt",
        "rept_cod", "occp_cod", "reporter_country", "occr_country",
    ]].copy()

    cases = cases.merge(outcomes_by_case, on="primaryid", how="left")
    cases["outcomes"] = cases["outcomes"].fillna("")
    cases = cases.merge(serious_by_case.rename("serious"), on="primaryid", how="left")
    cases["serious"] = cases["serious"].astype("boolean").fillna(False)
    cases = cases.merge(n_reactions_by_case, on="primaryid", how="left")
    cases["n_reactions"] = cases["n_reactions"].fillna(0).astype(int)

    return cases


def build_case_pt(deduped_demo: pd.DataFrame, reac: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    kept_primaryids = set(deduped_demo["primaryid"])
    caseid_by_primaryid = deduped_demo.set_index("primaryid")["caseid"]

    case_pt = reac[reac["primaryid"].isin(kept_primaryids)][["primaryid", "pt", "quarter"]].copy()
    case_pt["caseid"] = case_pt["primaryid"].map(caseid_by_primaryid)

    log = {"case_pt_rows": len(case_pt), "unique_pts": case_pt["pt"].nunique()}

    if MEDDRA_PT_SOC_PATH.exists():
        soc_map = pd.read_csv(MEDDRA_PT_SOC_PATH, dtype=str)
        soc_map.columns = [c.strip().lower() for c in soc_map.columns]
        case_pt = case_pt.merge(soc_map[["pt", "soc"]], on="pt", how="left")
        n_unmapped = case_pt["soc"].isna().sum()
        log["soc_mapping"] = "applied"
        log["pt_rows_unmapped_to_soc"] = int(n_unmapped)
    else:
        case_pt["soc"] = pd.NA
        log["soc_mapping"] = "SKIPPED -- no MedDRA PT->SOC crosswalk available"
        print(
            f"\nGAP: no MedDRA PT->SOC crosswalk found at {MEDDRA_PT_SOC_PATH}. "
            "MedDRA's hierarchy is MSSO-licensed and not freely redistributable, "
            "so this project cannot map PT to SOC without a licensed crosswalk. "
            "The 'soc' column in case_pt.csv is left null. State this gap in the "
            "manuscript methods if PT-level results are reported without SOC grouping."
        )

    return case_pt, log


def main():
    demo = pd.read_csv(RAW_DIR / "demo.csv", dtype=str)
    drug = pd.read_csv(RAW_DIR / "drug.csv", dtype=str)
    outc = pd.read_csv(RAW_DIR / "outc.csv", dtype=str)
    reac = pd.read_csv(RAW_DIR / "reac.csv", dtype=str)

    deduped_demo, dedup_log = dedup_cases(demo, drug)
    analysis_cases = build_analysis_cases(deduped_demo, outc, reac)
    case_pt, pt_log = build_case_pt(deduped_demo, reac)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    analysis_cases.to_csv(OUT_DIR / "analysis_cases.csv", index=False)
    case_pt.to_csv(OUT_DIR / "case_pt.csv", index=False)

    print("=== DATA PROVENANCE SUMMARY ===")
    print(f"Raw PS-matched report versions (data/raw/faers/demo.csv): {dedup_log['raw_report_versions']}")
    print(f"  -> verified all backed by a PS-role drug record: {dedup_log['ps_verified_report_versions']}")
    print(f"  -> deduped to unique cases (latest FDA_DT, tie-break highest PRIMARYID): {dedup_log['deduped_unique_cases']}")
    print(f"  -> removed as present on an FDA quarterly deleted-cases list: {dedup_log['removed_as_fda_deleted']}")
    print(f"  -> final analytic cases: {dedup_log['final_cases']}")
    print(f"  -> case-PT rows (data/processed/case_pt.csv): {pt_log['case_pt_rows']} ({pt_log['unique_pts']} unique PTs)")
    print(f"  -> SOC mapping: {pt_log['soc_mapping']}")
    print(f"\nWrote {OUT_DIR / 'analysis_cases.csv'} ({len(analysis_cases)} rows)")
    print(f"Wrote {OUT_DIR / 'case_pt.csv'} ({len(case_pt)} rows)")

    return {**dedup_log, **pt_log}


if __name__ == "__main__":
    main()
