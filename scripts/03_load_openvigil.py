"""
Path A (matches the reference papers' tooling): load and standardize an OpenVigil
FDA / OpenVigil 2.1 export for datopotamab deruxtecan.

This script does NOT query OpenVigil itself -- OpenVigil is a web tool with no
scriptable API for bulk export, so the query has to be run by hand:

    1. Go to the OpenVigil FDA / OpenVigil 2.1 web interface.
    2. Query drug name using the synonym list below (OR'd together, or run once
       per synonym if the tool only accepts one term at a time).
    3. Set Role = Primary Suspect (PS).
    4. Set the report window to 2025Q1 through the latest available quarter.
    5. Export BOTH the "Raw data" and "Frequency" outputs as CSV.
    6. Save them into data/raw/openvigil/ as:
         data/raw/openvigil/raw_data.csv
         data/raw/openvigil/frequency.csv

Run this script after that manual step. It loads whatever is in data/raw/openvigil/,
runs the same drug-name synonym verification as the Path B (programmatic) script so
the two paths are directly comparable, standardizes column names onto a common
schema, and writes standardized copies alongside the originals.

If data/raw/openvigil/ is empty, this script prints instructions and exits without
error -- this is the expected state when Path B is being used instead.
"""

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
OV_DIR = REPO_ROOT / "data" / "raw" / "openvigil"

SYNONYMS = [
    "datopotamab deruxtecan",
    "datopotamab deruxtecan-dlnk",
    "datroway",
    "dato-dxd",
    "ds-1062",
    "ds-1062a",
]

RAW_DATA_FILE = OV_DIR / "raw_data.csv"
FREQUENCY_FILE = OV_DIR / "frequency.csv"

# Best-effort column-name mapping for OpenVigil FDA / 2.1 exports onto a common
# schema shared with the Path B (FAERS quarterly ASCII) output. OpenVigil's exact
# column headers vary by version/config, so this mapping is deliberately loose
# (case/underscore/space-insensitive) and the script prints anything it could not
# map rather than silently dropping or mis-mapping columns.
RAW_DATA_COLUMN_ALIASES = {
    "case_id": ["caseid", "case_id", "primaryid", "case"],
    "drug_name": ["drugname", "drug_name", "medicinalproduct", "substance"],
    "role_cod": ["role_cod", "role", "drug_role"],
    "event_pt": ["pt", "event", "reaction", "adr", "meddra_pt"],
    "sex": ["sex", "gender"],
    "age": ["age"],
    "age_unit": ["age_cod", "age_unit"],
    "country": ["country", "reporter_country", "occr_country"],
    "event_date": ["event_dt", "event_date", "onset_date"],
    "report_date": ["fda_dt", "report_date", "receive_date"],
    "outcome": ["outc_cod", "outcome"],
    "serious": ["serious", "seriousness"],
}

FREQUENCY_COLUMN_ALIASES = {
    "event_pt": ["pt", "event", "reaction", "adr", "meddra_pt"],
    "count_drug": ["count_drug", "n_drug", "drug_count", "count(drug)"],
    "count_background": ["count_background", "n_background", "background_count"],
    "total_drug": ["total_drug", "n_total_drug"],
    "total_background": ["total_background", "n_total_background"],
    "ror": ["ror", "reporting_odds_ratio"],
    "ror_lower": ["ror_lower", "ror_ci_lower", "ror025"],
    "ror_upper": ["ror_upper", "ror_ci_upper", "ror975"],
    "prr": ["prr", "proportional_reporting_ratio"],
}


def _normalize(name: str) -> str:
    return name.strip().lower().replace(" ", "_").replace("-", "_")


def standardize_columns(df: pd.DataFrame, aliases: dict) -> pd.DataFrame:
    normalized_to_original = {_normalize(c): c for c in df.columns}
    rename_map = {}
    for standard_name, alias_list in aliases.items():
        for alias in alias_list:
            key = _normalize(alias)
            if key in normalized_to_original:
                rename_map[normalized_to_original[key]] = standard_name
                break
    df = df.rename(columns=rename_map)
    unmapped = [c for c in df.columns if c not in aliases]
    if unmapped:
        print(f"  NOTE: columns left as-is (not in standard schema): {unmapped}")
    return df


def synonym_match_counts(raw_df: pd.DataFrame) -> dict:
    if "drug_name" not in raw_df.columns:
        print("  WARNING: no drug_name column found after standardization; "
              "skipping per-synonym verification.")
        return {}
    name_col = raw_df["drug_name"].fillna("").str.lower()
    counts = {}
    for syn in SYNONYMS:
        n = int(name_col.str.contains(syn.lower(), regex=False).sum())
        counts[syn] = n
    return counts


def main():
    if not OV_DIR.exists() or not (RAW_DATA_FILE.exists() or FREQUENCY_FILE.exists()):
        print(f"No OpenVigil export found in {OV_DIR}.")
        print("This is expected if you're using Path B (programmatic FAERS pull) instead.")
        print("To use Path A, export from OpenVigil FDA / OpenVigil 2.1 and save as:")
        print(f"  {RAW_DATA_FILE}")
        print(f"  {FREQUENCY_FILE}")
        return None

    result = {"source": "openvigil_export"}

    if RAW_DATA_FILE.exists():
        print(f"Loading {RAW_DATA_FILE}...")
        raw_df = pd.read_csv(RAW_DATA_FILE, dtype=str, low_memory=False)
        raw_df = standardize_columns(raw_df, RAW_DATA_COLUMN_ALIASES)

        counts = synonym_match_counts(raw_df)
        print("\nSynonym match counts (OpenVigil raw_data export):")
        for syn, n in counts.items():
            print(f"  {syn!r}: {n} records")
        result["synonym_counts"] = counts
        result["n_raw_records"] = len(raw_df)

        out_path = OV_DIR / "raw_data_standardized.csv"
        raw_df.to_csv(out_path, index=False)
        print(f"Wrote standardized raw data: {out_path} ({len(raw_df)} rows)")
    else:
        print(f"NOTE: {RAW_DATA_FILE} not found, skipping raw-data standardization.")

    if FREQUENCY_FILE.exists():
        print(f"\nLoading {FREQUENCY_FILE}...")
        freq_df = pd.read_csv(FREQUENCY_FILE, dtype=str, low_memory=False)
        freq_df = standardize_columns(freq_df, FREQUENCY_COLUMN_ALIASES)
        out_path = OV_DIR / "frequency_standardized.csv"
        freq_df.to_csv(out_path, index=False)
        print(f"Wrote standardized frequency table: {out_path} ({len(freq_df)} rows)")
        result["n_frequency_rows"] = len(freq_df)
    else:
        print(f"NOTE: {FREQUENCY_FILE} not found, skipping frequency standardization.")

    print("\n=== SUMMARY ===")
    for k, v in result.items():
        print(f"{k}: {v}")
    return result


if __name__ == "__main__":
    main()
