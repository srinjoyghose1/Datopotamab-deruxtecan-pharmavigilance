"""
Path B (programmatic): acquire FAERS case-level data for datopotamab deruxtecan.

Primary source: FDA FAERS quarterly ASCII files (DEMO, DRUG, REAC, OUTC, THER, INDI),
2025Q1 through the latest quarter FDA has published. Records are matched on drug name
synonyms in the DRUG file and restricted to Role = Primary Suspect (role_cod == "PS")
for the analytic cohort; unrestricted (any role) counts are reported alongside for the
synonym-verification step.

Falls back to the openFDA API (https://api.fda.gov/drug/event.json) if the quarterly
ASCII files cannot be downloaded. The openFDA fallback cannot distinguish Primary
Suspect (PS) from Secondary Suspect (SS) -- openFDA's drugcharacterization field only
distinguishes Suspect (1, PS+SS combined) from Concomitant (2) and Interacting (3) --
so results pulled via the fallback are documented as "Suspect (PS+SS)", not "PS only".

Outputs land in data/raw/faers/ (git-ignored). Nothing here is committed.
"""

import subprocess
import zipfile
from pathlib import Path

import pandas as pd
import requests

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = REPO_ROOT / "data" / "raw" / "faers_quarterly"
ZIP_DIR = RAW_DIR / "zips"
EXTRACT_DIR = RAW_DIR / "extracted"
OUT_DIR = REPO_ROOT / "data" / "raw" / "faers"

SYNONYMS = [
    "datopotamab deruxtecan",
    "datopotamab deruxtecan-dlnk",
    "datroway",
    "dato-dxd",
    "ds-1062",
    "ds-1062a",
]

QUARTER_URL_TMPL = "https://fis.fda.gov/content/Exports/faers_ascii_{quarter}.zip"
START_YEAR, START_Q = 2025, 1
TABLES = ["DEMO", "DRUG", "REAC", "OUTC", "THER", "INDI"]


def quarter_sequence(start_year=START_YEAR, start_q=START_Q, max_quarters=12):
    year, q = start_year, start_q
    for _ in range(max_quarters):
        yield f"{year}q{q}"
        q += 1
        if q > 4:
            q = 1
            year += 1


def quarter_available(quarter: str) -> tuple[bool, int]:
    # fis.fda.gov times out under python-requests (TLS/handshake quirk observed
    # 2026-07-08) but is reliably reachable via curl, so shell out for all HTTP
    # traffic against this host.
    url = QUARTER_URL_TMPL.format(quarter=quarter)
    result = subprocess.run(
        ["curl", "-sL", "-D", "-", "-o", "/dev/null", url, "--max-time", "30", "-r", "0-0"],
        capture_output=True, text=True,
    )
    headers = result.stdout
    if "206" in headers.splitlines()[0] if headers else False:
        for line in headers.splitlines():
            if line.lower().startswith("content-range"):
                total = int(line.split("/")[-1].strip())
                return True, total
        return True, 0
    return False, 0


def detect_available_quarters() -> list[str]:
    available = []
    for quarter in quarter_sequence():
        ok, size = quarter_available(quarter)
        if ok:
            print(f"  {quarter}: available ({size / 1e6:.1f} MB)")
            available.append(quarter)
        else:
            print(f"  {quarter}: not available, stopping detection")
            break
    return available


def download_quarter(quarter: str) -> Path:
    ZIP_DIR.mkdir(parents=True, exist_ok=True)
    dest = ZIP_DIR / f"faers_ascii_{quarter}.zip"
    if dest.exists():
        print(f"  {quarter}: zip already present, skipping download")
        return dest
    url = QUARTER_URL_TMPL.format(quarter=quarter)
    tmp = dest.with_suffix(".zip.part")
    print(f"  {quarter}: downloading via curl (this host is slow, ~5-10 min/quarter)...")
    subprocess.run(
        ["curl", "-sL", "-o", str(tmp), url, "--max-time", "1800"],
        check=True,
    )
    tmp.rename(dest)
    return dest


def extract_quarter(quarter: str, zip_path: Path) -> Path:
    dest = EXTRACT_DIR / quarter
    if dest.exists() and any(dest.rglob("*.txt")):
        print(f"  {quarter}: already extracted, skipping")
        return dest
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(dest)
    return dest


def find_table_file(extract_root: Path, table: str, quarter: str) -> Path | None:
    candidates = list(extract_root.rglob(f"{table}*.txt")) + list(
        extract_root.rglob(f"{table}*.TXT")
    )
    if not candidates:
        return None
    return candidates[0]


def load_table(extract_root: Path, table: str, quarter: str) -> pd.DataFrame:
    path = find_table_file(extract_root, table, quarter)
    if path is None:
        print(f"    WARNING: {table} file not found for {quarter}")
        return pd.DataFrame()
    df = pd.read_csv(
        path, delimiter="$", low_memory=False, encoding="latin1", dtype=str
    )
    df.columns = [c.strip().lower() for c in df.columns]
    df["quarter"] = quarter
    return df


def synonym_match_mask(drug_df: pd.DataFrame) -> pd.DataFrame:
    """Return a DataFrame of per-synonym boolean match columns against drugname/prod_ai."""
    name_col = drug_df["drugname"].fillna("").str.lower() if "drugname" in drug_df else pd.Series([""] * len(drug_df))
    ai_col = drug_df["prod_ai"].fillna("").str.lower() if "prod_ai" in drug_df else pd.Series([""] * len(drug_df))
    masks = {}
    for syn in SYNONYMS:
        s = syn.lower()
        masks[syn] = name_col.str.contains(s, regex=False) | ai_col.str.contains(s, regex=False)
    return pd.DataFrame(masks, index=drug_df.index)


def run_quarterly_ascii_path() -> dict | None:
    print("Detecting available FAERS quarterly ASCII files...")
    quarters = detect_available_quarters()
    if not quarters:
        print("No quarterly ASCII files reachable.")
        return None

    print(f"\nAcquiring {len(quarters)} quarter(s): {quarters}")
    drug_frames, other_frames = [], {t: [] for t in TABLES if t != "DRUG"}

    for quarter in quarters:
        print(f"\n=== {quarter} ===")
        zip_path = download_quarter(quarter)
        extract_root = extract_quarter(quarter, zip_path)
        drug_df = load_table(extract_root, "DRUG", quarter)
        if drug_df.empty:
            continue
        drug_frames.append(drug_df)
        for table in other_frames:
            other_frames[table].append(load_table(extract_root, table, quarter))

    if not drug_frames:
        print("No DRUG records loaded from any quarter.")
        return None

    all_drug = pd.concat(drug_frames, ignore_index=True)

    # --- synonym verification (unrestricted by role) ---
    match_df = synonym_match_mask(all_drug)
    print("\nSynonym match counts (drug records, any role, 2025Q1-latest):")
    synonym_counts = {}
    for syn in SYNONYMS:
        n = int(match_df[syn].sum())
        synonym_counts[syn] = n
        print(f"  {syn!r}: {n} drug records")

    any_match = match_df.any(axis=1)
    matched_any_role = all_drug[any_match]
    print(f"\nTotal drug records matching ANY synonym (any role): {len(matched_any_role)}")
    print("Breakdown by role_cod:")
    if "role_cod" in matched_any_role:
        print(matched_any_role["role_cod"].value_counts().to_string())

    # --- analytic cohort: Primary Suspect only ---
    ps_mask = any_match & (all_drug.get("role_cod", "") == "PS")
    matched_ps = all_drug[ps_mask]
    print(f"\nTotal drug records matching ANY synonym AND role_cod == 'PS': {len(matched_ps)}")
    primaryids = set(matched_ps["primaryid"].unique())
    print(f"Unique primaryids (cases) in PS analytic cohort: {len(primaryids)}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    matched_ps.to_csv(OUT_DIR / "drug.csv", index=False)

    for table, frames in other_frames.items():
        if not frames:
            continue
        combined = pd.concat(frames, ignore_index=True)
        subset = combined[combined["primaryid"].isin(primaryids)]
        subset.to_csv(OUT_DIR / f"{table.lower()}.csv", index=False)
        print(f"  saved {table.lower()}.csv: {len(subset)} rows")

    verification_df = pd.DataFrame(
        [{"synonym": s, "matches_any_role": c} for s, c in synonym_counts.items()]
    )
    verification_df.to_csv(OUT_DIR / "synonym_verification.csv", index=False)

    return {
        "source": "faers_quarterly_ascii",
        "quarters": quarters,
        "n_ps_records": len(matched_ps),
        "n_cases": len(primaryids),
        "synonym_counts": synonym_counts,
    }


def run_openfda_fallback() -> dict:
    print("\nFalling back to openFDA API (patient-level, Suspect = PS+SS combined)...")
    base = "https://api.fda.gov/drug/event.json"
    synonym_counts = {}
    for syn in SYNONYMS:
        q = f'patient.drug.medicinalproduct:"{syn}"'
        r = requests.get(base, params={"search": q, "limit": 1}, timeout=30)
        if r.status_code == 200:
            n = r.json()["meta"]["results"]["total"]
        else:
            n = 0
        synonym_counts[syn] = n
        print(f"  {syn!r}: {n} records (patient.drug.medicinalproduct)")

    or_query = " OR ".join(f'patient.drug.medicinalproduct:"{s}"' for s in SYNONYMS)
    full_query = f'({or_query}) AND receivedate:[20250101 TO 20261231]'

    all_records = []
    skip = 0
    limit = 100
    while True:
        r = requests.get(
            base, params={"search": full_query, "limit": limit, "skip": skip}, timeout=60
        )
        if r.status_code == 404:
            break
        r.raise_for_status()
        payload = r.json()
        results = payload.get("results", [])
        if not results:
            break
        all_records.extend(results)
        skip += limit
        total = payload["meta"]["results"]["total"]
        print(f"  fetched {len(all_records)}/{total}")
        if skip >= total or skip >= 25000:  # openFDA hard skip limit is 25000
            break

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.json_normalize(all_records).to_csv(OUT_DIR / "openfda_events.csv", index=False)

    verification_df = pd.DataFrame(
        [{"synonym": s, "matches_any_role": c} for s, c in synonym_counts.items()]
    )
    verification_df.to_csv(OUT_DIR / "synonym_verification.csv", index=False)

    return {
        "source": "openfda_api",
        "query": full_query,
        "n_records": len(all_records),
        "synonym_counts": synonym_counts,
    }


def main():
    result = run_quarterly_ascii_path()
    if result is None:
        result = run_openfda_fallback()
    print("\n=== SUMMARY ===")
    for k, v in result.items():
        print(f"{k}: {v}")
    return result


if __name__ == "__main__":
    main()
