# data/raw/ contents

This directory is git-ignored. This README documents what was pulled, when, how, and
what it contains, so the pull is fully reproducible without needing to commit the raw
files themselves.

## Path B (used): FAERS quarterly ASCII files

**Source:** FDA FAERS quarterly ASCII data files, downloaded directly from
`https://fis.fda.gov/content/Exports/faers_ascii_{quarter}.zip`.

**Download date:** 2026-07-08.

**Quarters covered:** 2025Q1, 2025Q2, 2025Q3, 2025Q4, 2026Q1 (all quarters FDA had
published as of the download date; 2026Q2 returned HTTP 404 — not yet released).

**Acquisition script:** `scripts/03_pull_faers.py`.

**Query logic:** for each quarter's `DRUG` table, match `drugname` OR `prod_ai`
(case-insensitive substring) against the synonym list below, then restrict the
analytic cohort to `role_cod == "PS"` (Primary Suspect). Matching `primaryid`s are
then used to subset the `DEMO`, `REAC`, `OUTC`, `THER`, and `INDI` tables for the
same quarters.

### Drug-name synonym verification

Counted against `drugname`/`prod_ai` in the `DRUG` table, any role, across all 5
quarters combined:

| Synonym | Matching drug records (any role) |
|---|---|
| datopotamab deruxtecan | 619 |
| datopotamab deruxtecan-dlnk | 518 |
| datroway | 518 |
| dato-dxd | 0 |
| ds-1062 | 0 |
| ds-1062a | 0 |

Saved to `data/raw/faers/synonym_verification.csv`.

**Note on the zero-match synonyms:** "dato-dxd", "ds-1062", and "ds-1062a" never
appear as a submitted drug name in this data. A manual check of raw hits for the
substring "ds-1062" anywhere in the DRUG table (not just the name fields) turned up
exactly 2 rows, both false leads: free-text dosing annotations referencing
"(DS-1062a Study)" attached to entirely different drugs (epirubicin and
pembrolizumab, in what read as standard neoadjuvant chemo regimens for two cases,
2025Q3 primaryid 255412293 and 2026Q1 primaryid 264139722) — i.e. the case was
associated with the DS-1062a trial's free-text context, but datopotamab deruxtecan
itself was not a listed drug in either case. This confirms the zero counts are
genuine and not a matching bug: reporters use the brand or full generic name, never
the internal development code, when they name the drug.

Role breakdown for the 619 any-role matches: PS 460, SS 154, C 5 (no "I" /
interacting-role matches).

### Analytic cohort (role_cod == "PS")

| Table | Rows | Output file |
|---|---|---|
| DRUG (PS only) | 460 | `data/raw/faers/drug.csv` |
| DEMO | 444 | `data/raw/faers/demo.csv` |
| REAC | 957 | `data/raw/faers/reac.csv` |
| OUTC | 311 | `data/raw/faers/outc.csv` |
| THER | 466 | `data/raw/faers/ther.csv` |
| INDI | 753 | `data/raw/faers/indi.csv` |

444 unique cases (`primaryid`) make up the PS analytic cohort across the 5 quarters.
This is a modest N, consistent with the short post-marketing window — see the
data-window caveat in `CLAUDE.md` and the feasibility note in
`docs/candidate_justification.md`.

### Intermediate artifacts (not analytic outputs)

- `data/raw/faers_quarterly/zips/` — the 5 downloaded quarterly zip files (~333 MB
  total). Cached so re-running the script doesn't re-download.
- `data/raw/faers_quarterly/extracted/` — full extracted ASCII tables for all drugs,
  all quarters (~1.8 GB). Cached so re-running the script doesn't re-extract. Safe to
  delete to reclaim disk space; the script will re-download/re-extract as needed the
  next time it runs (downloads are skipped if the zip is already present; extraction
  is skipped if `.txt` files are already present).

Both are git-ignored (`data/raw/` is excluded wholesale) and were never intended to be
committed.

## Path A (not used this run): OpenVigil FDA / OpenVigil 2.1 export

No OpenVigil export was provided, so Path A was not exercised. It remains available:
`scripts/03_load_openvigil.py` loads and standardizes an OpenVigil export from
`data/raw/openvigil/raw_data.csv` and `data/raw/openvigil/frequency.csv`, and runs
the same synonym-verification logic against whatever drug-name column that export
uses. To use it: query OpenVigil FDA / OpenVigil 2.1 with the synonym list above,
Role = Primary Suspect, window 2025Q1 through the latest available quarter, export
both "Raw data" and "Frequency" as CSV into `data/raw/openvigil/`, then run
`python scripts/03_load_openvigil.py`.

## JADER sensitivity-analysis source

The independent Japanese analysis uses the PMDA Japanese Adverse Drug Event
Report database (JADER) public CSV release `pmdacasereport202606`:

| Source file | Role | Encoding |
|---|---|---|
| `demo202606.csv` | Case/report-version characteristics | CP932 |
| `drug202606.csv` | Drug role, generic name, and reason for use | CP932 |
| `reac202606.csv` | MedDRA/J Preferred Terms, outcome, and onset date | CP932 |
| `hist202606.csv` | Medical history/underlying disease; retained but not required by the current analysis | CP932 |

The release documentation was updated 2026-05-15 and identifies MedDRA/J 29.0.
The analysis treats `identifier + report_count` as one report version; PMDA warns
that the same clinical case may appear more than once when reported by multiple
reporters, so this key must not be described as a unique-patient identifier.

Raw JADER files are not committed. Obtain the release directly from PMDA, accept
the terms supplied with it, and run:

```bash
python scripts/07_jader_comparator.py --jader-dir /path/to/pmdacasereport202606
```

The script writes aggregate comparator tables to `outputs/tables/` and local
UTF-8 cohort extracts to `data/processed/`. The case-level extracts are git-ignored
to avoid secondary redistribution. The PMDA terms also require users publishing
results to notify PMDA in advance and to state explicitly that the results use the
PMDA Japanese Adverse Drug Event Report database (JADER). Repository documentation
is not legal advice; users must review the current terms distributed with their
release.

## Known limitation carried forward

openFDA's API (the fallback path in `scripts/03_pull_faers.py`, not used this run
since the quarterly ASCII files were reachable) cannot distinguish Primary Suspect
from Secondary Suspect — its `drugcharacterization` field only has "Suspect"
(PS+SS combined), "Concomitant", and "Interacting". If a future re-pull has to fall
back to openFDA, any "PS-only" claim in this project would actually be "PS+SS
combined" for that pull, and this must be flagged wherever it's used.
