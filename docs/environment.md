# Environment

Python: 3.14.5 (system `python3`, macOS/arm64)

Virtual environment created with `python3 -m venv .venv`. Activate with:

```
source .venv/bin/activate
```

## Resolved package versions

Installed from `requirements.txt` on 2026-07-08. Exact resolved versions (`pip freeze`):

```
autograd==1.9.1
autograd-gamma==0.5.0
certifi==2026.6.17
charset-normalizer==3.4.9
contourpy==1.3.3
cycler==0.12.1
et_xmlfile==2.0.0
fonttools==4.63.0
formulaic==1.2.2
idna==3.18
interface_meta==2.0.1
kiwisolver==1.5.0
lifelines==0.30.3
matplotlib==3.11.0
narwhals==2.23.0
numpy==2.5.1
openpyxl==3.1.5
packaging==26.2
pandas==2.3.3
patsy==1.0.2
pillow==12.3.0
pyparsing==3.3.2
python-dateutil==2.9.0.post0
pytz==2026.2
PyYAML==6.0.3
requests==2.34.2
scipy==1.18.0
six==1.17.0
statsmodels==0.14.6
tqdm==4.68.4
typing_extensions==4.16.0
tzdata==2026.2
urllib3==2.7.0
wrapt==2.2.2
```

To reproduce this environment exactly, install from this pinned list rather than
`requirements.txt` (which intentionally tracks top-level packages only).

## MedDRA version

FAERS quarterly ASCII extracts do not carry a per-record MedDRA version -- only the
XML extract's `<reactionmeddraversionpt>` tag does (confirmed against the FAERS
Readme.pdf shipped in each quarterly export). Sampled directly from the XML extracts
for the two ends of the analysis window (2026-07-08):

- **2025Q1** (first quarter of the analysis window): MedDRA **27.1**
- **2026Q1** (most recent quarter available at acquisition time): MedDRA **28.1**

MedDRA is upversioned biannually and FAERS re-codes to the current version "per CDER
guidelines" (per the FAERS Readme), so **the MedDRA version is not constant across
this study's analysis window** -- it advanced by two releases (27.1 -> 28.0 -> 28.1)
between 2025Q1 and 2026Q1. The intermediate quarters (2025Q2-Q4) were not individually
checked, since confirming each would require downloading and streaming a ~100+ MB XML
extract per quarter solely for a metadata tag; the two endpoints establish the range.
State this version range (27.1-28.1), not a single fixed version, in the manuscript
methods, and note it as a source of potential PT-coding inconsistency across the study
period if a reaction term was renamed/deprecated/added between MedDRA releases in that
range.
