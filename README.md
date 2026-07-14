# Datopotamab Deruxtecan Post-Marketing Pharmacovigilance Signal Study

A reproducible pipeline mining the FDA Adverse Event Reporting System (FAERS) for
disproportionality signals associated with **datopotamab deruxtecan** (Datroway;
Dato-DXd; DS-1062a), a TROP2-directed antibody-drug conjugate, over the first five
quarters of its US marketing history (2025Q1-2026Q1). It analyzes spontaneous
adverse event reports using four disproportionality algorithms, compares recovered
signals against the current FDA label, and characterizes reporting demographics,
serious outcomes, subgroup patterns, and time to onset.

## Background & rationale

Datopotamab deruxtecan couples a humanized IgG1 antibody against TROP2
(trophoblast cell-surface antigen 2) to DXd, a topoisomerase-I inhibitor payload,
via a cleavable linker. It received Japanese approval on 2024-12-27, US FDA
approval for HR-positive/HER2-negative breast cancer on **2025-01-17**,
accelerated approval for EGFR-mutated NSCLC on 2025-06-23, and a first-line
triple-negative breast cancer indication on 2026-05-22 — three US indications
within roughly sixteen months of first approval.

An iterative literature search (documented in full in
[`docs/candidate_justification.md`](docs/candidate_justification.md), conducted
2026-07-08 across PubMed, web search, and preprint servers, using all known name
variants of the drug and adjacent framings) found **no dedicated, full
disproportionality study of this drug** using FAERS, VigiBase, or any other
spontaneous-report database. The nearest overlap — a 2026 class-level FAERS study
of ADC payload components — lists the drug among fourteen approved ADCs but
reports no drug-level result for it, most plausibly because that study's data
window closed at or before this drug's US approval.

**This gap is a direct consequence of the drug's short market history, not
despite it.** The same short window limits this project's own case counts and
statistical power, and this is treated as the central, standing caveat of the
whole study rather than a footnote: every FAERS pull here covers at most
eighteen months since first approval, so disproportionality signals are
necessarily built on far fewer reports than a mature-drug analysis would have,
confidence intervals are correspondingly wide, and this study should be
re-run as additional post-marketing quarters accumulate.

## Methods at a glance

- **Data source & window:** FAERS quarterly ASCII data files (DEMO, DRUG, REAC,
  OUTC, THER, INDI), all quarters published as of 2026-07-08: **2025Q1-2026Q1**
  (2026Q2 not yet released).
- **Drug identification:** six name variants searched case-insensitively against
  `drugname`/`prod_ai`; two abbreviated development-code variants (Dato-DXd,
  DS-1062/DS-1062a) never appeared as a submitted drug name in this data.
- **Case selection:** restricted to `role_cod == "PS"` (Primary Suspect).
- **Deduplication:** FDA-standard rule — within duplicate `caseid`s, keep the
  record with the latest `FDA_DT`, tie-broken by highest `PRIMARYID`; cases on
  any FAERS quarterly deleted-case list are removed.
- **Event coding:** MedDRA Preferred Terms from the REAC table. **MedDRA version
  is not constant across the window** — confirmed 27.1 (2025Q1) advancing to
  28.1 (2026Q1) by direct inspection of the FAERS XML extract (the ASCII files
  used for the main pull carry no version tag at all). PT-to-System-Organ-Class
  mapping was not possible: MedDRA's hierarchy is MSSO-licensed and no crosswalk
  was available.
- **Disproportionality — four algorithms, whole-database background**
  (1,810,135 deduped cases, all drugs, all roles):
  - **ROR** (reporting odds ratio), log-method 95% CI.
  - **PRR** (proportional reporting ratio) with Yates-corrected chi-square.
  - **BCPNN Information Component (IC)** — an approximate Poisson-normal
    delta-method implementation, not the full hierarchical Bayes posterior.
  - **EBGM/EB05** — a single-component Gamma-Poisson shrinker fit by method of
    moments on the drug's own within-drug PT ratio distribution, a simplified
    analogue of the full two-component MGPS mixture.
  - **Consensus signal** = ROR CI-lower >1 (a≥3) **AND** Evans' PRR/chi-square
    (PRR≥2, chi²≥4, a≥3) **AND** IC025 >0. EB05 >2 is a separate, stricter tier,
    not folded into the base consensus rule.
- **Subgroup analysis:** ROR recomputed within sex and age-band strata for the
  consensus signals, against the same background stratified identically, with
  Haldane-Anscombe continuity correction where any cell is zero.
- **Time to onset:** `EVENT_DT` minus the drug's own `START_DT` (joined via
  `drug_seq`, not a concomitant drug's start date); negative or >3,650-day
  (10-year) onsets excluded and logged; Weibull fit (lifelines) for hazard
  classification (β<1 early-failure, β≈1 random, β>1 wear-out).

Full methodological detail, including exact formulas and every documented
approximation, is in [`docs/manuscript.md`](docs/manuscript.md) (Methods) and in
each script's docstring.

## How to reproduce

**1. Environment.**

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Exact resolved package versions used for this analysis are pinned in
[`docs/environment.md`](docs/environment.md).

**2. Raw data.** Raw FAERS pulls are **not** in this repository (see
[Data notes](#data-notes--medlra-version) below) — running `scripts/03_pull_faers.py` downloads
them fresh into `data/raw/faers_quarterly/` and `data/raw/faers/`. No manual
placement is needed for Path B (the path actually used for this analysis). If
you instead have an OpenVigil FDA / OpenVigil 2.1 export, place it at
`data/raw/openvigil/raw_data.csv` and `data/raw/openvigil/frequency.csv` before
running `scripts/03_load_openvigil.py` (Path A).

**3. Run the pipeline, in order.** Candidate selection and the literature
saturation check (conceptually "phases 1-2") were manual research, not scripts —
see `docs/candidate_justification.md`. The numbered pipeline itself is:

| Script | What it does | Key outputs |
|---|---|---|
| `scripts/03_pull_faers.py` | Downloads FAERS quarterly ASCII files, verifies drug-name synonyms, filters to Primary Suspect | `data/raw/faers/*.csv` |
| `scripts/03_load_openvigil.py` | (Alternative Path A) loads/standardizes an OpenVigil export | `data/raw/openvigil/*_standardized.csv` |
| `scripts/04_clean_dedup.py` | FDA-standard dedup, deleted-case removal, builds the analytic dataset | `data/processed/analysis_cases.csv`, `data/processed/case_pt.csv` |
| `scripts/05_disproportionality.py` | ROR/PRR/BCPNN/EBGM against the whole-database background, consensus signal rule | `outputs/tables/signals_all.csv`, `outputs/tables/signals_significant.csv` |
| `scripts/06_stratify_tto.py` | Demographics, serious outcomes, subgroup ROR, SOC rollup, Weibull time-to-onset | `outputs/tables/demo_*.csv`, `subgroup_ror_*.csv`, `tto_summary.csv`, etc. |
| `scripts/07_tables_figures.py` | Renders publication tables (CSV + xlsx) and 300dpi figures | `outputs/tables/T1-T5_*.csv`, `outputs/tables/tables.xlsx`, `outputs/figures/F1-F5_*.png/.svg` |

Each script is runnable independently (`python scripts/0N_*.py`) and prints a
data-provenance summary as it runs. `scripts/05` and `scripts/06` re-derive the
whole-database background from the raw quarterly files, so they take roughly
20-30 seconds each; `scripts/03`'s download step is network-bound (observed
~190 KB/s against FDA's export host).

## Repository map

```
data/
  raw/          # git-ignored FAERS/OpenVigil pulls (README.md tracked; see Data notes)
  processed/    # small, git-tracked analytic datasets (<5MB): analysis_cases.csv, case_pt.csv
scripts/        # numbered pipeline, 03-07 (see table above)
outputs/
  tables/       # all analysis CSVs + T1-T5 publication tables + tables.xlsx
  figures/      # F1-F5, 300dpi PNG + SVG
docs/
  candidate_justification.md   # drug selection + literature saturation check
  environment.md               # pinned package versions + MedDRA version notes
  manuscript.md                # full manuscript draft (Methods/Results/Discussion/Limitations)
refs/           # the two reference papers this project's voice/format follows
CLAUDE.md       # project rules, standing caveats, and a running methodological decisions log
```

## Key results summary

After FDA-standard deduplication, **416 unique cases** with datopotamab
deruxtecan as Primary Suspect were analyzed. Reporters were predominantly female
(76.7%, n=319); 54.1% of cases (n=225) carried a FAERS serious-outcome flag,
including 102 deaths (24.5%). Of **284 distinct PTs** tested against the
whole-database background, **8 met the consensus signal rule** (7 also met the
stricter EB05>2 tier) — split below by whether the PT is a genuine clinical
adverse event or a FAERS administrative/case-context code, not mixed into one
ROR-ranked list (see Discussion for why that split matters):

**Clinical AE signals (4)** — map directly onto the current label:

| PT | a | ROR (95% CI) | PRR | IC (IC025) | EBGM (EB05) | Strict signal |
|---|---|---|---|---|---|---|
| Stomatitis | 64 | 59.12 (45.23-77.26) | 50.18 | 5.17 (3.03) | 49.20 (39.85) | Yes |
| Interstitial lung disease | 16 | 15.20 (9.21-25.07) | 14.65 | 3.37 (1.03) | 14.18 (9.18) | Yes |
| Dry eye | 13 | 7.63 (4.39-13.26) | 7.42 | 2.58 (0.55) | 7.15 (4.40) | Yes |
| Infusion related reaction | 14 | 7.34 (4.31-12.51) | 7.13 | 2.56 (0.61) | 6.88 (4.32) | Yes |

**Administrative/case-context signals (4)** — statistically disproportionate but not adverse events:

| PT | a | ROR (95% CI) | PRR | IC (IC025) | EBGM (EB05) | Strict signal |
|---|---|---|---|---|---|---|
| Disease progression | 86 | 44.12 (34.78-55.98) | 35.21 | 4.87 (3.20) | 34.73 (28.97) | Yes |
| Prescribed underdose | 9 | 19.29 (9.95-37.40) | 18.90 | 3.28 (0.28) | 17.82 (9.86) | Yes |
| No adverse event | 42 | 10.23 (7.43-14.08) | 9.30 | 3.08 (1.75) | 9.18 (7.07) | Yes |
| Off label use | 66 | 2.56 (1.97-3.33) | 2.31 | 1.19 (0.57) | 2.30 (1.87) | No (EB05<2) |

*a = case count; full intervals and all 284 tested PTs, including the
"Signal category" column that drives this split, are in
`outputs/tables/signals_all.csv` and `outputs/tables/T3_top_signals_by_strength.csv`
(T3 is grouped the same way, Clinical AE block then Administrative/case-context
block, each ranked by ROR within).*

Time to onset was computable for 90/416 cases (21.6%); median 22.0 days (IQR
3.0-55.75), Weibull shape β=0.656 (95% CI 0.546-0.766, early-failure/decreasing
hazard). System Organ Class rollup was not possible (no licensed MedDRA
crosswalk); all 284 PTs fall into a single unmapped bucket. Full tables:
`outputs/tables/T1_demographics.csv` through `T5_time_to_onset.csv` (and
`tables.xlsx`). Full figures: `outputs/figures/F1`-`F5` (300dpi PNG + SVG) — F2
and F3 color/group by this same clinical-vs-administrative split, and F5 (the
subgroup heatmap) is restricted to the 4 clinical AE signals only, since a
sex/age breakdown of "Off label use" isn't a clinically meaningful question.

### FAERS breast-cancer comparator analysis

A sensitivity analysis restricted the Dato-DXd target to **196 deduplicated
cases with an explicitly linked breast-cancer indication**. DRUG and INDI were
joined using `primaryid + drug_seq`, preventing an indication belonging to a
different medication from qualifying the case. All explicitly reported breast-
cancer subtypes were included. Three comparator cohorts were evaluated:

| Cohort | N cases |
|---|---:|
| Dato-DXd, breast cancer | 196 |
| Full FAERS excluding all 416 Dato-DXd cases | 1,809,719 |
| Active breast-cancer comparator | 6,092 |
| Active comparator excluding sacituzumab govitecan and trastuzumab deruxtecan | 3,091 |

The active comparator included Primary-Suspect breast-cancer reports for
capecitabine, eribulin, gemcitabine, vinorelbine, paclitaxel/nab-paclitaxel,
carboplatin, sacituzumab govitecan, and trastuzumab deruxtecan. The following
table contains every PT that met the complete ROR + PRR + approximate BCPNN
consensus definition under at least one comparator. Each result is ROR (95% CI),
followed by consensus status.

| PT | Category | Dato cases | Full FAERS excluding Dato | Active breast comparator | Active comparator excluding SG/T-DXd |
|---|---|---:|---:|---:|---:|
| Stomatitis | Clinical AE | 41 | 86.00 (60.90-121.46); **Yes** | 42.14 (26.36-67.38); **Yes** | 38.67 (22.31-67.03); **Yes** |
| Interstitial lung disease | Clinical AE | 8 | 16.17 (7.96-32.83); **Yes** | 0.96 (0.47-1.97); No | 3.42 (1.57-7.43); No |
| Dry eye | Clinical AE | 9 | 11.38 (5.83-22.22); **Yes** | 15.38 (6.87-34.45); **Yes** | 10.58 (4.52-24.76); **Yes** |
| Infusion related reaction | Clinical AE | 8 | 8.97 (4.42-18.20); **Yes** | 4.01 (1.89-8.48); No | 2.44 (1.14-5.20); No |
| Disease progression | Case context | 78 | 111.92 (84.02-149.07); **Yes** | 8.75 (6.46-11.84); **Yes** | 26.58 (18.43-38.35); **Yes** |
| Prescribed underdose | Administrative | 9 | 41.99 (21.48-82.08); **Yes** | 73.25 (22.36-240.01); **Yes** | 148.72 (18.74-1,180.07); **Yes** |
| No adverse event | Administrative | 9 | 4.38 (2.25-8.56); No | 19.50 (8.43-45.13); **Yes** | 148.72 (18.74-1,180.07); **Yes** |
| Off label use | Administrative | 26 | 2.08 (1.37-3.14); No | 2.89 (1.88-4.44); **Yes** | 1.73 (1.12-2.67); No |

Administrative and case-context terms are displayed for transparency but are not
interpreted as adverse drug reactions. Among clinical events, stomatitis and dry
eye were robust across comparator definitions. ILD and infusion-related reaction
were strongly disproportionate against full FAERS but did not meet the complete
consensus rule against the active comparators.

#### ILD 2×2 contingency tables

These tables show the actual case counts underlying the ILD estimates. Each table
accounts for all 196 Dato-DXd breast-cancer cases and every case in its comparator.

**1. Full FAERS excluding all Dato-DXd cases (N=1,809,719)**

| Report group | ILD reported | ILD not reported | Total |
|---|---:|---:|---:|
| Dato-DXd, breast cancer | 8 (a) | 188 (b) | 196 |
| Full FAERS excluding Dato-DXd | 4,750 (c) | 1,804,969 (d) | 1,809,719 |

ROR 16.17 (95% CI 7.96-32.83); PRR 15.55; IC025 0.10; EB05 7.78;
consensus signal: **Yes**.

**2. Active breast-cancer comparator (N=6,092)**

| Report group | ILD reported | ILD not reported | Total |
|---|---:|---:|---:|
| Dato-DXd, breast cancer | 8 (a) | 188 (b) | 196 |
| Active breast-cancer drugs | 259 (c) | 5,833 (d) | 6,092 |

ROR 0.96 (95% CI 0.47-1.97); PRR 0.96; IC025 -1.41; EB05 0.53;
consensus signal: **No**.

**3. Active comparator excluding SG and T-DXd (N=3,091)**

| Report group | ILD reported | ILD not reported | Total |
|---|---:|---:|---:|
| Dato-DXd, breast cancer | 8 (a) | 188 (b) | 196 |
| Active drugs excluding SG/T-DXd | 38 (c) | 3,053 (d) | 3,091 |

ROR 3.42 (95% CI 1.57-7.43); PRR 3.32; IC025 -0.46; EB05 1.66;
consensus signal: **No**.

Complete outputs are in
`outputs/tables/faers_signals_comparator_all.csv` and
`outputs/tables/faers_comparator_cohort_sizes.csv`.

## Discussion & interpretation

Four of the eight consensus signals map directly onto datopotamab deruxtecan's
labeled toxicity profile: **stomatitis**, **interstitial lung disease**, **dry
eye**, and **infusion-related reaction**. Stomatitis was both the strongest
labeled signal (ROR 59.12) and among the most frequently reported PTs (n=64),
consistent with the label's description of stomatitis in 63% of patients in
pooled trial data. Interstitial lung disease reached a comparable ROR (15.20)
on a much smaller case count (n=16), reflecting ILD's far lower background
reporting rate across the whole database relative to stomatitis.

The ILD/pneumonitis signal is directionally consistent with the disproportionality
profile previously described for other deruxtecan-payload ADCs: in a FAERS-based
study of ILD across antibody-drug conjugates (Shi et al., 2024), trastuzumab
deruxtecan showed the strongest ILD correlation among ten ADCs studied
(IC025=5.75). A direct numeric comparison isn't appropriate here — that study used
a different background definition, a different window, and a fully hierarchical
BCPNN rather than this project's approximation — but the shared DXd payload
offers a plausible, unconfirmed mechanistic rationale for the pattern this
analysis also recovered.

The related PT "Pneumonitis" (a=6) did **not** itself reach consensus, failing
only the BCPNN criterion despite clearly passing ROR and Evans' PRR/chi-square —
investigated by hand and confirmed not a computation error, but a known
conservative property of the BCPNN approximation at low counts, compounded here
by the short window. Four of the eight consensus signals (Disease progression,
Off label use, No adverse event, Prescribed underdose) are FAERS administrative
codes rather than adverse events; their statistical disproportionality reflects
reporting conventions in a heavily pretreated oncology population, not a safety
concern, and this is a useful illustration of why automated disproportionality
output requires clinical review before being read as a set of safety findings.
This distinction is now a first-class `signal_category` column in
`outputs/tables/signals_all.csv` (curated against all 284 tested PTs, not just
these four), not just a note in prose — see Key results summary above.

The sex-concentration of several signals in female patients (Stomatitis, Disease
progression, Prescribed underdose, Infusion related reaction) parallels this
cohort's overall sex skew (76.7% female), itself consistent with the drug's
currently breast-cancer-dominated indication mix; whether this reflects a true
sex-differential risk or is fully explained by indication distribution cannot be
determined from these data. The infusion-related-reaction median onset of 0 days
matches the well-established clinical expectation of reactions presenting during
or immediately after infusion, and serves as a useful internal check on the
time-to-onset methodology.

**Clinically**, these findings are directionally consistent with, and do not
suggest revising, the current label's emphasis on early monitoring for
interstitial lung disease/pneumonitis, ocular toxicity, stomatitis, and
infusion-related reactions. This analysis does not identify a candidate
off-label signal warranting a new, specific monitoring recommendation beyond
what the label already specifies.

## Surveillance implications

A natural question for a study like this is whether it actually supports early
safety surveillance for a newly deployed drug, or whether it's only useful in
hindsight. The evidence says yes, with a specific structural caveat.

**How this compares to when dedicated FAERS studies normally appear for a new
ADC:**

| Drug | Approved | First dedicated FAERS disproportionality study | Lag |
|---|---|---|---|
| Trastuzumab deruxtecan (Enhertu) | 2019-12-20 | 2022-10-06, *J Clin Pharm Ther* | ~34 months |
| Sacituzumab govitecan (Trodelvy) | 2020-04-22 | 2023-11-10, *Front Pharmacol* (explicitly billed as "the first real-world safety profile of SG") | ~43 months |
| Datopotamab deruxtecan | 2025-01-17 | This analysis | ~18 months |

Both comparator ADCs are the same drug class and era of FAERS methodology, and
both waited roughly three years for a data window their authors considered
mature before publishing. This project analyzed a new ADC's post-marketing
profile in about half that time by treating the short window as the object of
study rather than a reason to wait — and still recovered four signals
(stomatitis, interstitial lung disease, dry eye, infusion-related reaction)
that map directly onto the eventual label (Discussion, above).

**But this speed exposes a real design problem, not just a benefit.** This
study's own results show it directly: "Pneumonitis" (a=6) passed ROR and Evans'
PRR/chi-square cleanly but failed the approximate BCPNN criterion purely
because its variance is dominated by a tiny expected count at low N — a known
conservative property of the approximation, worse the earlier you look. A
surveillance system that demands full four-algorithm consensus from day one
will systematically under-call real signals during exactly the window when
catching them early matters most.

**The practical implication is a tiered alarm, not a single threshold:**

1. **Watch list** — ROR + Evans' PRR/chi-square alone (cheap, stable even at
   low N) routes a PT to human review, not automatic action.
2. **Priority signal** — full four-algorithm consensus (or the stricter EB05
   tier) escalates faster and with more confidence.
3. **Administrative-artifact triage runs before either tier.** Four of this
   study's eight consensus "signals" were FAERS coding noise (Disease
   progression, Off label use, No adverse event, Prescribed underdose), not
   adverse events — see Discussion. An early-surveillance system built on this
   pattern needs that filter first, or real signals drown in noise exactly
   when reviewer attention is scarcest, right after launch.
4. **The output is a trend, not a verdict.** The value of this approach isn't
   any single snapshot — it's rerunning the same pipeline every quarter and
   watching which watch-list items graduate to priority-signal status as case
   counts grow.

## Limitations

**The short ~18-month observation window is the primary limitation and should
be read first.** Case counts for individual PTs are frequently low (several
strong ROR point estimates rest on fewer than ten reports), confidence intervals
are wide relative to what a mature-drug dataset would produce, and the four
per-AESI time-to-onset subgroups (n=3 to n=8) are too small to support confident
hazard-pattern classification. Every signal here is hypothesis-generating and
should be re-examined as additional post-marketing quarters accumulate.

Beyond the window: FAERS has no denominator of drug exposure, so no incidence or
absolute risk can be estimated, and disproportionality is not causality — a
statistical association does not establish that the drug caused the event.
Reports are voluntary and subject to substantial, non-uniform under-reporting,
and a newly approved, high-profile oncology drug is subject to the Weber effect
and notoriety bias (elevated reporting independent of true event rates).
Indication confounding is directly visible in this study's own results:
"Disease progression" is among the most frequent and most disproportionate
reported terms despite reflecting disease course, not drug toxicity, in a
heavily pretreated oncology population. MedDRA version was not constant across
the window (27.1 to 28.1), which could in principle introduce coding
inconsistency for terms renamed or revised between releases, though this was
not directly assessed; System Organ Class rollup could not be performed at all
for lack of a licensed MedDRA crosswalk. Finally, the BCPNN and EBGM
implementations here are documented approximations of the full published
methods (see Methods above), and time-to-onset results describe only the 21.6%
of the cohort with a computable onset time, not the full 416-case cohort.

## References

Full reference list in [`docs/manuscript.md`](docs/manuscript.md#references).
Key sources: FDA approval notices and DailyMed label for datopotamab deruxtecan;
the 2026 class-level ADC FAERS study that nominally includes but does not report
drug-level results for this drug; Shi et al. (2024), the FAERS-based ADC-class
ILD/pneumonitis study cited above; and the two reference papers in `refs/`
(OpenVigil FDA methods paper; FAERS-based GLP-1RA ADR study) whose presentation
conventions this project's tables, figures, and manuscript follow. The
surveillance-timing comparison above cites two additional sources: the first
dedicated FAERS study of trastuzumab deruxtecan (*J Clin Pharm Ther*,
2022-10-06, data window 2019-12-21 to 2022-06-28, PMC9827941) and of
sacituzumab govitecan (*Front Pharmacol*, 2023-11-10, data window April
2020-March 2023, PMC10667432).

## Data notes & MedDRA version

Raw FAERS pulls are excluded from this repository by size (`data/raw/` is
git-ignored; only its own `README.md`, documenting exact source/date/counts, is
tracked). To re-download: FAERS quarterly ASCII files are published at
`https://fis.fda.gov/content/Exports/faers_ascii_{quarter}.zip` (e.g.
`faers_ascii_2025q1.zip`); `scripts/03_pull_faers.py` automates this download,
extraction, and filtering end to end. See `data/raw/README.md` for the exact
quarters, query strings, and record counts from the pull this analysis used.

MedDRA version is **27.1 to 28.1** across the study window (not a single fixed
version) — see [Methods](#methods-at-a-glance) above and `docs/environment.md`
for how this was confirmed and why the ASCII-only pull can't determine it
per-record.

This analysis is fully reproducible from the scripts in `scripts/` in the order
listed above; every number in this README, the manuscript, and the tables/figures
traces to a script reading a data file (`CLAUDE.md`'s standing rule of full
traceability).
