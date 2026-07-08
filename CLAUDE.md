# Datopotamab Deruxtecan Post-Marketing Pharmacovigilance Signal Study

This project mines the FDA Adverse Event Reporting System (FAERS), accessed via
OpenVigil, for disproportionality signals associated with datopotamab
deruxtecan, an early-stage antibody-drug conjugate. The goal is a reproducible,
fully-traceable pipeline — from raw data pull through cleaned analysis dataset,
statistical signal detection, and figures/tables — suitable for a manuscript
draft written in the register of the two reference papers in `refs/`.

## Drug under study

**Datopotamab deruxtecan** (Datroway; Dato-DXd; DS-1062a) — a TROP2-directed
antibody-drug conjugate (ADC) with a topoisomerase-I inhibitor (DXd) payload.

## Standing rules

1. **Full traceability.** Every number that appears in prose (manuscript,
   protocol, docs) must be traceable to a script reading a data file. No
   numbers typed by hand, no "approximately," no figures pulled from memory
   or from a source other than this repo's own pipeline output.
2. **Academic voice.** Writing matches the register of the two reference
   papers in `refs/` (the OpenVigil FDA methods paper and the FAERS-based
   GLP-1RA ADR study): plain, precise, no AI-isms, no filler, no
   hedging-as-padding. Say what the data shows and stop.
3. **Candor about weak findings.** Report low-N, non-significant, or
   otherwise weak/negative results as such. Do not inflate, spin, or bury
   them. A null or underpowered result is a valid result and must be stated
   plainly.

## Data-window caveat

Datopotamab deruxtecan's first approval was **17 January 2025**. Every FAERS
pull in this project therefore covers a short post-marketing window. This is
an **early post-marketing signal study**:

- All disproportionality signals are hypothesis-generating, not confirmatory.
- State the short observation window as a limitation everywhere it matters:
  protocol, results, discussion, figure/table captions where relevant.
- Do not describe signals using language that implies causal or confirmed
  safety findings.

## Decisions log

Append methodological decisions here as they are made, in the format:

```
### YYYY-MM-DD — short title
Decision and rationale.
```

### 2026-07-08 — Candidate confirmed: datopotamab deruxtecan
Kept datopotamab deruxtecan as the study drug after an iterative saturation check (see
`docs/candidate_justification.md`) found no dedicated FAERS/VigiBase disproportionality
study of the drug; nearest overlap is a 2026 class-level ADC FAERS paper that lists it in
a roster but reports zero drug-level results for it. Re-check for saturation before
manuscript submission given how fast this literature is moving.

### 2026-07-08 — Dedup rule and MedDRA version
FDA-standard dedup applied in `scripts/04_clean_dedup.py`: within duplicate CASEIDs,
keep the record with the latest FDA_DT, tie-break by highest PRIMARYID; cases on any
FDA quarterly deleted-cases list are dropped. 444 raw PS-matched report versions ->
416 unique analytic cases (0 removed as FDA-deleted). Primary Suspect filtering is
applied at acquisition (script 03), before dedup, not after -- see
`data/processed/README.md` note 1 for why this ordering is fine in practice but should
be stated if asked.

MedDRA version is **not constant** across the analysis window: confirmed 27.1 at
2025Q1 up to 28.1 at 2026Q1 by inspecting the FAERS XML extract's
`<reactionmeddraversionpt>` tag directly (the ASCII files used for the actual pull
carry no version tag at all). State this range, not a single version, in the
manuscript methods. MedDRA PT->SOC mapping is not applied (MSSO-licensed, no free
crosswalk available) -- `case_pt.csv`'s `soc` column is null; see
`data/processed/README.md` note 2.

### 2026-07-08 — Disproportionality: consensus rule and method approximations
`scripts/05_disproportionality.py` computes ROR, PRR/chi-square, BCPNN IC, and
MGPS-style EBGM per PT against a whole-database background (all drugs, no role
restriction; N_total = 1,810,135 deduped non-deleted cases, 2025Q1-2026Q1).
**Consensus signal** = ROR CI-lower>1 (a>=3) AND Evans PRR/chi2 (a>=3) AND IC025>0;
EB05>2 is reported as a separate, stricter tier (`strict_signal`), not part of base
consensus. All thresholds are configurable at the top of the script.

Both BCPNN and MGPS are **documented approximations**, not the full published
methods: BCPNN uses a Poisson-normal delta-method approximation to IC/IC025 (not the
full hierarchical Bayes posterior), and MGPS/EBGM uses a single-component
Gamma-Poisson shrinker fit by method-of-moments on this drug's own within-drug
PT ratio distribution (not the full two-component mixture fit by EM across the
whole multi-drug database). State both as approximations in the manuscript methods,
not as full FDA/WHO-grade BCPNN/MGPS output.

Sanity check: of the 4 known label AESIs checked, ILD and stomatitis reached
consensus; pneumonitis and nausea did not. Both were investigated by hand (arithmetic
confirmed correct, not a bug) and are real, expected statistical results: pneumonitis
fails BCPNN only, because the approximate variance is dominated by a tiny expected
count at just 6 cases (a low-count/short-window effect); nausea fails ROR itself,
because it's reported in ~3.8% of ALL FAERS cases regardless of drug and doesn't
separate from that background rate -- a known limitation of whole-database (vs.
custom-comparator) designs for non-specific AEs. Report both candidly, do not adjust
thresholds to force a match.

Also note: 4 of the 8 consensus-signal PTs ("Disease progression", "Off label use",
"No adverse event", "Prescribed underdose") are FAERS administrative/coding
artifacts common in oncology data, not genuine adverse-event signals. They pass the
statistical criteria as computed but should not be presented as safety signals
without saying so -- decide explicitly (and document the decision here) whether to
exclude them via a PT stoplist in a later phase, rather than silently dropping or
silently keeping them.

### 2026-07-08 — Stratification, subgroup ROR, SOC rollup, Weibull TTO
`scripts/06_stratify_tto.py`. Key choices:

- **AESI PT groupings are hand-curated against the label**, not a keyword match --
  e.g. infectious pneumonia variants and lung-cancer/indication terms are excluded
  from ILD/pneumonitis; "Dry mouth" is excluded from Stomatitis (different
  mechanism, xerostomia not mucositis); "Intraocular pressure increased" is excluded
  from Ocular (distinct glaucoma-type finding, not the corneal/surface toxicity the
  label describes). See the AESI_GROUPS comment in the script for the full
  rationale.
- **Age is unknown for 220/416 cases (52.9%)** -- state this plainly whenever an age
  band breakdown is presented; it is not a small-print caveat.
- **Subgroup ROR concentration flags involving "Unknown" sex/age are not
  interpretable as real demographic effects** -- with over half the cohort missing
  age, some PT will concentrate in "Unknown" by base rate alone. The script prints
  this caveat automatically; carry it into the manuscript rather than reporting an
  "Unknown"-concentrated PT as if it were a real subgroup finding.
- **TTO is computed for only 90/416 cases (21.6%)** -- 207 cases have no THER row at
  all for this drug's own drug_seq (as opposed to a concomitant drug's), 36 more have
  only a partial START_DT (year or year-month only, not day-level), and 149 lack
  EVENT_DT. State this completeness rate explicitly everywhere TTO results appear;
  do not present the Weibull fits as if computed on the full cohort.
- **Weibull shape/CI for the 4 AESI subgroups (n=3,4,8,6) are underpowered** --
  flagged via `low_n_caveat` in `tto_summary.csv`. Only the "Overall" TTO estimate
  (n=90) and the Infusion-related-reaction estimate (n=6, but median=0 days matches
  known clinical reality of immediate infusion reactions, a useful sanity check) are
  worth citing with any confidence; the ILD/pneumonitis, Ocular, and Stomatitis
  per-AESI Weibull fits should be reported as exploratory only.
- SOC rollup remains a single "Unmapped" bucket (284/284 PTs) -- same MedDRA
  licensing gap as scripts/04. No change in status.

### 2026-07-08 — Publication tables/figures: presentation conventions
`scripts/07_tables_figures.py`. Table/figure style follows the two reference papers
(plain bordered tables, light-gray header/section bands, N stated in every title;
simple bar/scatter charts with data labels, not dense chart chrome). Colors follow
the project's dataviz-skill palette: one categorical hue per series, status-red
highlight (not a generic series color) for "signal vs not" in the F2 volcano plot,
diverging blue/red centered at ROR=1 for the F5 heatmap.

"Label-listed Y/N" in T2 is the AESI_LABEL_PTS/COMMON_AR_LABEL_PTS sets from the
script -- the label's named AESIs (ILD/pneumonitis, ocular, stomatitis, infusion
reactions) plus its >=20% NSCLC-pool common-AR list, restricted to PTs that actually
string-match our data (verified by hand; several label lab-abnormality terms like
"decreased hemoglobin"/"increased AST" never appear as PTs in this cohort at all --
not forced into the set).

F1 (SOC bar chart) is currently a single bar (Unmapped) -- a literal "one-bar bar
chart" would be a dataviz anti-pattern, so the single bar carries an in-figure
annotation explaining the MedDRA crosswalk gap rather than pretending there's a
real categorical comparison. The code renders a normal multi-bar chart the moment
more than one SOC exists.

Caught and fixed before shipping: F1/F2/F5's x-axis labels initially overlapped the
mandatory N/data-window caption (fixed via a dedicated bottom margin in `savefig`);
F2's "Dry eye"/"Infusion related reaction" labels initially collided (fixed by
alternating label offsets by proximity). Always visually inspect rendered figures,
not just confirm the script exits 0 -- the anti-pattern list in the dataviz skill
does not catch layout collisions, only encoding mistakes.

### 2026-07-08 — Full manuscript draft
`docs/manuscript.md`. Every number in Abstract/Results/Discussion/Limitations was
checked by hand against the T1-T5 tables (or, where a claim needed a number not in
the polished T1-T5 set, against the underlying `signals_all.csv`/`top_pts_*.csv`
output it was cited to explicitly) -- not just written from memory of earlier
turns. This caught two real drafting errors before they shipped: the T2
label-listed count was written backwards (said "seven listed, eight not"; actual
split is 8 Y / 7 N per `T2_top_pts_by_frequency.csv`), and a first-draft Pneumonitis
citation pointed at Table 2 even though Pneumonitis (a=6) isn't in that table's top
15 -- refined to cite `signals_all.csv` directly instead of inventing a Table 2
row that doesn't exist. Lesson: re-verify every number against its cited source at
draft time, even numbers already established earlier in the conversation.

External comparator numbers (Shi et al. 2024 T-DXd ILD/pneumonitis IC025 values, used
qualitatively in Discussion) are cited to that paper, not to our own tables --
confirmed the exact figures via a targeted re-fetch of the source rather than
relying on a paraphrase from earlier phase-05 research notes.

Wrote the Methods section directly from the actual executed pipeline (scripts
03-06 docstrings, prior CLAUDE.md entries, docs/environment.md), not from
`docs/protocol.md` -- that file still does not exist despite being referenced by
name in four separate task requests across this project. It should be written
(retroactively, from this manuscript's Methods section, which is now the closest
thing to a canonical protocol this project has) so future phases have one
authoritative source instead of parameters scattered across chat turns.

American spelling used throughout (characterize, analyze, behavior), per the
project's "consistent, not mixed" spelling rule -- chosen over British because this
is a US-FDA-database study and the OpenVigil reference paper (PLOS ONE) is American,
even though the GLP-1RA reference paper is British-English.
