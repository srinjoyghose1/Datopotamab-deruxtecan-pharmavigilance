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
