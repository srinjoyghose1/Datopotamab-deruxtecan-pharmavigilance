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

(No entries yet — repository just scaffolded.)
