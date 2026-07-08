# Contributing / Reproduction Notes

This is a research pipeline, not a general-purpose library. This file covers how
to work in it; see [`README.md`](README.md) for what it does and its results.

## Quick reproduction

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/03_pull_faers.py        # downloads FAERS data fresh (not in this repo)
python scripts/04_clean_dedup.py
python scripts/05_disproportionality.py
python scripts/06_stratify_tto.py
python scripts/07_tables_figures.py
```

Each script is idempotent and prints a data-provenance summary as it runs. See
the table in `README.md` ("How to reproduce") for what each one consumes and
produces.

## Standing rules for any new work in this repo

These are enforced by convention, not tooling, so apply them by hand (see
`CLAUDE.md` for the full rationale behind each):

1. **Full traceability.** Every number in prose (manuscript, README, docs) must
   trace to a script reading a data file — no hand-typed figures, no numbers
   pulled from memory of an earlier finding. Before citing a number, re-derive
   or re-read it from the actual output CSV, even if you're confident you
   remember it correctly.
2. **Academic voice.** Match the register of the two papers in `refs/` — plain,
   precise, no filler, no hedging-as-padding.
3. **Candor about weak findings.** Report low-N, non-significant, or
   conflicting results as such. Do not adjust thresholds or reframe language to
   make a weak result look stronger.
4. **State the short observation window wherever it matters.** This is an early
   post-marketing signal study; every disproportionality result is
   hypothesis-generating, not causal or confirmatory.

## Adding a new pipeline phase

Follow the existing numbering convention in `scripts/` (`0N_description.py`).
New scripts should:
- Print a data-provenance summary (inputs consumed, rows in/out, outputs written).
- Write outputs to `outputs/tables/` or `outputs/figures/`, not back into `data/raw/`.
- Log any new methodological choice (threshold, approximation, exclusion rule) in
  `CLAUDE.md`'s decisions log, dated, with the *why* — not just the *what*.

## Data

Raw FAERS pulls are not committed (`data/raw/` is git-ignored — see
`data/raw/README.md` for exact source/date/counts and `README.md`'s "Data notes"
section for the re-download URL). `data/processed/` holds small, git-tracked
derived datasets only (<5MB each); verify size by hand before committing anything
new there, since `.gitignore` can't filter by file size.

## License / reuse

No LICENSE file is included yet — this is a private, in-progress research
repository, not a released package. Do not assume permission to reuse, copy, or
redistribute any part of this repo (code, data, or manuscript draft) until a
LICENSE is added.
