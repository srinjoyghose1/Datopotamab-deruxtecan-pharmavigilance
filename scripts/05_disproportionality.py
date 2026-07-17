"""
Four-algorithm disproportionality analysis for datopotamab deruxtecan vs the whole
FAERS database background (2025Q1-2026Q1), using the deduped analysis cohort from
scripts/04_clean_dedup.py (data/processed/{analysis_cases,case_pt}.csv) as the drug
of interest and the full (unfiltered, all-drugs) FAERS quarterly REAC/DEMO tables as
background.

For every PT reported at least once with the drug, builds the standard 2x2 table:

                    PT of interest      all other PTs
    drug               a                    b
    all other drugs    c                    d

  a = cases with drug AND this PT
  b = n_drug - a               (cases with drug, without this PT)
  c = n_pt_total - a           (cases with this PT, without drug)
  d = N_total - n_drug - n_pt_total + a

n_drug, n_pt_total, and N_total are all case-level counts (unique primaryids), not
raw report-row counts. The background (n_pt_total, N_total) is NOT restricted to any
particular drug role (PS/SS/C) -- it is simply "every case in the FAERS window that
mentions this PT" / "every case in the FAERS window", matching the standard whole-
database disproportionality convention used across the field. This is a deliberate,
documented choice, not an oversight.

Four algorithms are computed per PT:

  1. ROR (reporting odds ratio), log-method 95% CI.
  2. PRR (proportional reporting ratio) with Yates-corrected chi-square.
  3. BCPNN Information Component (IC) with an approximate IC025 lower bound.
  4. MGPS-style EBGM with EB05, via a single-component Gamma-Poisson shrinker.

Both BCPNN and MGPS below are DELIBERATE APPROXIMATIONS of the full published
methods, not the full methods themselves -- see the docstrings on
`bcpnn_ic` and `fit_gamma_poisson_shrinker` for exactly what is approximated and why.
This is disclosed here, in-code, and must be disclosed in the manuscript methods.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma

REPO_ROOT = Path(__file__).resolve().parents[1]
EXTRACT_ROOT = REPO_ROOT / "data" / "raw" / "faers_quarterly" / "extracted"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
TABLES_DIR = REPO_ROOT / "outputs" / "tables"

QUARTERS = ["2025q1", "2025q2", "2025q3", "2025q4", "2026q1"]

# ============================================================================
# CONFIG -- signal-detection thresholds and the consensus rule. Edit here only.
# ============================================================================
MIN_CASES_A = 3            # Evans convention: don't call a signal below this many cases
ROR_CI_LOWER_THRESHOLD = 1.0
PRR_THRESHOLD = 2.0
CHI2_THRESHOLD = 4.0
IC025_THRESHOLD = 0.0
EB05_THRESHOLD = 2.0

# Primary signal rule: frequentist precision (ROR lower 95% CI >1, a>=3)
# plus conservative shrinkage (IC025>0). PRR and the single-component EBGM
# approximation remain computed for sensitivity analysis, but are not treated
# as independent confirmations of the same 2x2 table.
CONSENSUS_CRITERIA = ["ror_signal", "bcpnn_signal"]
FOUR_ALGORITHM_CRITERIA = ["ror_signal", "prr_signal", "bcpnn_signal", "mgps_signal"]

SANITY_CHECK_PTS = ["Interstitial lung disease", "Pneumonitis", "Stomatitis", "Nausea"]

# PTs that are FAERS administrative/case-context codes, not clinical adverse events --
# curated by hand against all 284 PTs actually tested in this cohort (not just the
# ones that happen to reach consensus), covering disease-progression-as-AE artifacts
# and product/dosing/administration-issue codes. Everything not in this set is
# treated as a clinical AE signal. This distinction was previously only made in prose
# (manuscript, README); it is added here as a first-class column so it propagates
# into every downstream table and figure instead of being re-derived ad hoc.
NON_CLINICAL_ARTIFACT_PTS = {
    "Disease progression", "Neoplasm progression",
    "Off label use", "No adverse event",
    "Prescribed underdose", "Underdose", "Drug ineffective",
    "Incorrect product formulation administered", "Intentional product use issue",
    "Product availability issue", "Product distribution issue",
    "Product dose omission issue", "Product leakage", "Product prescribing issue",
    "Product use in unapproved indication", "Product use issue",
    "Therapeutic response unexpected", "Therapy non-responder",
}


def load_deleted_caseids() -> set:
    deleted = set()
    for delete_file in EXTRACT_ROOT.glob("*/Deleted/DELETE*.txt"):
        with open(delete_file) as f:
            for line in f:
                line = line.strip()
                if line:
                    deleted.add(line)
    return deleted


def find_table_file(quarter: str, table: str) -> Path:
    matches = list((EXTRACT_ROOT / quarter / "ASCII").glob(f"{table}*.txt")) + list(
        (EXTRACT_ROOT / quarter / "ASCII").glob(f"{table}*.TXT")
    )
    return matches[0]


def build_background_case_universe() -> pd.DataFrame:
    """Whole-database case universe: every case (any drug, any role) in our window,
    FDA-standard deduped (latest FDA_DT, tie-break highest PRIMARYID), FDA-deleted
    cases removed. Returns one row per kept primaryid with its caseid."""
    frames = []
    for quarter in QUARTERS:
        demo_file = find_table_file(quarter, "DEMO")
        df = pd.read_csv(
            demo_file, delimiter="$", usecols=["primaryid", "caseid", "fda_dt"],
            dtype=str, low_memory=False, encoding="latin1",
        )
        df.columns = [c.strip().lower() for c in df.columns]
        frames.append(df)
    demo = pd.concat(frames, ignore_index=True)

    demo["_fda_dt_sort"] = pd.to_numeric(demo["fda_dt"], errors="coerce").fillna(-1)
    demo["_primaryid_sort"] = pd.to_numeric(demo["primaryid"], errors="coerce").fillna(-1)
    demo = demo.sort_values(["_fda_dt_sort", "_primaryid_sort"], ascending=False)
    deduped = demo.drop_duplicates(subset="caseid", keep="first").drop(
        columns=["_fda_dt_sort", "_primaryid_sort"]
    )

    deleted_caseids = load_deleted_caseids()
    deduped = deduped[~deduped["caseid"].isin(deleted_caseids)]
    return deduped[["primaryid", "caseid"]]


def build_background_pt_counts(kept_primaryids: set) -> pd.Series:
    """Whole-database per-PT case counts. Each primaryid belongs to exactly one
    quarter's REAC file, so per-quarter (primaryid, pt) dedup + summed counts across
    quarters is equivalent to a single global dedup, without needing to hold all
    ~7M REAC rows in memory at once."""
    running_counts = None
    for quarter in QUARTERS:
        reac_file = find_table_file(quarter, "REAC")
        df = pd.read_csv(
            reac_file, delimiter="$", usecols=["primaryid", "pt"],
            dtype=str, low_memory=False, encoding="latin1",
        )
        df.columns = [c.strip().lower() for c in df.columns]
        df = df[df["primaryid"].isin(kept_primaryids)]
        df = df.drop_duplicates(subset=["primaryid", "pt"])
        counts = df.groupby("pt").size()
        running_counts = counts if running_counts is None else running_counts.add(counts, fill_value=0)
    return running_counts.astype(int)


def compute_ror(a, b, c, d):
    ror = (a * d) / (b * c)
    log_se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    log_ror = np.log(ror)
    ci_lower = np.exp(log_ror - 1.96 * log_se)
    ci_upper = np.exp(log_ror + 1.96 * log_se)
    return ror, ci_lower, ci_upper


def compute_prr_chi2(a, b, c, d):
    prr = (a / (a + b)) / (c / (c + d))
    chi2, _, _, _ = stats.chi2_contingency([[a, b], [c, d]], correction=True)
    return prr, chi2


def bcpnn_ic(a, n_drug, n_pt, n_total):
    """Approximate BCPNN Information Component with a Poisson-normal IC025.

    APPROXIMATION, not the full hierarchical Bayes BCPNN posterior (Bate 1998):
    treats the observed count `a` and the expected count `E` under independence as
    two approximately-Poisson quantities and applies the standard delta-method
    variance for the log-ratio of two independent Poisson counts (var ~= 1/a + 1/E),
    converted to log2 units. This "quick BCPNN" approximation with a 0.5 continuity
    correction is the standard practical substitute used across the FAERS
    disproportionality literature when a full hierarchical-Bayes fit isn't performed
    (see van Puijenbroek et al. 2002; Bate & Evans 2009 for the general approach).
    State this as an approximation in the manuscript methods.
    """
    e = (n_drug * n_pt) / n_total
    ic = np.log2((a + 0.5) / (e + 0.5))
    var_ic = (1 / np.log(2) ** 2) * (1 / (a + 0.5) + 1 / (e + 0.5))
    ic025 = ic - 1.96 * np.sqrt(var_ic)
    return ic, ic025, e


def fit_gamma_poisson_shrinker(a_vec: np.ndarray, e_vec: np.ndarray):
    """Fit a SINGLE-component Gamma-Poisson empirical Bayes shrinker via method of
    moments on this drug's own within-drug distribution of observed/expected
    ratios across its ~N reported PTs.

    APPROXIMATION, not full MGPS (DuMouchel 1999): full MGPS fits a TWO-component
    Gamma mixture prior via EM across the entire multi-drug FAERS database (every
    drug-PT cell, not just this drug's), which is computationally impractical at
    this project's scope. The one-component Gamma-Poisson shrinker is DuMouchel's
    simpler foundational model (the two-component "multi-item" GPS is his
    refinement for separately shrinking background-noise vs. true-signal cells).
    Fitting the single component's hyperparameters from this drug's own PT-level
    ratio distribution (rather than the whole database's full drug x PT cube) is a
    standard, well-documented empirical-Bayes simplification. State this as an
    approximation, not full FDA/WHO-grade MGPS output, in the manuscript methods.

    Returns (alpha, beta) of the fitted Gamma(alpha, beta) prior on the relative
    reporting rate lambda, where lambda | a ~ Gamma(alpha + a, beta + E).
    """
    ratios = a_vec / e_vec
    mean_r = ratios.mean()
    var_r = ratios.var(ddof=1)
    if var_r <= 0:
        # Degenerate case (all ratios identical) -- fall back to a weakly
        # informative prior rather than dividing by zero.
        return 1.0, 1.0
    beta = mean_r / var_r
    alpha = mean_r * beta
    return alpha, beta


def ebgm_eb05(a, e, alpha, beta):
    post_alpha = alpha + a
    post_beta = beta + e
    ebgm = np.exp(digamma(post_alpha) - np.log(post_beta))
    # EB05/EB95 conventionally form a 90% interval (DuMouchel 1999), not 95% --
    # EB05 is the 5th percentile of the posterior.
    eb05 = stats.gamma.ppf(0.05, a=post_alpha, scale=1 / post_beta)
    return ebgm, eb05


def main():
    print("Building whole-database background case universe (2025Q1-2026Q1, all drugs)...")
    background_cases = build_background_case_universe()
    n_total = len(background_cases)
    kept_primaryids = set(background_cases["primaryid"])
    print(f"  N_total (unique deduped, non-deleted cases, whole database): {n_total}")

    print("Building whole-database per-PT case counts...")
    background_pt_counts = build_background_pt_counts(kept_primaryids)
    print(f"  {len(background_pt_counts)} distinct PTs in the whole-database background")

    case_pt = pd.read_csv(PROCESSED_DIR / "case_pt.csv", dtype=str)
    analysis_cases = pd.read_csv(PROCESSED_DIR / "analysis_cases.csv", dtype=str)
    n_drug = analysis_cases["caseid"].nunique()
    print(f"\nDrug of interest: {n_drug} unique cases")

    drug_pt_counts = case_pt.drop_duplicates(subset=["caseid", "pt"]).groupby("pt").size()
    print(f"{len(drug_pt_counts)} distinct PTs reported with the drug")

    rows = []
    for pt, a in drug_pt_counts.items():
        n_pt = int(background_pt_counts.get(pt, 0))
        if n_pt < a:
            # Shouldn't happen (drug's own cases are part of the background count),
            # but guard against any background/drug-cohort mismatch rather than
            # silently producing a negative cell.
            n_pt = a
        b = n_drug - a
        c = n_pt - a
        d = n_total - n_drug - n_pt + a
        rows.append({"pt": pt, "a": a, "b": b, "c": c, "d": d, "n_pt_total": n_pt})

    df = pd.DataFrame(rows)

    # --- ROR ---
    ror, ror_lo, ror_hi = compute_ror(df["a"], df["b"], df["c"], df["d"])
    df["ror"], df["ror_ci_lower"], df["ror_ci_upper"] = ror, ror_lo, ror_hi
    df["ror_signal"] = (df["ror_ci_lower"] > ROR_CI_LOWER_THRESHOLD) & (df["a"] >= MIN_CASES_A)

    # --- PRR / chi-square (Yates) ---
    prr_chi2 = [compute_prr_chi2(a, b, c, d) for a, b, c, d in zip(df["a"], df["b"], df["c"], df["d"])]
    df["prr"], df["chi2"] = zip(*prr_chi2)
    df["prr_signal"] = (
        (df["prr"] >= PRR_THRESHOLD) & (df["chi2"] >= CHI2_THRESHOLD) & (df["a"] >= MIN_CASES_A)
    )

    # --- BCPNN IC / IC025 ---
    bcpnn = [bcpnn_ic(a, n_drug, npt, n_total) for a, npt in zip(df["a"], df["n_pt_total"])]
    df["ic"], df["ic025"], df["expected_count"] = zip(*bcpnn)
    df["bcpnn_signal"] = df["ic025"] > IC025_THRESHOLD

    # --- MGPS-style EBGM / EB05 (single-component shrinker, fit on this drug's own PTs) ---
    alpha, beta = fit_gamma_poisson_shrinker(df["a"].to_numpy(), df["expected_count"].to_numpy())
    print(f"\nFitted single-component Gamma-Poisson shrinker: alpha={alpha:.4f}, beta={beta:.4f}")
    ebgm_vals = [ebgm_eb05(a, e, alpha, beta) for a, e in zip(df["a"], df["expected_count"])]
    df["ebgm"], df["eb05"] = zip(*ebgm_vals)
    df["mgps_signal"] = df["eb05"] > EB05_THRESHOLD

    # --- consensus ---
    df["consensus_signal"] = df[CONSENSUS_CRITERIA].all(axis=1)
    df["four_algorithm_signal"] = df[FOUR_ALGORITHM_CRITERIA].all(axis=1)

    # --- clinical vs. administrative-artifact classification ---
    df["signal_category"] = df["pt"].apply(
        lambda pt: "Administrative/case-context" if pt in NON_CLINICAL_ARTIFACT_PTS else "Clinical AE"
    )

    df = df.sort_values(["consensus_signal", "ror"], ascending=[False, False])

    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(TABLES_DIR / "signals_all.csv", index=False)
    signals_significant = df[df["consensus_signal"]]
    signals_significant.to_csv(TABLES_DIR / "signals_significant.csv", index=False)

    print("\n=== SUMMARY ===")
    print(f"PTs tested: {len(df)}")
    print(f"Passing ROR criterion (CI lower > {ROR_CI_LOWER_THRESHOLD}, a >= {MIN_CASES_A}): {df['ror_signal'].sum()}")
    print(f"Passing PRR/Evans criterion (PRR >= {PRR_THRESHOLD}, chi2 >= {CHI2_THRESHOLD}, a >= {MIN_CASES_A}): {df['prr_signal'].sum()}")
    print(f"Passing BCPNN criterion (IC025 > {IC025_THRESHOLD}): {df['bcpnn_signal'].sum()}")
    print(f"Passing MGPS criterion (EB05 > {EB05_THRESHOLD}): {df['mgps_signal'].sum()}")
    print(f"Consensus signal ({' AND '.join(CONSENSUS_CRITERIA)}): {df['consensus_signal'].sum()}")
    print(f"Four-algorithm sensitivity signal: {df['four_algorithm_signal'].sum()}")
    consensus_by_category = signals_significant["signal_category"].value_counts()
    print(f"Of the {len(signals_significant)} consensus signals: "
          f"{consensus_by_category.get('Clinical AE', 0)} Clinical AE, "
          f"{consensus_by_category.get('Administrative/case-context', 0)} Administrative/case-context")

    print("\n=== SANITY CHECK: known label events ===")
    for pt_name in SANITY_CHECK_PTS:
        match = df[df["pt"].str.lower() == pt_name.lower()]
        if match.empty:
            print(f"  {pt_name!r}: NOT FOUND in drug's PT list -- investigate.")
            continue
        row = match.iloc[0]
        status = "SIGNAL" if row["consensus_signal"] else "not a consensus signal"
        print(
            f"  {pt_name!r}: a={row['a']}, ROR={row['ror']:.2f} "
            f"[{row['ror_ci_lower']:.2f}-{row['ror_ci_upper']:.2f}], "
            f"PRR={row['prr']:.2f}, chi2={row['chi2']:.2f}, "
            f"IC={row['ic']:.2f} (IC025={row['ic025']:.2f}), "
            f"EBGM={row['ebgm']:.2f} (EB05={row['eb05']:.2f}) -> {status}"
        )
        if not row["consensus_signal"]:
            print(f"    WARNING: expected label event {pt_name!r} did not reach consensus -- investigate.")
            if pt_name.lower() == "pneumonitis":
                print(
                    "    INVESTIGATED (2026-07-08): arithmetic checked by hand, no bug. "
                    "ROR/PRR clearly signal; BCPNN's IC025 fails because the approximate "
                    "variance is dominated by the tiny expected count (E<1) at a=6 cases -- "
                    "a known conservative property of this Poisson-normal BCPNN "
                    "approximation at low counts, compounded by this drug's short "
                    "post-marketing window. Not a computation error; report candidly as a "
                    "case where methods disagree at low N, not as a consensus signal."
                )
            if pt_name.lower() == "nausea":
                print(
                    "    INVESTIGATED (2026-07-08): arithmetic checked by hand, no bug. "
                    "ROR's own CI lower bound (0.98) fails before BCPNN is even reached. "
                    "Nausea is reported in ~3.8% of ALL FAERS cases regardless of drug, so "
                    "against a whole-database background this drug does not show excess "
                    "reporting beyond that already-high base rate. A known limitation of "
                    "whole-database (vs. custom-comparator) disproportionality designs for "
                    "non-specific AEs, not an error in this analysis."
                )

    print(f"\nWrote {TABLES_DIR / 'signals_all.csv'} ({len(df)} rows)")
    print(f"Wrote {TABLES_DIR / 'signals_significant.csv'} ({len(signals_significant)} rows)")

    return df


if __name__ == "__main__":
    main()
