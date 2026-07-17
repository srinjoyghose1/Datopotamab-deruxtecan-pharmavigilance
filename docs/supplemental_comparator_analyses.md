# Supplemental comparator and algorithm analyses

This supplement records comparator designs and algorithmic sensitivity analyses that are important for auditability but are not treated as co-primary results. All FAERS analyses use the same 182 Dato-DXd Primary-Suspect cases with an explicitly linked breast-cancer indication compatible with HR-positive/HER2-negative disease. Explicit TNBC, HER2-positive, and hormone-receptor-negative indications were excluded; nonspecific, metastatic/recurrent, HER2-low, and incompletely specified breast-cancer terms were retained to avoid discarding otherwise eligible reports.

## Signal rule and interpretation

The primary rule is deliberately limited to two complementary quantities:

- ROR lower 95% confidence limit greater than 1, with at least three Dato-DXd cases; and
- approximate BCPNN lower credibility bound IC025 greater than 0.

PRR with Yates-corrected chi-square and the approximate Gamma-Poisson EBGM/EB05 remain computed for sensitivity analysis. They are not additional votes in the primary rule. ROR and PRR are highly correlated transformations of the same 2×2 table, so requiring both creates the appearance of independent corroboration without providing it. The local EBGM implementation is also the least fully validated component: it is a single-component, within-drug approximation rather than the production two-component MGPS model. The four-algorithm column below therefore answers a sensitivity question only: would a primary ROR+IC signal also clear Evans' PRR/chi-square criteria and EB05 >2?

Comparator restriction can reduce confounding by indication, but it also changes the clinical question and can introduce masking when comparator drugs share the event of interest. These trade-offs are described by Alkabbani and Gamble (2023), Gravel and Douros (2023), Gravel, Bai, and Douros (2024), and Bai, Douros, and Gravel (2025). Stratified or subgroup analyses can improve relevance but lose power rapidly when cell counts are small, as shown by Seabroke et al. (2016). Accordingly, no single comparator is presented as universally correct.

## Prespecified primary comparator hierarchy

| Tier | Comparator | N reports | Role |
|---|---|---:|---|
| 1 | Full FAERS excluding all 416 Dato-DXd reports | 1,809,719 | Broad screening background and continuity with the primary analysis |
| 2 | Trial-aligned chemotherapy: capecitabine, eribulin, vinorelbine, gemcitabine | 940 | Primary active comparator; closest observable analogue to investigator's-choice chemotherapy in the pivotal trial |
| 3a | Trastuzumab deruxtecan (T-DXd) alone | 1,630 | Drug-specific ADC contextual analysis |
| 3b | Sacituzumab govitecan (SG) alone | 348 | Drug-specific TROP2-ADC contextual analysis |

The active-drug cohorts require Primary-Suspect status and a breast-cancer indication linked to that drug through `primaryid + drug_seq`. The hierarchy is interpretive rather than a sequence of hypothesis tests. Tier 2 is the principal clinical benchmark; Tier 1 shows broad reporting disproportionality; Tiers 3a and 3b show whether findings change against individual ADC backgrounds.

### Why no joint or pooled ADC comparator

T-DXd (N=1,630) and SG (N=348) are not pooled as a primary comparator. The solo arms isolate two competing mechanistic hypotheses: T-DXd tests context shared through the DXd payload, whereas SG tests context shared through the TROP2 target. A pooled ADC arm would conflate those hypotheses. At most it could bound the magnitude of masking without attributing it to payload or target, which is strictly less informative than retaining the two solo estimates. The drugs also differ in treatment setting, market tenure, utilization, and adverse-event profile. The older nine-drug pool is retained below only as a sensitivity analysis.

## Four-algorithm sensitivity across primary tiers

The tables include every PT meeting the primary ROR+IC rule in each analysis. `Four algorithms` means that the PT also met Evans' PRR/chi-square criteria and EB05 >2. Administrative and disease-course terms are retained for transparency and are not interpreted as adverse drug reactions. Full all-PT results, including every non-signal row, are supplied in the cited machine-readable CSVs.

### Primary 416-case whole-database analysis

| PT | a | c | ROR (95% CI) | PRR; chi-square | IC (IC025) | EBGM (EB05) | Four algorithms |
|---|---:|---:|---|---|---|---|---|
| Stomatitis | 64 | 5,549 | 59.12 (45.23-77.26) | 50.17; 3,010.18 | 5.17 (3.03) | 49.19 (39.85) | Yes |
| Disease progression | 86 | 10,626 | 44.12 (34.78-55.98) | 35.21; 2,818.26 | 4.87 (3.20) | 34.73 (28.97) | Yes |
| Prescribed underdose | 9 | 2,072 | 19.29 (9.95-37.40) | 18.90; 134.74 | 3.28 (0.28) | 17.82 (9.86) | Yes |
| Interstitial lung disease | 16 | 4,750 | 15.20 (9.21-25.07) | 14.65; 189.98 | 3.37 (1.03) | 14.18 (9.18) | Yes |
| No adverse event | 42 | 19,649 | 10.23 (7.43-14.08) | 9.30; 305.50 | 3.08 (1.75) | 9.18 (7.06) | Yes |
| Dry eye | 13 | 7,622 | 7.63 (4.39-13.26) | 7.42; 66.10 | 2.58 (0.55) | 7.15 (4.40) | Yes |
| Infusion related reaction | 14 | 8,548 | 7.34 (4.31-12.51) | 7.12; 67.93 | 2.55 (0.61) | 6.88 (4.32) | Yes |
| Off label use | 66 | 124,227 | 2.56 (1.97-3.33) | 2.31; 51.29 | 1.19 (0.57) | 2.29 (1.86) | No |

The complete 284-PT table is `outputs/tables/signals_all.csv`.

### Indication-compatible comparator tiers

| Comparator | PT | a | c | ROR (95% CI) | PRR; chi-square | IC (IC025) | EBGM (EB05) | Four algorithms |
|---|---|---:|---:|---|---|---|---|---|
| Full FAERS | Disease progression | 78 | 10,626 | 126.98 (94.62-170.42) | 72.99; 5,458.99 | 5.64 (3.36) | 71.99 (59.51) | Yes |
| Full FAERS | Stomatitis | 40 | 5,549 | 91.59 (64.42-130.21) | 71.68; 2,706.34 | 5.25 (2.47) | 70.26 (53.72) | Yes |
| Full FAERS | Prescribed underdose | 9 | 2,072 | 45.39 (23.19-88.83) | 43.19; 328.88 | 3.74 (0.26) | 40.72 (22.53) | Yes |
| Full FAERS | Dry eye | 8 | 7,622 | 10.87 (5.35-22.09) | 10.44; 59.34 | 2.75 (0.05) | 9.85 (5.24) | Yes |
| Full FAERS | Infusion related reaction | 8 | 8,548 | 9.69 (4.77-19.69) | 9.31; 51.49 | 2.64 (0.03) | 8.78 (4.67) | Yes |
| Trial chemotherapy | Prescribed underdose | 9 | 0 | 102.99 (5.97-1,777.83) | 97.70; 40.23 | 2.28 (0.06) | 5.15 (3.09) | Yes |
| Trial chemotherapy | Stomatitis | 40 | 10 | 26.20 (12.81-53.56) | 20.66; 151.77 | 2.23 (1.17) | 4.78 (3.69) | Yes |
| Trial chemotherapy | Disease progression | 78 | 51 | 13.07 (8.70-19.64) | 7.90; 206.29 | 1.87 (1.18) | 3.71 (3.08) | Yes |
| Trial chemotherapy | Off label use | 19 | 11 | 9.84 (4.60-21.07) | 8.92; 46.84 | 1.86 (0.48) | 3.81 (2.63) | Yes |
| T-DXd alone | Stomatitis | 40 | 4 | 114.51 (40.39-324.63) | 89.56; 317.26 | 3.04 (1.69) | 8.62 (6.63) | Yes |
| T-DXd alone | Prescribed underdose | 9 | 1 | 84.75 (10.67-672.91) | 80.60; 62.53 | 2.66 (0.18) | 7.48 (4.35) | Yes |
| T-DXd alone | Disease progression | 78 | 149 | 7.45 (5.31-10.46) | 4.69; 166.78 | 1.75 (1.09) | 3.42 (2.83) | Yes |
| T-DXd alone | Off label use | 19 | 32 | 5.82 (3.23-10.50) | 5.32; 39.96 | 1.79 (0.44) | 3.68 (2.52) | Yes |
| SG alone | Stomatitis | 40 | 1 | 97.75 (13.31-717.88) | 76.48; 75.76 | 1.47 (0.61) | 2.71 (2.11) | Yes |
| SG alone | Off label use | 19 | 2 | 20.17 (4.64-87.61) | 18.16; 28.02 | 1.34 (0.14) | 2.47 (1.74) | No |
| SG alone | Disease progression | 78 | 51 | 4.37 (2.88-6.63) | 2.92; 50.09 | 0.81 (0.28) | 1.77 (1.47) | No |

The SG comparison demonstrates why four-algorithm unanimity is not the primary definition: Off label use and Disease progression meet the ROR+IC rule but fail the approximate EB05 threshold. This is a sensitivity discrepancy, not evidence that either case-context term is a clinical toxicity.

Complete per-PT results, including negative PTs and every algorithm, are in `outputs/tables/faers_signals_full_faers_excluding_all_dato.csv`, `faers_signals_trial_aligned_chemo.csv`, `faers_signals_trastuzumab_deruxtecan_alone.csv`, and `faers_signals_sacituzumab_govitecan_alone.csv`.

## Leave-one-drug-out analysis of the trial-aligned pool

Leave-one-drug-out analyses were limited to the four PTs positive under the base trial-aligned comparator. A result is stable when the ROR lower limit remains above 1 and IC025 remains above 0 after each component drug is removed.

| PT | Drug removed | Comparator N | a | ROR (95% CI) | IC025 | Primary signal |
|---|---|---:|---:|---|---:|---|
| Disease progression | Capecitabine | 457 | 78 | 33.52 (16.78-66.98) | 0.97 | Yes |
| Stomatitis | Capecitabine | 457 | 40 | 31.90 (11.22-90.71) | 0.74 | Yes |
| Off label use | Capecitabine | 457 | 19 | 13.20 (4.43-39.38) | 0.23 | Yes |
| Prescribed underdose | Capecitabine | 457 | 9 | 50.10 (2.90-865.47) | -0.23 | No |
| Stomatitis | Eribulin | 558 | 40 | 25.92 (10.77-62.33) | 0.84 | Yes |
| Disease progression | Eribulin | 558 | 78 | 7.46 (4.94-11.25) | 0.69 | Yes |
| Off label use | Eribulin | 558 | 19 | 5.80 (2.70-12.43) | 0.11 | Yes |
| Prescribed underdose | Eribulin | 558 | 9 | 61.16 (3.54-1,056.25) | -0.14 | No |
| Prescribed underdose | Gemcitabine | 872 | 9 | 95.55 (5.54-1,649.38) | 0.03 | Yes |
| Stomatitis | Gemcitabine | 872 | 40 | 24.28 (11.87-49.65) | 1.11 | Yes |
| Disease progression | Gemcitabine | 872 | 78 | 15.20 (9.89-23.36) | 1.20 | Yes |
| Off label use | Gemcitabine | 872 | 19 | 14.40 (5.96-34.82) | 0.55 | Yes |
| Prescribed underdose | Vinorelbine | 933 | 9 | 102.23 (5.92-1,764.60) | 0.06 | Yes |
| Stomatitis | Vinorelbine | 933 | 40 | 26.00 (12.72-53.16) | 1.17 | Yes |
| Disease progression | Vinorelbine | 933 | 78 | 12.97 (8.63-19.49) | 1.18 | Yes |
| Off label use | Vinorelbine | 933 | 19 | 9.77 (4.56-20.91) | 0.48 | Yes |

Stomatitis, Disease progression, and Off label use remain positive in every deletion. Prescribed underdose is unstable because its Dato-DXd count is only nine and the trial comparator contains no such reports; its IC025 becomes negative when capecitabine or eribulin is removed. The machine-readable table is `outputs/tables/faers_signals_leave_one_out.csv`, and Figure 8 shows the same estimates.

## Demoted comparator analyses

### Expanded chemotherapy tier

The expanded chemotherapy comparator adds paclitaxel, nab-paclitaxel, and carboplatin to the four trial drugs (N=2,129). It is useful for checking sensitivity to a wider chemotherapy background, but is less directly aligned with the pivotal trial and mixes treatment settings.

| PT meeting primary rule | a | c | ROR (95% CI) | IC025 | Four algorithms |
|---|---:|---:|---|---:|---|
| Prescribed underdose | 9 | 0 | 233.20 (13.52-4,023.78) | 0.24 | Yes |
| Stomatitis | 40 | 12 | 49.69 (25.50-96.84) | 1.75 | Yes |
| Disease progression | 78 | 65 | 23.82 (16.23-34.95) | 1.85 | Yes |

### Pooled nine-drug tier

The former active comparator pooled the seven chemotherapy drugs with T-DXd and SG (N=4,107). It is demoted because the mixture obscures which drug or class drives the denominator and creates substantial opportunity for event-specific masking.

| PT meeting primary rule | a | c | ROR (95% CI) | IC025 | Four algorithms |
|---|---:|---:|---|---:|---|
| Prescribed underdose | 9 | 1 | 213.61 (26.91-1,695.53) | 0.28 | Yes |
| Stomatitis | 40 | 17 | 67.77 (37.51-122.46) | 2.08 | Yes |
| Dry eye | 8 | 15 | 12.54 (5.25-29.98) | 0.00 | Yes |
| Disease progression | 78 | 265 | 10.87 (7.90-14.96) | 1.59 | Yes |
| Off label use | 19 | 161 | 2.86 (1.73-4.71) | 0.08 | No |

After excluding the two ADCs, this pool is mathematically identical to the expanded chemotherapy tier (N=2,129); it is therefore not a separate analysis layer.

The complete all-PT demoted tables are `outputs/tables/faers_signals_expanded_chemo.csv`, `outputs/tables/faers_signals_active_breast_comparator.csv`, and `outputs/tables/faers_signals_active_breast_comparator_class_exclusion.csv`.

## ILD contingency tables for demoted comparators

These are the former README tables, retained here so the exact 2×2 evidence remains auditable.

### Full FAERS excluding all Dato-DXd reports

| Report group | ILD | No ILD | Total |
|---|---:|---:|---:|
| Dato-DXd, compatible breast cancer | 7 | 175 | 182 |
| Full FAERS excluding Dato-DXd | 4,750 | 1,804,969 | 1,809,719 |

ROR 15.20 (95% CI 7.14-32.37); PRR 14.65; IC025 -0.10; EB05 6.94; primary signal: **No**.

### Pooled nine-drug comparator

| Report group | ILD | No ILD | Total |
|---|---:|---:|---:|
| Dato-DXd, compatible breast cancer | 7 | 175 | 182 |
| Pooled active drugs | 178 | 3,929 | 4,107 |

ROR 0.88 (95% CI 0.41-1.91); PRR 0.89; IC025 -1.58; EB05 0.49; primary signal: **No**.

### Expanded chemotherapy / pooled ADC-excluded comparator

| Report group | ILD | No ILD | Total |
|---|---:|---:|---:|
| Dato-DXd, compatible breast cancer | 7 | 175 | 182 |
| Expanded chemotherapy | 28 | 2,101 | 2,129 |

ROR 3.00 (95% CI 1.29-6.97); PRR 2.92; IC025 -0.67; EB05 1.42; primary signal: **No**.

## References

1. Alkabbani W, Gamble JM. Active-comparator restricted disproportionality analysis for pharmacovigilance signal detection studies of chronic disease medications: an example using sodium/glucose cotransporter 2 inhibitors. *Br J Clin Pharmacol.* 2023;89:431-439. doi:10.1111/bcp.15178.
2. Gravel CA, Douros A. Considerations on the use of different comparators in pharmacovigilance: a methodological review. *Br J Clin Pharmacol.* 2023;89:2671-2676. doi:10.1111/bcp.15802.
3. Gravel CA, Bai W, Douros A. Comparators in pharmacovigilance: a quasi-quantification bias analysis. *Drug Saf.* 2024;47:809-819. doi:10.1007/s40264-024-01433-5.
4. Bai W, Douros A, Gravel CA. Masking in active comparator designs in pharmacovigilance: a retrospective bias analysis on the spontaneous reporting of thiazolidinediones and cardiovascular events. *Pharmacoepidemiol Drug Saf.* 2025;34:e70102. doi:10.1002/pds.70102.
5. Seabroke S, Candore G, Juhlin K, et al. Performance of stratified and subgrouped disproportionality analyses in spontaneous databases. *Drug Saf.* 2016;39:355-364. doi:10.1007/s40264-015-0388-3.
