# Bias inventory

This inventory records material threats to interpretation and the analysis layer used to address them. “Addressed” means a planned diagnostic was implemented; it does not mean the underlying bias was eliminated.

| Bias or design risk | Status | Analysis response | Residual limitation |
|---|---|---|---|
| Confounding by indication and oncology disease severity | Addressed | Restricted target and active comparators to drug-linked, HR-positive/HER2-negative-compatible breast-cancer indications | Nonspecific breast-cancer terms remain, and FAERS lacks staging, line, and reliable disease-severity variables |
| Explicitly discordant TNBC, HER2-positive, or HR-negative target reports | Addressed | Excluded at case level in FAERS and report-version level in JADER | Unreported receptor status cannot be inferred |
| Choice of active comparator | Addressed | Prespecified full-FAERS, trial-aligned chemotherapy, and individual T-DXd/SG tiers | No spontaneous-report comparator can balance prescribing channel, market tenure, or exposure |
| One trial-aligned drug dominating the pooled result | Addressed | Implemented leave-one-drug-out analysis for every PT positive in the base trial pool; see `docs/supplemental_comparator_analyses.md` and Gravel, Bai, and Douros (2024), doi:10.1007/s40264-024-01433-5 | Deletions reduce comparator N and can destabilize Bayesian lower bounds |
| Masking from pooled ADCs | Addressed | Analyzed T-DXd and SG separately; demoted the joint nine-drug pool | Individual ADC cohorts, especially SG, remain small |
| Correlated ROR and PRR counted as independent confirmation | Addressed | Primary rule uses ROR lower CI and IC025 only; PRR/chi-square is sensitivity analysis | ROR and IC still arise from the same underlying report table |
| Approximate EBGM treated as production MGPS | Addressed | EBGM/EB05 demoted to sensitivity status and explicitly labeled approximate | A validated two-component, whole-database MGPS model was not fitted |
| Sparse-event instability | Partly addressed | Minimum a≥3, confidence/credibility lower bounds, shrinkage outputs, and leave-one-out analysis | Wide intervals and conservative IC025 values persist in low-count strata |
| Duplicate reports within FAERS | Addressed | FDA-standard CASEID/FDA_DT/PRIMARYID deduplication and deleted-case removal | Undetectable clinical duplicates submitted under different CASEIDs may remain |
| Possible FAERS-JADER overlap | Open | Databases analyzed separately | No shared identifier exists to verify cross-database duplication |
| Under-reporting, stimulated reporting, and Weber/notoriety effects | Open | Explicitly treated results as reporting associations and displayed quarterly reporting trend (Figure 6) | Exposure denominators and the counterfactual number of unreported events are unavailable |
| Missing or incorrectly linked indication | Partly addressed | Required INDI linkage through `primaryid + drug_seq` | Reports without a usable indication are excluded; submitted indication can be incomplete or miscoded |
| Concomitant therapy and comorbidity confounding | Open | No defensible adjustment possible with current fields | FAERS does not provide complete longitudinal treatment or clinical histories |
| MedDRA version changes across quarters | Open | Version range documented; PT-level source terms preserved | No licensed hierarchy/cross-version normalization was available |
| Missing time-to-onset data | Open | Reported completeness and analyzed only validated nonnegative intervals ≤10 years | Onset subset may not represent all Dato-DXd reports |
| Cross-database coding and reporting differences | Partly addressed | Kept JADER and FAERS separate and retained Japanese JADER PTs as authoritative | ROR differences cannot be interpreted as geographic risk differences |

## Leave-one-out decision

The leave-one-drug-out item moved from **Open** to **Addressed** on 2026-07-17. Stomatitis, Disease progression, and Off label use retained the primary ROR+IC signal under all four deletions. Prescribed underdose did not: IC025 fell below zero when capecitabine or eribulin was removed. Because Prescribed underdose is an administrative term based on nine target cases, this instability does not alter the clinical safety interpretation.
