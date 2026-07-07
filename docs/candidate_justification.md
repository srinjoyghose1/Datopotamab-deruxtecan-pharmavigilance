# Candidate Justification: Datopotamab Deruxtecan

Search conducted 2026-07-08 against PubMed, Google Scholar (via web search), Frontiers,
ScienceDirect, medRxiv, and ResearchGate. Searches were run iteratively across all known
name variants of the drug and several adjacent framings (see Section 2). All web-derived
facts below are cited to a primary or near-primary source; none are asserted from memory.

## 1. Drug profile

**Datopotamab deruxtecan** (Datroway; Dato-DXd; DS-1062a; datopotamab deruxtecan-dlnk) is
a TROP2 (trophoblast cell-surface antigen 2)-directed antibody-drug conjugate. A humanized
IgG1 monoclonal antibody against TROP2 is joined by a cleavable linker to DXd, an
exatecan-derivative topoisomerase I inhibitor payload. After the antibody binds
TROP2-expressing tumor cells and the conjugate is internalized, DXd is released
intracellularly, causing DNA damage and cell death; bystander killing of adjacent
TROP2-low/negative cells is part of the proposed mechanism ([FDA, NSCLC accelerated
approval notice](https://www.fda.gov/drugs/resources-information-approved-drugs/fda-grants-accelerated-approval-datopotamab-deruxtecan-dlnk-egfr-mutated-non-small-cell-lung-cancer)).
Recommended dose is 6 mg/kg IV every 3 weeks (capped at 540 mg for patients ≥90 kg).

### Approval timeline

| Date | Action | Indication | Basis |
|---|---|---|---|
| 2024-12-27 | Japan approval | HR+/HER2− unresectable/recurrent breast cancer | [Datopotamab Deruxtecan: First Approval, PubMed 40323341](https://pubmed.ncbi.nlm.nih.gov/40323341/) |
| 2025-01-17 | US FDA approval (regular) | Unresectable/metastatic HR+, HER2− breast cancer | [FDA approval notice](https://www.fda.gov/drugs/resources-information-approved-drugs/fda-approves-datopotamab-deruxtecan-dlnk-unresectable-or-metastatic-hr-positive-her2-negative-breast) |
| 2025-06-23 | US FDA accelerated approval | Locally advanced/metastatic EGFR-mutated NSCLC after prior EGFR-directed therapy and platinum chemo | [FDA accelerated approval notice](https://www.fda.gov/drugs/resources-information-approved-drugs/fda-grants-accelerated-approval-datopotamab-deruxtecan-dlnk-egfr-mutated-non-small-cell-lung-cancer) |
| 2026-05-22 | US FDA approval | First-line, unresectable/metastatic TNBC, PD-1/PD-L1-ineligible | [FDA approval notice](https://www.fda.gov/drugs/resources-information-approved-drugs/fda-approves-datopotamab-deruxtecan-dlnk-unresectable-or-metastatic-triple-negative-breast-cancer) |

Two of three US indications (NSCLC, TNBC) were added within the ~16 months following
first approval — the drug's utilization base is still expanding as of this writing
(2026-07-08).

### Label warnings and precautions

Per the current DATROWAY US label (DailyMed, [setid
2950227c-6230-4ca4-a135-46e44d9424a0](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=2950227c-6230-4ca4-a135-46e44d9424a0)),
there is **no boxed warning**. Section 5 (Warnings and Precautions) contains four
numbered subsections:

- **5.1 Interstitial lung disease (ILD) / pneumonitis** — severe, life-threatening, or
  fatal cases can occur; reported in 7% of NSCLC trial patients and 3% of breast cancer
  trial patients.
- **5.2 Ocular adverse reactions** — dry eye, keratitis, conjunctivitis and related
  findings, affecting 38% of patients across pooled trials.
- **5.3 Stomatitis** — including mouth ulcers and oral mucositis, in 63% of patients in
  pooled safety data.
- **5.4 Embryo-fetal toxicity** — no human pregnancy data; risk inferred from mechanism
  (genotoxic topoisomerase inhibitor payload).

Infusion-related reactions are addressed in the label's dosing/administration guidance
(extended monitoring for the first two cycles) rather than as a fifth numbered
Warnings and Precautions subsection — worth noting precisely, since this project's
protocol should not overstate infusion reactions as an FDA-designated "warning" if we
cite the label structure directly. It is nonetheless a recognized adverse event of
special interest (AESI) in the clinical literature and will be tracked as such.

Most common adverse reactions (≥20%, NSCLC pool): stomatitis, nausea, alopecia, fatigue,
decreased hemoglobin, decreased lymphocytes, constipation, increased calcium, increased
AST, decreased WBC, increased LDH, musculoskeletal pain, decreased appetite, increased
ALT, rash.

## 2. Novelty / saturation check

This check was run iteratively, not as a single query, specifically to avoid asserting a
novelty gap that a broader search would falsify. Rounds:

1. Direct: "datopotamab deruxtecan FAERS disproportionality"
2. Alternate names: "Dato-DXd FAERS", "DS-1062" / "DS-1062a" FAERS/VigiBase, "Datroway
   FAERS VigiBase"
3. Adjacent framings: TROP2 ADC safety disproportionality, ADC-class pneumonitis/ILD
   VigiBase/FAERS, ADC ocular toxicity FAERS, ADC cardiac AE FAERS, ADC stomatitis
   disproportionality
4. Time-boxed re-checks: explicit 2026 date qualifiers, medRxiv/preprint search,
   VigiBase/Uppsala Monitoring Centre-specific search, real-world claims-database search

**Finding: no dedicated, full disproportionality study of datopotamab deruxtecan using
FAERS, VigiBase, JADER, or any other spontaneous-report pharmacovigilance database
currently exists in the indexed literature.** Specifically:

- **PubMed** search for datopotamab deruxtecan returns pharmacokinetics, clinical-trial
  safety summaries, approval summaries, and adverse-event *management* papers (e.g.
  [Lisberg et al., "Datopotamab deruxtecan-associated select adverse events," Oncologist
  2025](https://pubmed.ncbi.nlm.nih.gov/40700616/); [clinical management/prophylaxis of
  AESIs, Cancer Treat Rev
  2024](https://pubmed.ncbi.nlm.nih.gov/38502995/)) — these are expert-opinion and
  institutional-protocol pieces synthesizing clinician experience and trial data, **not**
  disproportionality analyses of spontaneous-report data. No FAERS/VigiBase methodology
  (ROR, PRR, IC, EBGM) is used in any of them.
- The stomatitis-specific literature on Dato-DXd (Delphi consensus, prophylaxis
  guidelines — e.g. [US expert Delphi consensus, PMC12325563](https://pmc.ncbi.nlm.nih.gov/articles/PMC12325563/))
  is similarly consensus/guideline-based, not database-derived.
- **Near-miss, reported candidly:** one 2026 Frontiers paper, ["A real-world study of
  adverse event profiles associated with the four key components of antibody-drug
  conjugates based on the FAERS
  database"](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2026.1702195/full),
  compiles a roster of the 14 FDA-approved ADCs as of 2024-12-31 and explicitly lists
  datopotamab deruxtecan in that roster (FAERS window: 2004Q1–2024Q4). However, on direct
  inspection **no drug-level case counts, ROR, PRR, IC, or EBGM values are reported for
  datopotamab deruxtecan anywhere in the paper** — the paper's positive-signal results
  table covers only 13 of the 14 ADCs, and datopotamab deruxtecan is absent from it. The
  paper does not state why. The most likely explanation is that this study's FAERS
  cutoff (2024Q4) falls at or just before the drug's 2025-01-17 US approval date, leaving
  too few US reports to clear a signal-detection case-count threshold. Regardless of the
  reason, **no drug-specific result for Dato-DXd exists even in the one paper that
  nominally included it.**
- **Class-level ADC-ILD studies** using FAERS/VigiBase/JADER/CVAR were checked directly
  for inclusion of Dato-DXd and confirmed **not to include it**, in every case because
  their data windows end at or before Dato-DXd's approval:
  - [Shi et al., "Interstitial lung disease with antibody-drug conjugates... FAERS
    2014–2023," Ther Adv Med Oncol
    2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11635890/) — 10 ADCs, window ends
    2023-03-31, confirmed absent.
  - Cardiac-AE ADC FAERS study (PMC12825936) — 9 ADCs (Enhertu, Polivy, Padcev, Trodelvy,
    Tivdak, Akalux, Zynlonta, RC48, Elahere), window 2019-01 to 2023-09, confirmed
    absent.
  - Breast-cancer ILD disproportionality study (PMC11957809) — window 2004–2023, three
    ADCs (T-DM1, T-DXd, sacituzumab govitecan), confirmed absent.
  - Three-database (FAERS/CVAR/JADER) ovarian-cancer ILD study (PMC12823856) — windows
    through 2024–2025, only mirvetuximab soravtansine among ADCs, confirmed absent.
- **VigiBase-specific and claims-database searches** for datopotamab, Dato-DXd, DS-1062,
  DS-1062a, and Datroway returned no matching studies as of 2026-07-08.

**Conclusion: the novelty gap is real, current, and narrow.** No one has published a
dedicated multi-algorithm FAERS or VigiBase disproportionality analysis of datopotamab
deruxtecan. The only paper that even lists the drug in a database-derived roster reports
no results for it. This is consistent with — not despite — the fact that the drug has
only been on the US market since January 2025: existing FAERS-based ADC studies are
mechanically excluded from covering it by their own data windows. That same fact is this
project's central limitation (Section 3) and must be treated as such, not talked around.

This gap should be re-checked before manuscript submission, since the literature moves
quickly in this space (three papers found here were published in 2026 alone).

## 3. Justification

**Rationale.** Datopotamab deruxtecan is a newly approved ADC with a rapidly expanding
label (three indications added within ~16 months) and a labeled toxicity profile (ILD/
pneumonitis, ocular toxicity, stomatitis) that is qualitatively well characterized from
trials but not yet characterized in real-world spontaneous-report data at the
population level. This is precisely the situation in which early spontaneous-report
signal detection has the most marginal value: trial populations are selected and
monitored under protocol; post-marketing populations are not. No dedicated multi-
algorithm FAERS study of this drug exists (Section 2), so this project fills an
identifiable, current gap rather than replicating existing work.

**Feasibility note.** Because the observation window is short (US approval
2025-01-17 to the present), raw FAERS report counts for datopotamab deruxtecan will be
modest relative to mature drugs. This has two direct methodological consequences for the
protocol:
- Bayesian shrinkage methods (BCPNN, MGPS/EBGM) should be prioritized over, or at minimum
  reported alongside, frequentist disproportionality measures (ROR, PRR), since shrinkage
  estimators are specifically designed to stabilize signal estimates under low report
  counts and reduce false-positive signals driven by small cells.
- Multi-algorithm consensus (requiring agreement across ≥2 methods before calling a
  signal, following the convention used in the class-level ADC literature reviewed
  above) should be the default signal-calling rule, not a single-method threshold.
- All signals produced under these conditions are hypothesis-generating. None should be
  described as confirmatory, and the manuscript should not use causal language for any
  individual disproportionality result.

**Limitation, stated plainly.** The short post-marketing window is below what would be
considered a mature-data threshold for spontaneous-report pharmacovigilance. Small
report counts widen confidence intervals substantially and increase the risk that any
signal reflects the Weber effect / notoriety bias — heightened reporting attention paid
to a new, high-profile oncology drug — rather than a true elevated adverse event rate.
This is not a caveat to mention once in a limitations paragraph; it constrains how every
signal in this study should be interpreted and must be restated wherever a specific
signal is discussed (results, discussion, abstract, and any table/figure caption
presenting a signal estimate).

## References

- FDA. [FDA approves datopotamab deruxtecan-dlnk for unresectable or metastatic, HR-positive, HER2-negative breast cancer](https://www.fda.gov/drugs/resources-information-approved-drugs/fda-approves-datopotamab-deruxtecan-dlnk-unresectable-or-metastatic-hr-positive-her2-negative-breast).
- FDA. [FDA grants accelerated approval to datopotamab deruxtecan-dlnk for EGFR-mutated non-small cell lung cancer](https://www.fda.gov/drugs/resources-information-approved-drugs/fda-grants-accelerated-approval-datopotamab-deruxtecan-dlnk-egfr-mutated-non-small-cell-lung-cancer).
- FDA. [FDA approves datopotamab deruxtecan-dlnk for unresectable or metastatic triple-negative breast cancer](https://www.fda.gov/drugs/resources-information-approved-drugs/fda-approves-datopotamab-deruxtecan-dlnk-unresectable-or-metastatic-triple-negative-breast-cancer).
- DailyMed. [DATROWAY - datopotamab deruxtecan injection, powder, lyophilized, for solution](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=2950227c-6230-4ca4-a135-46e44d9424a0).
- [Datopotamab Deruxtecan: First Approval. PubMed 40323341](https://pubmed.ncbi.nlm.nih.gov/40323341/).
- Lisberg A, et al. [Datopotamab deruxtecan-associated select adverse events. PubMed 40700616](https://pubmed.ncbi.nlm.nih.gov/40700616/).
- [Clinical management, monitoring, and prophylaxis of adverse events of special interest associated with datopotamab deruxtecan. PubMed 38502995](https://pubmed.ncbi.nlm.nih.gov/38502995/).
- [US expert Delphi consensus on the prevention and management of stomatitis in patients treated with datopotamab deruxtecan. PMC12325563](https://pmc.ncbi.nlm.nih.gov/articles/PMC12325563/).
- [A real-world study of adverse event profiles associated with the four key components of antibody-drug conjugates based on the FAERS database. Front Pharmacol 2026](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2026.1702195/full).
- Shi J, et al. [Interstitial lung disease with antibody-drug conjugates: a real-world pharmacovigilance study based on the FAERS database during the period 2014-2023. Ther Adv Med Oncol. PMC11635890](https://pmc.ncbi.nlm.nih.gov/articles/PMC11635890/).
- [Systematic analysis and mechanistic investigation of cardiac adverse events associated with antibody-drug conjugates using FAERS database. PMC12825936](https://pmc.ncbi.nlm.nih.gov/articles/PMC12825936/).
- [Disproportionality analysis of interstitial lung disease associated with novel antineoplastic agents during breast cancer treatment: a pharmacovigilance study. PMC11957809](https://pmc.ncbi.nlm.nih.gov/articles/PMC11957809/).
- [Signals of interstitial lung disease with novel antineoplastic agents in ovarian cancer: a three-database disproportionality study. PMC12823856](https://pmc.ncbi.nlm.nih.gov/articles/PMC12823856/).
