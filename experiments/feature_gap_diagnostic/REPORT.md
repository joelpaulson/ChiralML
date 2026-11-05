# Feature-gap diagnostic of the completed Astra experiment

**The fixed RFs do not broadly rescue the problematic chemistries.** B32 transfers better under leave-one-chemistry-out evaluation. RFs soften some extreme mistakes and improve a few chemistries, but introduce more errors elsewhere. These results support localized limitations of B32's representation and substantial remaining transfer/information limitations; they do not establish that the descriptors are sufficient or that a better mathematical structure cannot work.

## Matched comparison

Reused Astra's exact saved partitions: 631 development observations in 47 chemistries, 47 LOCO folds, 20 repeats of five-fold polymer–solvent grouped CV, and 142 holdout observations in 12 chemistries. Molecular-weight variants remain together in LOCO. Target: `abs_g_factor >= 0.001`; errors below use the unchanged `P >= 0.5` cutoff.

B32 uses **saved fold-specific and development-only predictions**, never the all-data fit. LOCO/holdout training means, sample SDs, and predictions were reconstructed and checked. Both RFs use the original configuration: 1,100 trees, depth 12, leaf/split minimum 2, `max_features=0.5`, balanced training class weights, bootstrap sampling. Seeds follow `chiral_sr.protocol`; holdout seed is 20261104. RFs use raw descriptors and native missing-value routing, with no scaling/imputation. No equation fitting, search, feature selection, calibration, or tuning was performed.

| Evaluation | Model | AP ↑ | ROC-AUC ↑ | Log loss ↓ | Equal-chemistry loss ↓ |
|---|---|---:|---:|---:|---:|
| Development LOCO | B32 | **0.7869** | **0.8638** | **0.4472** | **0.3820** |
| | RF15 | 0.6289 | 0.6767 | 0.6754 | 0.6244 |
| | RF41 | 0.7011 | 0.7563 | 0.5960 | 0.5784 |
| Chemistry holdout | B32 | **0.6557** | **0.7594** | 0.4846 | 0.9699 |
| | RF15 | 0.5105 | 0.7238 | **0.4143** | 0.6337 |
| | RF41 | 0.5462 | 0.7168 | 0.4487 | **0.6098** |
| Pair CV, mean of 20 repeats | B32 | 0.7929 | **0.8693** | **0.4391** | — |
| | RF15 | **0.8098** | 0.8631 | 0.4575 | — |
| | RF41 | 0.8025 | 0.8665 | 0.4506 | — |

Pair-CV scores average repeat-level pooled OOF metrics. **84.7%** of its test rows share a chemistry with training; the RF AP advantage does not transfer to unseen chemistries. Compare AP within partitions because prevalence differs (development 39.8%, holdout 14.1%). Loss uses stable B32 logits and RF probabilities clipped only at machine epsilon.

RF15/RF41 improve LOCO loss in only **12/47 and 14/47** chemistries. Their paired row-weighted loss increases versus B32 are +0.228 [0.119, 0.378] and +0.149 [0.062, 0.288]. RF41 improves over RF15 by −0.079 [−0.147, −0.012]. Intervals resample entire chemistries 2,000 times with fixed predictions; they omit training and historical model-selection uncertainty.

Holdout loss gains are concentrated in IDTBT's three rows. Excluding IDTBT descriptively, losses are **0.333/0.401/0.431** for B32/RF15/RF41. Holdout intervals for RF-minus-B32 loss span zero: RF15 [−0.668, 0.310], RF41 [−0.622, 0.291]. This is reduced overconfidence on an extreme case, not a general holdout improvement.

## Chemistry errors and overlap

The [59-chemistry table](per_chemistry_comparison.csv) includes all metrics and failure overlaps; single-class AP/AUC remain undefined. Selected Astra cases and the clearest RF success follow. Entries are **log loss / errors**; H denotes holdout.

| Chemistry | n / positive | B32 | RF15 | RF41 |
|---|---:|---:|---:|---:|
| N2200, all MW variants | 102 / 62 | 0.588 / 22 | 0.697 / 51 | 0.595 / 37 |
| PII-2T (C10-C10), all MW variants | 50 / 20 | 0.748 / 19 | 0.762 / 24 | 0.912 / 26 |
| N2200-OEG (H) | 17 / 5 | 1.223 / 5 | 0.642 / 4 | 0.624 / 5 |
| IDTBT (H) | 3 / 2 | 7.528 / 2 | 1.045 / 2 | 1.269 / 2 |
| AcDOT-[OE3]2 | 3 / 0 | 1.320 / 2 | 0.667 / 2 | 0.850 / 3 |
| D18 | 5 / 5 | 1.010 / 5 | 2.349 / 5 | 2.104 / 5 |
| PE2-biOE2 | 5 / 2 | 0.781 / 3 | 0.709 / 2 | 0.613 / 2 |
| ProDOT(OE3)-co-ProDOT(EtHx) | 4 / 1 | 0.728 / 1 | 0.539 / 1 | 0.465 / 1 |
| PTQ-10 | 6 / 0 | 1.357 / 6 | 1.036 / 6 | 0.967 / 6 |
| PII-3T | 10 / 0 | 0.645 / 6 | 0.512 / 3 | 0.306 / 0 |

- **N2200:** the four positive 9.1k/anisole measurements at 100–200 mg/mL remain missed by both RFs. Probabilities rise from B32's 0.016–0.066 to RF15's 0.283–0.286 and RF41's 0.365–0.369, still below 0.5. RF41 improves over RF15 overall but remains worse than B32.
- **PII-2T:** overprediction persists: false positives increase from 16 to 20/23. Seventeen rows are wrong in all three models.
- **N2200-OEG:** RF15 rescues four of five positives, RF41 three, but each adds three false positives. Both rescue chlorobenzene at 300 mg/mL; only RF15 rescues o-dichlorobenzene there. Neither rescues the 70 °C anisole positive.
- **IDTBT:** positive probabilities rise from approximately 0.000076/0.000002 to 0.224/0.228 (RF15) and 0.158/0.161 (RF41). Both positives remain missed. D18's five positives and PTQ-10's six negatives also remain wrong in every model.
- **Smaller gains:** both RFs correct one PE2-biOE2 positive. For ProDOT, they exchange a false positive for a false negative. AcDOT-[OE3]2 has no classification gain. PII-3T is a local success on ten all-negative observations.

| Partition | B32 / RF15 / RF41 errors | Wrong in all three | RF15 corrects B32 / adds new | RF41 corrects B32 / adds new |
|---|---:|---:|---:|---:|
| Development LOCO | 138 / 227 / 175 | 82 | 43 / 132 | 40 / 77 |
| Holdout | 17 / 22 / 22 | 9 | 8 / 13 | 7 / 12 |

The RFs share **21/22 holdout errors**. See [pairwise overlap](failure_overlap.csv) and [all eight error patterns](failure_patterns.csv).

## Which descriptors plausibly contribute?

Frozen-tree path contributions allocate each parent-to-child probability change to its split feature. Contributions plus the training-root baseline reconstruct predictions within 10⁻¹². These are **model associations, not SHAP, ablation tests, or causal effects**; correlations and split order affect attribution. Numbers below average probability contributions on the specified rows, not differences from B32.

- **PII-3T:** on RF41's six corrected B32 errors, solvent logP contributes −0.170, `R_ring_2` −0.165, and backbone logP −0.121. The latter two are RF41-only descriptors. Feature-generation code defines `R_ring_2` as side-chain branch points / backbone rings. RF15's four corrected errors likewise show a large solvent-logP contribution (−0.379), a descriptor available in the 15-input pool but unused in B32. This supplies concrete evidence for useful solvent and architecture distinctions in this local success.
- **N2200-OEG and PE2-biOE2:** corrected positives receive positive contributions from concentration, side-chain unit density, and branching, already available to B32. N2200-OEG concentration contributions are +0.067/+0.084 for RF15/RF41 corrections. Alternative interactions plausibly explain these rescues. RF41-only side-chain asphericity contributes negatively here and should not simply be credited for the rescue.
- **N2200:** RF41 corrects 27 RF15 errors. Side-chain asphericity contributes +0.051 on those rows, Mn −0.119, and PDI −0.032: a changed balance of associations, not proof of any descriptor's incremental benefit.
- **Confidence-only gains:** side-chain unit density and branching strongly lower RF probabilities for AcDOT and IDTBT. ProDOT's corrected negative also receives a negative RF41 side-chain-solubility contribution (−0.055), but its positive is lost. None supports a claim of a resolved assembly mechanism.

The feature sets are **not nested**: 12 descriptors are shared; RF41 omits side-chain carbon, backbone carbon, and conjugated-unit counts and adds 29 descriptors, including family indicators. RF41−RF15 therefore does not isolate “26 extra features.” B32 uses ten of the 15 inputs.

Training support also matters. With N2200 held out, 12 observations have Mn below the training minimum of 8,333 g/mol. PII-2T includes solvent/fragment descriptors outside training ranges. D18 exceeds training maxima for monomer mass and monomer/backbone logP. Holdout IDTBT exceeds training monomer logP and its family indicator is entirely absent from training. These checks identify limited coverage, not sufficient explanations of failure. See [contributions](descriptor_contributions.csv) and RF15/RF41 descriptor-support tables.

## Outcome and replicate audit

All 773 rows have absolute g-factor labels; 271 pass the threshold, including ten exactly at 0.001. The model workbook has no signed g, measurement uncertainty, batch identifier, or preparation-history field. The available assembly labels are **c = chiral, a = aggregated, i = isotropic**, as defined in the local [Supporting Information, Section S1](../../experimental_data/Supporting_Information.docx).

There are 307 c rows (271 positive, **36 below threshold**), 279 a rows, and 187 i rows. All a/i g values are stored as zero; four c rows also have stored zero. These values alone cannot distinguish a measured exact zero from an encoding or detection-limit convention. The negative target combines weak chiral, aggregated, and isotropic outcomes. Phase labels therefore help categorize errors: all-model shared failures comprise 55 c, 15 a, and 21 i rows across the two evaluations. They are observed assembly outcomes, however, so using them as predictors would require information unavailable at the original prediction stage.

Exact recorded polymer/solvent/temperature/concentration conditions, either feature set, and even all 77 supplied features each yield **635 distinct conditions, with 14 conflicting conditions covering 49 rows** (development 13/46; holdout 1/3). Any deterministic classifier on these inputs must misclassify at least **14/773 rows**, a lower bound rather than an explanation of all errors. Every conflicting condition is phase c: conditioning on phase leaves this bound unchanged. Ten conflicts involve N2200, including a 9.1k/p-xylene/200 mg/mL condition spanning |g| = 0.0001–0.131. Disagreement is not confined to values just beside the cutoff.

The supporting information states that some 100–300 mg/mL ranges became three ML rows and that separate CD measurements were retained. The workbook's `source_record` numbering restarts between datasets: 475 unique numbers correspond to 507 `(dataset_split, source_record)` pairs. Neither these pairs nor 773 rows should be equated with independent experimental replicates. No rows were relabeled, deduplicated, or reweighted. [Condition audit](condition_audit.csv), [conflict details](conflicting_conditions.csv), and [phase audit](phase_audit.csv) preserve the evidence.

## Scope and reproducibility

RF failure cannot prove descriptor insufficiency; exact-input label conflicts do demonstrate a limit of recorded information. This is retrospective: B32's structure was discovered on all 773 rows, RF settings/features are inherited rather than selected within folds, and Astra's holdout was already examined. Small chemistry samples and correlated/expanded observations limit generalization claims.

Run `PYTHONDONTWRITEBYTECODE=1 .venv/bin/python experiments/feature_gap_diagnostic/run_diagnostic.py --summarize` to regenerate tables **without fitting**. The initial run performed 296 fixed RF fits and refuses to overwrite completed RF archives. [Provenance](provenance.json) records sources, features, settings, seeds, and versions; [verification](verification.json) records checks and output hashes. Original code and completed Astra artifacts are unchanged.
