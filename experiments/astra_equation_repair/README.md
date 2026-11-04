# Reasoning-driven equation repair

**Result: retain B32.** Three reasoning-driven structural edits did not provide a convincing improvement in transfer to unseen chemistries. All proposals, including unsuccessful ones, are preserved.

- [Complete study in one reading file](CHATGPT_READING.md): report, final equation, chronological log, numerical verification, and an artifact index.
- [Scientific report](REPORT.md), [complete fitted equation and experimental tests](FINAL_EQUATION.md), and [chronological log](EXPERIMENT_LOG.md).
- [Development concentration figure](figures/development_concentration.png) and [final holdout figure](figures/holdout_concentration.png).

For ChatGPT or another web-enabled reader, provide the [plain-text reading URL](https://raw.githubusercontent.com/joelpaulson/ChiralML/astra-equation-repair/experiments/astra_equation_repair/CHATGPT_READING.md). It contains the study's conclusions, equations, results, limitations, and proposal history in one document, with public links to every artifact. If a chat cannot open web links, download this Markdown file and attach it. Individual reports remain available separately.

All experiment code and artifacts live here. The existing `chiral_sr` implementation and workbook were read, not modified. The original experiment was completed locally before publication was requested. Python uses the repository `.venv`; exact package versions, source hashes, the workbook SHA-256, grouping columns and split seed are in `split.json` and `environment.json`.

## Run from a fresh clone

```bash
git clone --branch astra-equation-repair https://github.com/joelpaulson/ChiralML.git
cd ChiralML
python3.12 -m venv .venv
.venv/bin/python -m pip install -r experiments/astra_equation_repair/requirements.txt
.venv/bin/python experiments/astra_equation_repair/reproduce.py
```

The branch includes the source workbook and the existing symbolic-regression implementation required by the experiment. The pinned requirements record the environment used; CPU PyTorch is sufficient.

## Verify completed results

From the repository root:

```bash
.venv/bin/python experiments/astra_equation_repair/reproduce.py
.venv/bin/python experiments/astra_equation_repair/diagnostics.py
.venv/bin/python experiments/astra_equation_repair/post_evaluation.py
.venv/bin/python experiments/astra_equation_repair/figures.py
```

`reproduce.py` performs a fixed numerical replay, not another model-selection exercise. It reproduces the label-independent split and all 147 training folds per model (47 LOCO plus 100 polymer–solvent folds), checks disjoint groups, recomputes final all-data coefficients, independently solves the same objective, and compares the printed equation with implemented predictions. It verifies holdout metrics from the saved single-evaluation logits/labels. It overwrites only `verification.json`. The other three commands regenerate numerical diagnostic summaries and two figures without changing structures or fits.

Reports and chronology are deliberately authored interpretations, not generated selection decisions. The per-proposal registrars contain exactly the equations and rationales that were logged before evaluation. No automated candidate-generation code is present.

## Original staged sequence

These commands document how the study was run. **Do not rerun setup/final over this completed record**: they intentionally refuse to redraw a split or reopen the one-shot holdout. To conduct a clean replication from the workbook, copy these scripts into an equivalently located experiment folder, keep the repository import paths valid, and execute the following there. The retained B32 decision is a documented development judgment, not a general automatic selection algorithm.

```bash
python experiment.py setup
python experiment.py evaluate --model B32
python proposal_A.py
python experiment.py evaluate --model A_overlap
python proposal_B.py
python experiment.py evaluate --model B_contrast
python proposal_C.py
python experiment.py evaluate --model C_selectivity
python diagnostics.py
python select_structure.py
python experiment.py final
python post_evaluation.py
python figures.py
```

`setup` reads only features and identifier metadata. `load_partition('development')` excludes held-aside rows before parsing metadata outcomes. `final` requires `selection.json`, checks the immutable registry hash, and writes `final_started.json` before loading holdout outcomes. Only B32 is evaluated there because it is also the selected structure. Full-data refitting occurs after that assessment. Raw probability calculations use SciPy's stable sigmoid and log loss uses logits through `logaddexp`, avoiding clipping extreme mistakes out of the score.

## Artifacts

- `EXPERIMENT_LOG.md`: timestamped protocol, proposals before feedback, all results, reasoning and frozen decision.
- `REPORT.md`, `FINAL_EQUATION.md`: interpretation and complete equation.
- `split.json`, `split_manifest.csv`, `folds.json`: source provenance, row/group assignments and identical development partitions. Fold indices are **local to development row order**, with original row indices in the manifest/prediction archives.
- `registry.json`, `selection.json`, `final_started.json`: exact four structures, frozen selection and holdout-opening record.
- `results/<model>/`: raw/standardized development coefficients, all LOCO coefficients, OOF logits (not probabilities) in NPZ, repeat-level and chemistry-level metrics. `pair` arrays have shape 20×631; `loco` has length 631. `y`, `chem` and original zero-based workbook rows accompany them.
- `development_comparison.csv`, `development_predictions.csv`, `development_paired_comparisons.json`: full score comparison, labeled development predictions and paired bootstrap uncertainty.
- `holdout_metrics.json`, `holdout_predictions.npz`, `holdout_predictions.csv`, `holdout_B32_per_chemistry.csv`: the single frozen assessment. The NPZ `B32` field holds logits.
- `final_B32_fit.json`, `final_equation.py`, `verification.json`: final 773-row fitted equation and numerical audit.
- `figures/`: two concentration plots, each in PNG and vector PDF, plus exact plotted-row provenance.

`final_equation.predict(X)` expects the unsuffixed, already mass-normalized columns from the analysis loader; it applies no additional scaling. Original Excel count column names end in `_Norm_MolWt`. Do not pass unnormalized atom counts.

The bootstrap resamples chemistry groups, seed 20261105, 2,000 draws. It does not refit models or account for adaptive structure selection. Single-class per-chemistry AP/AUC are represented as blank CSV cells / null, not scored as zero. The substantial chemical and measurement dependencies are discussed in the report.
