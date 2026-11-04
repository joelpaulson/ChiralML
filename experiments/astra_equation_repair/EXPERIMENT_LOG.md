# Chronological equation-repair experiment

Date: 2026-11-04 (America/Chicago); entries use UTC. No outcome data were inspected before the split.

## 2026-11-04T16:01:51.305778+00:00 — Protocol and partition fixed before outcome inspection

Seed 20261104; 59 chemistry groups, 12 held aside. Development n=631; holdout n=142. The manifest contains identifiers only. Group definitions and all defining columns are in split.json. All molecular-weight variants sharing the polymer descriptors stay together. No stratification or redraw.

Primary: leave-one-chemistry-out (LOCO), pooled AP, ROC-AUC and stable log loss, plus equal-chemistry loss and mixed-class within-chemistry metrics. Complement: existing 20 repeats of 5-fold StratifiedGroupKFold, polymer-name || solvent, seeds 1000–1019. All models use identical saved folds, train-only term standardization (sample SD), summed logistic loss + 0.05 sum(beta_standardized^2), unpenalized intercept, SyMANTIC fit_glm, no chemistry offsets, no hyperparameter tuning. At most three new structures, each registered before evaluation. Prefer a five-term structure unless a sixth term brings consistent LOCO benefit. Assess transfer using AP together with log loss and chemistry-wise consistency; do not select on training scores or complementary CV alone. Bootstrap intervals resample chemistries of fixed OOF predictions, so omit selection and refitting uncertainty. Freeze selection before final labels are read. Final evaluates only selected and B32 once, then refits selected on all data. B32 was originally discovered on the full dataset: this is not a pristine prospective holdout.

Original B32: f=b0+b1*(S-B)*c*R+b2*H*U*F+b3*c^2*F/Q+b4*(S-B)*U*B+b5*c*M*B/E. S=SC_C_atoms, B=BB_C_atoms, R=SC_n_branch_points, U=SC_n_units, H=Solvent_NumHAcceptors, F=Backbone_fastsolv, Q=BB_conjugated_units, E=BB_hetero_atoms, M=MW_Mn, c=Concentration_mg_mL. All count descriptors are normalized mol/kg. This is the original structure with numerical coefficients deliberately refitted within each training partition.

## 2026-11-04T16:02:14.468433+00:00 — Feedback for B32

Development only.

```json
{
  "model": "B32",
  "n": 631,
  "positives": 251,
  "terms": 5,
  "loco": {
    "ap": 0.7868677086620305,
    "auc": 0.8638341371356678,
    "log_loss": 0.4471500442559523,
    "brier": 0.14695590648658535
  },
  "chemistry_macro": {
    "ap": 0.8654290362549262,
    "auc": 0.8608862073049255,
    "log_loss": 0.38198710972724254,
    "brier": 0.12602683719200547
  },
  "mixed_class_chemistries": 26,
  "pair_mean": {
    "ap": 0.7929411928073359,
    "auc": 0.8693143216607254,
    "log_loss": 0.43906587302607464,
    "brier": 0.1433871340060193
  },
  "pair_sd": {
    "ap": 0.01096930370605735,
    "auc": 0.0062924974135190265,
    "log_loss": 0.01149545694711452,
    "brier": 0.003907558452565901
  },
  "training": {
    "ap": 0.8192004737789333,
    "auc": 0.8865642692388341,
    "log_loss": 0.4074076906959214,
    "brier": 0.1337202831043698
  },
  "evaluated_utc": "2026-11-04T16:02:14.466503+00:00"
}
```

Per-chemistry results, all OOF predictions, and fitted coefficients saved under results/B32/.

## 2026-11-04T16:03:58.388986+00:00 — B32 development diagnosis before proposal A

Development fit: raw coefficients [0.00115835706274, 0.809790263368, 0.000116097591378, -0.00294776220184, 6.2988722917e-8], intercept -0.206676429413. At fixed polymer/solvent/T, f=A+Lc+Qc^2. Curvature is negative in 603/631 development rows, positive in 28; the concentration slope is negative at 244 observed rows. An interior maximum lies between 0 and 300 mg/mL for 484 row-associated conditions (not 484 independent series). A decline is not inherently impossible for strong chirality: aggregation or phase changes may suppress it. But forcing its sign through the numerical zero of log10 molar fragment solubility has no universal mechanistic justification. Largest errors include N2200 9.1k in anisole (two positive repeats at c=100 have LOCO P=0.0163), AcDOT-[C10]2-co-[OE3]2 in chloroform, and PII-2T(C10-C10) false positives. The positive c*Mn term is neither chain number concentration (c/Mn) nor a generally established overlap variable. Independent chemistry errors include PTQ-10 and D18.

## 2026-11-04T16:03:58.389411+00:00 — PROPOSAL A_overlap — before any fit or evaluation

Hypothesis A: B32 overstates the molecular-weight dependence of concentration-induced ordering. Replace only c*M*B/E with c*sqrt(M/M0)*B/E, M0=1000 g/mol; retain the other four terms. For an ideal coil with size Rg proportional to M^(1/2), overlap concentration c* scales as M/Rg^3 proportional to M^(-1/2); c/c* therefore scales as c*sqrt(M). This motivates a fixed exponent 1/2, not an exponent fitted or searched from these outcomes. B/E remains the original empirical backbone-composition modifier. These conjugated chains can be semiflexible and concentrations are high, so the overlap analogy is explicitly speculative. Expected: weaker increase with molecular weight, less disparity between low- and high-M samples, and improved low-M N2200 transfer without adding a coefficient. Concentration curvature and the nonzero dilute-limit probability remain unresolved. No monotonicity or coefficient-sign constraint is imposed; all coefficients are refitted with identical ridge penalty. Reject if improvement is confined to training or a few heavily sampled chemistries.

Exact terms (with free intercept and one fitted coefficient each):

- `((SC_C_atoms-BB_C_atoms)*(Concentration_mg_mL*SC_n_branch_points))`
- `(Solvent_NumHAcceptors*(SC_n_units*Backbone_fastsolv))`
- `((Concentration_mg_mL*Backbone_fastsolv)*(Concentration_mg_mL/BB_conjugated_units))`
- `((SC_C_atoms-BB_C_atoms)*(SC_n_units*BB_C_atoms))`
- `Concentration_mg_mL*sqrt(MW_Mn/1000)*(BB_C_atoms/BB_hetero_atoms)`

## 2026-11-04T16:04:01.720475+00:00 — Feedback for A_overlap

Development only.

```json
{
  "model": "A_overlap",
  "n": 631,
  "positives": 251,
  "terms": 5,
  "loco": {
    "ap": 0.7777905365928819,
    "auc": 0.8609509331096666,
    "log_loss": 0.4545791046037808,
    "brier": 0.14917692202589075
  },
  "chemistry_macro": {
    "ap": 0.8663321152407687,
    "auc": 0.8611416839071762,
    "log_loss": 0.38427138942465705,
    "brier": 0.12801021411645308
  },
  "mixed_class_chemistries": 26,
  "pair_mean": {
    "ap": 0.7914802415331047,
    "auc": 0.8705708743971481,
    "log_loss": 0.4402376683402035,
    "brier": 0.14345967677622118
  },
  "pair_sd": {
    "ap": 0.011163404360936139,
    "auc": 0.006010775898240957,
    "log_loss": 0.010551598141432061,
    "brier": 0.0036498121435323363
  },
  "training": {
    "ap": 0.8203401992569661,
    "auc": 0.8887659886768714,
    "log_loss": 0.40819918072097,
    "brier": 0.13370630467113742
  },
  "evaluated_utc": "2026-11-04T16:04:01.719108+00:00"
}
```

Per-chemistry results, all OOF predictions, and fitted coefficients saved under results/A_overlap/.

## 2026-11-04T16:05:23.255719+00:00 — Interpretation of A feedback and reason to revise

A worsened LOCO AP by 0.00908 and loss by 0.00743. Chemistry loss improved for 25/47 but macro loss worsened (0.38427 versus 0.38199). It did not repair the targeted N2200-anisole errors (P at c=100 fell from 0.01625 to 0.01392); N2200 chemistry loss worsened by 0.01378. PII-2T(C10-C10) loss worsened by 0.06108. Thus the ideal-coil analogy is not supported as a useful repair here. Restore original Mn term. The largest opposing errors also differ in full-repeat-unit versus backbone solubility: high-confidence PII-2T false positives have Fmono-Fbb around 0.78–1.11, whereas missed N2200-anisole positives have gap 0.23. This suggests testing differential solvation, not further adjusting the MW exponent.

## 2026-11-04T16:05:23.256343+00:00 — PROPOSAL B_contrast — before any fit or evaluation

Hypothesis B: concentration-induced suppression of helical order depends on the solvent selectivity between backbone and full repeat unit, not solely absolute fragment solubility. Replace only c^2*Fbb/Q with c^2*(Fbb-Fmono)/Q, and restore B32's original c*Mn*B/E term after A failed. Five terms, no added coefficients. Fbb-Fmono = log10(Sbb/Smono) is a dimensionless log ratio with a meaningful zero (equal predicted fragment solubilities), unlike the original curvature zero at 1 mol/L. For Fbb < Fmono, side chains may sustain dissolution while the backbone preferentially associates; overly strong backbone association could compete with torsional/helical assembly. A positive fitted coefficient would encode increasing suppression with concentration for that situation. This sign is an expectation, not a constraint. Expected: stronger high-c suppression for large positive Fmono-Fbb, less suppression for small contrasts; help PII-2T false positives while sparing N2200. The denominator Q and quadratic concentration remain inherited empirical choices; no microscopic free-energy derivation is asserted. FASTSOLV predicts small-molecule fragment solubilities, not polymer saturation concentrations; this difference is a proxy, not a measured solvent selectivity or Flory-Huggins parameter. A clear failure would indicate that absolute backbone solubility carries useful information lost by a pure contrast.

Exact terms (with free intercept and one fitted coefficient each):

- `((SC_C_atoms-BB_C_atoms)*(Concentration_mg_mL*SC_n_branch_points))`
- `(Solvent_NumHAcceptors*(SC_n_units*Backbone_fastsolv))`
- `(Concentration_mg_mL**2*(Backbone_fastsolv-Monomer_fastsolv)/BB_conjugated_units)`
- `((SC_C_atoms-BB_C_atoms)*(SC_n_units*BB_C_atoms))`
- `((Concentration_mg_mL*MW_Mn)*(BB_C_atoms/BB_hetero_atoms))`

## 2026-11-04T16:05:26.166258+00:00 — Feedback for B_contrast

Development only.

```json
{
  "model": "B_contrast",
  "n": 631,
  "positives": 251,
  "terms": 5,
  "loco": {
    "ap": 0.6931102848371479,
    "auc": 0.8000681484587965,
    "log_loss": 0.5303537374163343,
    "brier": 0.18009401728177207
  },
  "chemistry_macro": {
    "ap": 0.8367390577994485,
    "auc": 0.8356638003704359,
    "log_loss": 0.43878956398378405,
    "brier": 0.14644931049174953
  },
  "mixed_class_chemistries": 26,
  "pair_mean": {
    "ap": 0.7199830127835621,
    "auc": 0.8179503040469701,
    "log_loss": 0.5064352697027765,
    "brier": 0.17024824657134624
  },
  "pair_sd": {
    "ap": 0.011056714658391142,
    "auc": 0.005967131576159901,
    "log_loss": 0.009045627867002448,
    "brier": 0.0029816567250762357
  },
  "training": {
    "ap": 0.750529036633339,
    "auc": 0.8377175508492346,
    "log_loss": 0.47776162340023415,
    "brier": 0.16048655412751295
  },
  "evaluated_utc": "2026-11-04T16:05:26.164829+00:00"
}
```

Per-chemistry results, all OOF predictions, and fitted coefficients saved under results/B_contrast/.

## 2026-11-04T16:06:15.491665+00:00 — Interpretation of B feedback and reason for final revision

Replacing absolute solubility by a contrast failed substantially: LOCO AP 0.69311 versus B32 0.78687; loss 0.53035 versus 0.44715; training loss also deteriorated to 0.47776. Thus absolute backbone solubility cannot be removed on the evidence available. The relative-solvation hypothesis is not established by the diagnostic errors. One final controlled test will preserve all five B32 terms and ask whether the contrast is complementary, rather than a replacement. No further structural trials will follow.

## 2026-11-04T16:06:15.492810+00:00 — PROPOSAL C_selectivity — before any fit or evaluation

Hypothesis C: absolute backbone solubility and relative solvation encode different aspects of assembly. Keep B32 intact and add c*(Fmono-Fbb) as one sixth interaction. The existing quadratic term describes empirical high-concentration curvature; the new linear term allows the concentration slope to depend on whether the side-chain-bearing repeat unit is preferentially solvated relative to its backbone. The new term vanishes as c tends to zero and at equal predicted fragment solubilities. It changes the location of a concentration maximum without changing the original curvature sign. Its low-order form is the leading concentration-dependent correction and avoids asserting that selectivity must enter with exactly the same quadratic scaling and coefficient as absolute backbone solubility (the restriction tested in B). Expected coefficient is negative if backbone-selective association competes with helical ordering: high positive Fmono-Fbb would lower high-c PII-2T false positives. The sign remains unconstrained. Solubility contrast is still a computed proxy, not a measured transfer free energy. Accept added complexity only if primary LOCO ranking and probability quality improve with reasonable chemistry-wise consistency; retain B32 if benefits are small or inconsistent. This is proposal 3 of 3, with no exponent, scale, or regularization tuning.

Exact terms (with free intercept and one fitted coefficient each):

- `((SC_C_atoms-BB_C_atoms)*(Concentration_mg_mL*SC_n_branch_points))`
- `(Solvent_NumHAcceptors*(SC_n_units*Backbone_fastsolv))`
- `((Concentration_mg_mL*Backbone_fastsolv)*(Concentration_mg_mL/BB_conjugated_units))`
- `((SC_C_atoms-BB_C_atoms)*(SC_n_units*BB_C_atoms))`
- `((Concentration_mg_mL*MW_Mn)*(BB_C_atoms/BB_hetero_atoms))`
- `Concentration_mg_mL*(Monomer_fastsolv-Backbone_fastsolv)`

## 2026-11-04T16:06:18.774424+00:00 — Feedback for C_selectivity

Development only.

```json
{
  "model": "C_selectivity",
  "n": 631,
  "positives": 251,
  "terms": 6,
  "loco": {
    "ap": 0.7882538973598187,
    "auc": 0.8652180750681485,
    "log_loss": 0.44876797241226984,
    "brier": 0.14691784793122872
  },
  "chemistry_macro": {
    "ap": 0.865389022603151,
    "auc": 0.8603763882709129,
    "log_loss": 0.39448323889479636,
    "brier": 0.12790050591197266
  },
  "mixed_class_chemistries": 26,
  "pair_mean": {
    "ap": 0.7941844740169923,
    "auc": 0.8710846089326904,
    "log_loss": 0.44129556100178596,
    "brier": 0.14286985527092586
  },
  "pair_sd": {
    "ap": 0.011762700874849085,
    "auc": 0.006594140505835352,
    "log_loss": 0.012627094951507429,
    "brier": 0.004187042396680783
  },
  "training": {
    "ap": 0.8251020897852349,
    "auc": 0.891460473893898,
    "log_loss": 0.40305031761801224,
    "brier": 0.1311743231757315
  },
  "evaluated_utc": "2026-11-04T16:06:18.772707+00:00"
}
```

Per-chemistry results, all OOF predictions, and fitted coefficients saved under results/C_selectivity/.

## 2026-11-04T16:07:43.900531+00:00 — Literature-search exposure caveat

A broad background search after the split returned a published numerical concentration/chirality example for DPP-T2M2 (a held-aside chemistry) in a search snippet. No workbook holdout outcomes, predictions, scores or errors were loaded, and that example was not used as candidate feedback or a selection criterion. Nevertheless, blindness to all external published information was not perfect. This is an additional exploratory-study limitation, beyond B32 having been discovered on the full dataset. Subsequent reasoning and every numerical comparison used development data only. Background source: https://pmc.ncbi.nlm.nih.gov/articles/PMC10683494/ .

## 2026-11-04T16:08:53.337121+00:00 — Final development feedback and selection reasoning

Retain five-term B32. A failed its targeted low-M repair and worsened primary scores; B lost substantial predictive information. C has the expected negative selectivity coefficient in all 47 LOCO fits and improves loss in 32/47 chemistries (median delta -0.00659), which makes the hypothesis scientifically interesting. However, its pooled AP gain is only 0.001386 (paired chemistry-bootstrap 95% interval [-0.02208, 0.01812]); LOCO loss worsens by 0.001618 (interval [-0.01395, 0.01839]), macro chemistry loss worsens from 0.38199 to 0.39448, and complementary grouped-CV loss worsens in 16/20 repeats. The losses for AcDOT-[OE3]2, D18, and several DOT chemistries outweigh many small gains. Training improvement and a stable coefficient sign do not establish transfer benefit. These small and inconsistent gains do not justify a sixth term. This is a negative equation-repair result, not evidence that B32 is a physical law. All three allowed new structures have been used; no new candidates or structure decisions will follow the holdout.

## 2026-11-04T16:08:53.340201+00:00 — STRUCTURE FROZEN before holdout labels or predictions

Selected B32. Retain five-term B32. A failed its targeted low-M repair and worsened primary scores; B lost substantial predictive information. C has the expected negative selectivity coefficient in all 47 LOCO fits and improves loss in 32/47 chemistries (median delta -0.00659), which makes the hypothesis scientifically interesting. However, its pooled AP gain is only 0.001386 (paired chemistry-bootstrap 95% interval [-0.02208, 0.01812]); LOCO loss worsens by 0.001618 (interval [-0.01395, 0.01839]), macro chemistry loss worsens from 0.38199 to 0.39448, and complementary grouped-CV loss worsens in 16/20 repeats. The losses for AcDOT-[OE3]2, D18, and several DOT chemistries outweigh many small gains. Training improvement and a stable coefficient sign do not establish transfer benefit. These small and inconsistent gains do not justify a sixth term. This is a negative equation-repair result, not evidence that B32 is a physical law. All three allowed new structures have been used; no new candidates or structure decisions will follow the holdout.

No further modifications or changes to selection are permitted. Registry and fold hashes stored in selection.json.

## 2026-11-04T16:08:54.866161+00:00 — Final holdout opened

Frozen model and B32 only; coefficients previously fitted on development. This is the first outcome read for held-aside rows.

## 2026-11-04T16:08:56.139389+00:00 — Final holdout results (selection unchanged)

```json
{
  "B32": {
    "ap": 0.6557159672063086,
    "auc": 0.7594262295081967,
    "log_loss": 0.4845644413971913,
    "brier": 0.1033133346963426
  }
}
```

## 2026-11-04T16:22:35.647265+00:00 — Completed all-data refit, figures and numerical audit

Preferred structure remains B32. All 773 rows (271 positives) were used only after the frozen holdout evaluation to fit the final deployable equation. Raw intercept -0.6118869804607501; raw coefficients [0.0009209449877986874, 0.7134842553850236, 0.00007998264880261217, -0.002067536219278662, 5.0231627977273944e-8].

All 588 development-fold fits reproduced exactly, with no structure changes. The final all-data fit reproduced exactly. The actual printed Markdown equation and independent direct implementation agree with the fitted model within 2.220446049250313e-15 in probability. Independent SciPy BFGS, starting at zero with the identical ridge objective, converged and agreed within 5.677256442737644e-10. Saved holdout scores were arithmetically verified from the original predictions; rejected models were never tested on holdout.

The first independent-optimizer check emitted NumPy matmul runtime warnings despite finite matching results. Explicit summation in that independent check removed the warnings; the optimizer then converged without warnings. This was implementation debugging, not a structure change. Both final figures were rendered and visually inspected. No existing SyMANTIC files or workbook were modified and no remote push was performed. Final scientific conclusion: no convincing improvement from the three attempted reasoning-driven structural edits.
