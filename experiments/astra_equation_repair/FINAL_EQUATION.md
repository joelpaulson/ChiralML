# Final fitted equation

The selected structure is **the original five-term B32**, retained before holdout evaluation. None of the three reasoning-driven edits justified replacement. The coefficients below were subsequently refitted on all **773 observations (271 strong)**. They are for future use; the final holdout scores used the earlier development-only coefficients.

## Complete expression

\[
P(|g|\ge10^{-3}\mid x)=\sigma(f(x)),\qquad \sigma(z)=\frac{1}{1+e^{-z}},
\]

\[
\begin{aligned}
f(x)={}&-0.611886980460750\\
&+0.000920944987798687\,(S-B)cR\\
&+0.713484255385024\,HUF\\
&+0.0000799826488026122\,\frac{c^2F}{Q}\\
&-0.00206753621927866\,(S-B)UB\\
&+0.0000000502316279772739\,\frac{cMB}{E}.
\end{aligned}
\]

**Insert numerical values in the units below, without further centering or scaling.** Equivalently, each dimensional symbol denotes its numerical ratio to one of the listed units; all fitted coefficients in that convention are dimensionless. Substituting Mn in kg/mol instead of g/mol would give incorrect predictions unless its coefficient were converted.

| Symbol | Repository descriptor | Definition and numerical unit |
|---|---|---|
| c | `Concentration_mg_mL` | Polymer mass concentration, mg/mL (numerically equal to g/L and kg/m³). |
| M | `MW_Mn` | Number-average polymer molar mass, g/mol. |
| S | `SC_C_atoms_Norm_MolWt` | Side-chain carbon count per repeat unit divided by repeat-unit molar mass in kg/mol; mol/kg. |
| B | `BB_C_atoms_Norm_MolWt` | Backbone carbon count normalized in the same way; mol/kg. |
| R | `SC_n_branch_points_Norm_MolWt` | Curated side-chain branch-point count normalized in the same way; mol/kg. This is not a stereocenter count. |
| U | `SC_n_units_Norm_MolWt` | Curated number of side-chain units per repeat, normalized in the same way; mol/kg. |
| Q | `BB_conjugated_units_Norm_MolWt` | Curated number of backbone conjugated units per repeat, normalized in the same way; mol/kg. |
| E | `BB_hetero_atoms_Norm_MolWt` | Curated backbone heteroatom count per repeat, normalized in the same way; mol/kg. |
| H | `Solvent_NumHAcceptors` | RDKit solvent hydrogen-bond acceptor descriptor, dimensionless; retain the repository's mixture convention. |
| F | `Backbone_fastsolv` | Predicted log10 backbone-fragment solubility relative to 1 mol/L, dimensionless, in the specified solvent at the experimental temperature. |

Precisely, every normalized count is `1000 * count_per_repeat / MolWt_monomer`, with `MolWt_monomer` in g/mol. The analysis loader removes `_Norm_MolWt` from column names; it does **not** undo normalization. Repeat-unit masses use the uncapped repeat-unit SMILES with dummy polymerization atoms; the feature-generation notebook separately cleans/caps fragments for molecular descriptors. Use the existing curated conjugated-unit and side-chain-unit definitions rather than recounting them as rings or bonds.

For a solvent mixture, the feature-generation convention is `H = sum(w_j H_j)/sum(w_j)` and `F = log10(sum(w_j 10**F_j)/sum(w_j))`, using the supplied positive solvent mixing ratios w_j. This averages solubility in linear space, not log space, and is a descriptor convention rather than a thermodynamic mixture model. The same construction defines Fmono for the rejected proposals. No mixture coefficients were refitted here.

FASTSOLV is a small-molecule organic-solubility model. F is a calculated fragment descriptor, not an experimental polymer solubility, activity, interaction parameter, or assembly free energy. Its domain and accuracy do not establish validity for polymer phase boundaries. [FASTSOLV primary paper](https://www.nature.com/articles/s41467-025-62717-7)

## Executable transcription and numerical verification

This actual printed transcription is parsed and evaluated by `reproduce.py` against the standardized SyMANTIC implementation on every observation:

```python
f = (
    -0.611886980460750
    + 0.000920944987798687 * (S-B)*c*R
    + 0.713484255385024 * H*U*F
    + 0.0000799826488026122 * c**2*F/Q
    - 0.00206753621927866 * (S-B)*U*B
    + 0.0000000502316279772739 * c*M*B/E
)
```

`final_equation.py` provides `logit(X)` and `predict(X)` without the symbolic engine. `final_B32_fit.json` stores full-precision raw and standardized parameters, means and sample standard deviations. `verification.json` records the maximum differences and independent optimization check. Extra printed digits support reproducibility, not inferential precision.

Verified maximum error over all 773 observations: **1.33×10⁻¹⁴ in logit and 2.22×10⁻¹⁵ in probability**, including the transcription above. The independent SciPy optimizer converged and agreed in probability within **5.68×10⁻¹⁰**.

## Comparison with original B32

The mathematical structure is unchanged. Both the retained model and original B32 have five terms plus one intercept. The experiment changed coefficients through refitting, not the selected structure.

| Coefficient | Development fit (631 rows; used for holdout) | Final fit (773 rows) |
|---|---:|---:|
| Intercept | -0.206676429413 | -0.611886980461 |
| (S−B)cR | 0.00115835706274 | 0.000920944987799 |
| HUF | 0.809790263368 | 0.713484255385 |
| c²F/Q | 0.000116097591378 | 0.0000799826488026 |
| (S−B)UB | -0.00294776220184 | -0.00206753621928 |
| cMB/E | 6.29887229170e-8 | 5.02316279773e-8 |

Development-only LOCO coefficient signs were stable in all 47 folds: positive for terms 1, 2, 3 and 5; negative for term 4. The appreciable coefficient changes after adding holdout observations nevertheless show that precise magnitudes depend on the sampled chemistries. They are conditional associations among correlated descriptors.

## Mathematical consequences and physical interpretation

Write the coefficients in the displayed order as b0,…,b5. Holding polymer, solvent and temperature fixed gives

\[
f(c)=A+Lc+Kc^2,
\quad A=b_0+b_2HUF+b_4(S-B)UB,
\quad L=b_1(S-B)R+b_5MB/E,
\quad K=b_3F/Q.
\]

Thus `dP/dc = P(1−P)(L+2Kc)`. For F<0 and L>0, the probability maximum is at `c_peak = −L/(2K)`; it is a maximum of predicted strong-chirality probability, not a measured critical concentration. The equation can have zero, one or two positive crossings of P=0.5, determined by `A+Lc+Kc²=0`.

- **Molecular weight:** at fixed remaining descriptors, `∂f/∂M = b5*c*B/E > 0`. Greater Mn raises predicted probability and shifts the maximum to greater concentration when K<0. This does not imply more polymer chains per volume: their number concentration scales as c/M. Interchain connectivity or chain overlap is a speculative interpretation; the fixed ideal-coil square-root replacement was not supported by transfer performance.
- **Polymer structure:** the side-chain/backbone carbon contrast multiplies `b1*c*R + b4*U*B`. For S>B, branching increases the fitted concentration slope, while the constant structural term is negative. These terms can encode a concentration-dependent change in association with side-chain architecture. Since all counts share a mass denominator and cannot be varied independently in a real synthesis, these partial effects are not causal substitution rules.
- **Solvent:** `∂f/∂F = b2*H*U+b3*c²/Q ≥ 0`; better predicted backbone solubility raises the logit at fixed other inputs. For F<0, the direct HUF contribution is negative. Preferential solvation and competition between helical ordering and nonchiral aggregation are possible explanations, not identified mechanisms. Published experiments support distinct solvent-dependent assembly pathways, but do not validate this polynomial. [Primary solvent/assembly study](https://www.nature.com/articles/s41467-025-62915-3)
- **Temperature:** no explicit T term survives. Temperature can act only through the precomputed F descriptor. Kinetic history, dispersity, cooling/evaporation rates, orientational artifacts, and measurement uncertainty are absent.

The all-data fit has negative curvature at 733/773 row-associated conditions and negative concentration slope at 295 observed rows. For F<0 it inevitably approaches P=0 as c→∞, while its c→0 probabilities span approximately 0.000066–0.738 across the observed polymer/solvent descriptors. Neither limit is a universal physical law. The development diagnoses already exposed these problems; selection of B32 means the attempted repairs did not provide a better validated compact substitute.

Use this as an empirical classifier over the sampled domain (c=10–300 mg/mL, Mn=3,600–177,000 g/mol, T=298.15–463.15 K), not a globally valid constitutive equation. Those marginal ranges do not imply coverage of every combination. Require Q>0 and E>0. The model predicts the binary threshold event, not handedness or the magnitude of g.

## Three discriminating experiments

1. **Molecular-weight scaling within one chemistry.** For fractionated N2200 at approximately 3.6, 9.1, 22, 50 and 130 kg/mol, run matched concentration series (10–300 mg/mL) in chlorobenzene and anisole at 25 °C with matched dispersity and preparation history. Replicate near transitions. Measure scattering/viscosity overlap independently. Test whether response shifts follow M, sqrt(M), or neither; the failed A proposal does not establish M as the true scaling.
2. **Concentration turnover versus continued ordering.** Repeat the N2200 9.1 kg/mol chlorobenzene series densely between 150 and 300 mg/mL, where development LOCO curves decline despite strong high-c measurements. Pair CD/absorption with microscopy and scattering, using independent preparations and standardized thermal history. A reproducible downturn accompanied by nonchiral aggregation supports competition; a sustained high response falsifies this particular quadratic turnover in that system. Dilution/concentration cycles can distinguish equilibrium behavior from hysteresis.
3. **Absolute versus selective solvation.** Choose solvents or solvent mixtures giving similar predicted Fbb and H but differing Fmono−Fbb for PII-2T(C10-C10), and test matched c=50–250 mg/mL series. Verify actual dissolution and characterize backbone aggregation spectroscopically rather than assuming FASTSOLV predicts polymer saturation. B32 has no independent contrast effect, whereas C predicts a negative concentration-dependent contrast effect. These measurements would test that hypothesis without fitting another term to the current data.

Pre-register these comparisons and reserve new chemistries or independent batches for evaluation. The present experiment neither validates B32 prospectively nor establishes a useful reasoning-driven repair.
