"""Second and final revision: supplementary concentration/selectivity interaction."""
import experiment as e

if __name__ == '__main__':
    e.log('Interpretation of B feedback and reason for final revision',
          'Replacing absolute solubility by a contrast failed substantially: LOCO AP 0.69311 versus B32 0.78687; loss 0.53035 versus 0.44715; training loss also deteriorated to 0.47776. '
          'Thus absolute backbone solubility cannot be removed on the evidence available. The relative-solvation hypothesis is not established by the diagnostic errors. '
          'One final controlled test will preserve all five B32 terms and ask whether the contrast is complementary, rather than a replacement. No further structural trials will follow.')
    terms = list(e.laws.LAWS['B32']['terms'])
    terms.append('Concentration_mg_mL*(Monomer_fastsolv-Backbone_fastsolv)')
    e.register('C_selectivity', terms,
        'Hypothesis C: absolute backbone solubility and relative solvation encode different aspects of assembly. Keep B32 intact and add c*(Fmono-Fbb) as one sixth interaction. '
        'The existing quadratic term describes empirical high-concentration curvature; the new linear term allows the concentration slope to depend on whether the side-chain-bearing repeat unit is preferentially solvated relative to its backbone. '
        'The new term vanishes as c tends to zero and at equal predicted fragment solubilities. It changes the location of a concentration maximum without changing the original curvature sign. '
        'Its low-order form is the leading concentration-dependent correction and avoids asserting that selectivity must enter with exactly the same quadratic scaling and coefficient as absolute backbone solubility (the restriction tested in B). '
        'Expected coefficient is negative if backbone-selective association competes with helical ordering: high positive Fmono-Fbb would lower high-c PII-2T false positives. '
        'The sign remains unconstrained. Solubility contrast is still a computed proxy, not a measured transfer free energy. '
        'Accept added complexity only if primary LOCO ranking and probability quality improve with reasonable chemistry-wise consistency; retain B32 if benefits are small or inconsistent. '
        'This is proposal 3 of 3, with no exponent, scale, or regularization tuning.')
