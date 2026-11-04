"""First revision: replace the curvature descriptor with a solubility contrast."""
import experiment as e

if __name__ == '__main__':
    e.log('Interpretation of A feedback and reason to revise',
          'A worsened LOCO AP by 0.00908 and loss by 0.00743. Chemistry loss improved for 25/47 but macro loss worsened (0.38427 versus 0.38199). '
          'It did not repair the targeted N2200-anisole errors (P at c=100 fell from 0.01625 to 0.01392); N2200 chemistry loss worsened by 0.01378. '
          'PII-2T(C10-C10) loss worsened by 0.06108. Thus the ideal-coil analogy is not supported as a useful repair here. Restore original Mn term. '
          'The largest opposing errors also differ in full-repeat-unit versus backbone solubility: high-confidence PII-2T false positives have Fmono-Fbb around 0.78–1.11, '
          'whereas missed N2200-anisole positives have gap 0.23. This suggests testing differential solvation, not further adjusting the MW exponent.')
    terms = list(e.laws.LAWS['B32']['terms'])
    terms[2] = '(Concentration_mg_mL**2*(Backbone_fastsolv-Monomer_fastsolv)/BB_conjugated_units)'
    e.register('B_contrast', terms,
        'Hypothesis B: concentration-induced suppression of helical order depends on the solvent selectivity between backbone and full repeat unit, not solely absolute fragment solubility. '
        'Replace only c^2*Fbb/Q with c^2*(Fbb-Fmono)/Q, and restore B32\'s original c*Mn*B/E term after A failed. Five terms, no added coefficients. '
        'Fbb-Fmono = log10(Sbb/Smono) is a dimensionless log ratio with a meaningful zero (equal predicted fragment solubilities), unlike the original curvature zero at 1 mol/L. '
        'For Fbb < Fmono, side chains may sustain dissolution while the backbone preferentially associates; overly strong backbone association could compete with torsional/helical assembly. '
        'A positive fitted coefficient would encode increasing suppression with concentration for that situation. This sign is an expectation, not a constraint. '
        'Expected: stronger high-c suppression for large positive Fmono-Fbb, less suppression for small contrasts; help PII-2T false positives while sparing N2200. '
        'The denominator Q and quadratic concentration remain inherited empirical choices; no microscopic free-energy derivation is asserted. '
        'FASTSOLV predicts small-molecule fragment solubilities, not polymer saturation concentrations; this difference is a proxy, not a measured solvent selectivity or Flory-Huggins parameter. '
        'A clear failure would indicate that absolute backbone solubility carries useful information lost by a pure contrast.')
