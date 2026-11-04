"""First human/reasoning-authored structure; no search or adjustable exponent."""
import experiment as e

if __name__ == '__main__':
    e.log('B32 development diagnosis before proposal A',
          'Development fit: raw coefficients [0.00115835706274, 0.809790263368, 0.000116097591378, -0.00294776220184, 6.2988722917e-8], intercept -0.206676429413. '
          'At fixed polymer/solvent/T, f=A+Lc+Qc^2. Curvature is negative in 603/631 development rows, positive in 28; the concentration slope is negative at 244 observed rows. '
          'An interior maximum lies between 0 and 300 mg/mL for 484 row-associated conditions (not 484 independent series). '
          'A decline is not inherently impossible for strong chirality: aggregation or phase changes may suppress it. But forcing its sign through the numerical zero of log10 molar fragment solubility has no universal mechanistic justification. '
          'Largest errors include N2200 9.1k in anisole (two positive repeats at c=100 have LOCO P=0.0163), AcDOT-[C10]2-co-[OE3]2 in chloroform, and PII-2T(C10-C10) false positives. '
          'The positive c*Mn term is neither chain number concentration (c/Mn) nor a generally established overlap variable. Independent chemistry errors include PTQ-10 and D18.')
    terms = list(e.laws.LAWS['B32']['terms'])
    terms[4] = 'Concentration_mg_mL*sqrt(MW_Mn/1000)*(BB_C_atoms/BB_hetero_atoms)'
    e.register('A_overlap', terms,
        'Hypothesis A: B32 overstates the molecular-weight dependence of concentration-induced ordering. Replace only c*M*B/E with c*sqrt(M/M0)*B/E, M0=1000 g/mol; retain the other four terms. '
        'For an ideal coil with size Rg proportional to M^(1/2), overlap concentration c* scales as M/Rg^3 proportional to M^(-1/2); c/c* therefore scales as c*sqrt(M). '
        'This motivates a fixed exponent 1/2, not an exponent fitted or searched from these outcomes. B/E remains the original empirical backbone-composition modifier. '
        'These conjugated chains can be semiflexible and concentrations are high, so the overlap analogy is explicitly speculative. '
        'Expected: weaker increase with molecular weight, less disparity between low- and high-M samples, and improved low-M N2200 transfer without adding a coefficient. '
        'Concentration curvature and the nonzero dilute-limit probability remain unresolved. No monotonicity or coefficient-sign constraint is imposed; all coefficients are refitted with identical ridge penalty. '
        'Reject if improvement is confined to training or a few heavily sampled chemistries.')
