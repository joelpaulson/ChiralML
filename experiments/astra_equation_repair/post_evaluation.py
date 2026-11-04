"""Describe the already frozen predictions; never select or alter a structure."""
import json
import numpy as np
import pandas as pd
from scipy.special import expit
import experiment as e


def main():
    assert (e.HERE / 'holdout_metrics.json').exists()
    h = np.load(e.HERE / 'holdout_predictions.npz')
    y, eta = h['y'], h['B32']
    w = e.protocol.cluster_weights(h['chem'], 2000, 20261105)
    # Explicit sums avoid unnecessary platform BLAS dispatch on short matrices.
    valid = w[((w*y).sum(1)>0)&((w*(1-y)).sum(1)>0)]
    vals = [e.metrics(y, eta, v) for v in valid]
    ci = {k: np.quantile([v[k] for v in vals], [.025, .975]).tolist() for k in vals[0]}
    per = pd.read_csv(e.HERE / 'holdout_B32_per_chemistry.csv')
    out = dict(ci95=ci, bootstrap_valid=len(valid), bootstrap_total=2000,
               macro_log_loss=float(per.log_loss.mean()), true_positive_0_5=int(((expit(eta)>=.5)&(y==1)).sum()),
               false_positive_0_5=int(((expit(eta)>=.5)&(y==0)).sum()))
    e.save(e.HERE / 'holdout_diagnostics.json', out)
    print(json.dumps(out, indent=2))
    d = e.load_partition('holdout')
    err = d.meta.copy()
    err['chemistry'] = d.chem
    err['p_frozen'] = expit(eta)
    err['log_loss'] = np.logaddexp(0, eta)-y*eta
    err.to_csv(e.HERE / 'holdout_predictions.csv', index=False)
    d = e.load_partition('all')
    fit = e.unpack(e.read(e.HERE / 'final_B32_fit.json'))
    b = np.array(e.pack(fit)['raw_coefficients'])
    x = d.X
    L = b[0]*(x.SC_C_atoms-x.BB_C_atoms)*x.SC_n_branch_points+b[4]*x.MW_Mn*x.BB_C_atoms/x.BB_hetero_atoms
    Q = b[2]*x.Backbone_fastsolv/x.BB_conjugated_units
    A = e.pack(fit)['raw_intercept']+b[1]*x.Solvent_NumHAcceptors*x.SC_n_units*x.Backbone_fastsolv+b[3]*(x.SC_C_atoms-x.BB_C_atoms)*x.SC_n_units*x.BB_C_atoms
    out = dict(all_n=len(d.y), all_positive=int(d.y.sum()),
               training=e.metrics(d.y, fit.logit(e.laws.term_values(fit.terms, x))),
               negative_curvature=int((Q<0).sum()), negative_slope_observed=int(((L+2*Q*x[e.laws.C])<0).sum()),
               p_at_zero_range=[float(expit(A).min()), float(expit(A).max())],
               ranges={c: [float(x[c].min()), float(x[c].max())] for c in x})
    e.save(e.HERE / 'final_diagnostics.json', out)
    co = np.array([p['raw_coefficients'] for p in e.read(e.HERE / 'results/B32/loco_fits.json')])
    e.save(e.HERE / 'coefficient_stability.json', dict(loco_raw_min=co.min(0).tolist(), loco_raw_max=co.max(0).tolist(), positive_counts=(co>0).sum(0).tolist()))


if __name__ == '__main__':
    main()
