"""Development-only comparisons; no candidate generation or holdout access."""
import json
import numpy as np
import pandas as pd
from scipy.special import expit
import experiment as e


def paired_intervals(y, chem, a, b, n_boot=2000):
    W = e.protocol.cluster_weights(chem, n_boot=n_boot, seed=20261105)
    vals = []
    for w in W:
        if 0 < np.dot(w, y) < w.sum():
            ma, mb = e.metrics(y, a, w), e.metrics(y, b, w)
            vals.append([ma[k]-mb[k] for k in ['ap', 'auc', 'log_loss', 'brier']])
    return {k: np.quantile(np.array(vals)[:, j], [.025, .975]).tolist()
            for j, k in enumerate(['ap', 'auc', 'log_loss', 'brier'])}


def main():
    d = e.load_partition()
    registry = e.read(e.HERE / 'registry.json')
    baseline = np.load(e.HERE / 'results/B32/predictions.npz')
    base_per = pd.read_csv(e.HERE / 'results/B32/per_chemistry.csv')
    baseline_fit = e.read(e.HERE / 'results/B32/development_fit.json')
    table, comparisons, slopes, errors = [], {}, {}, []
    for name in registry:
        out = e.HERE / 'results' / name
        m = e.read(out / 'metrics.json')
        pred = np.load(out / 'predictions.npz')
        per = pd.read_csv(out / 'per_chemistry.csv')
        table.append(dict(model=name, terms=m['terms'], **{f'loco_{k}': v for k, v in m['loco'].items()},
                          **{f'pair_{k}': v for k, v in m['pair_mean'].items()}, macro_log_loss=m['chemistry_macro']['log_loss']))
        if name != 'B32':
            delta = per.log_loss-base_per.log_loss
            pair_base = pd.read_csv(e.HERE / 'results/B32/pair_repeat_metrics.csv')
            pair_new = pd.read_csv(out / 'pair_repeat_metrics.csv')
            delta_table = per[['chemistry', 'polymers', 'n', 'positives']].copy()
            for k in ['ap', 'auc', 'log_loss', 'brier']:
                delta_table[f'delta_{k}'] = per[k]-base_per[k]
            delta_table.to_csv(out / 'paired_chemistry_deltas.csv', index=False)
            comparisons[name] = dict(delta={k: m['loco'][k]-e.metrics(d.y, baseline['loco'])[k] for k in m['loco']},
                                     ci95=paired_intervals(d.y, d.chem, pred['loco'], baseline['loco']),
                                     chemistry_loss_wins=int((delta < 0).sum()), chemistry_loss_losses=int((delta > 0).sum()),
                                     median_chemistry_loss_delta=float(delta.median()),
                                     pair_AP_wins=int((pair_new.ap > pair_base.ap).sum()),
                                     pair_loss_wins=int((pair_new.log_loss < pair_base.log_loss).sum()))
        fit = e.unpack(e.read(out / 'development_fit.json'))
        x0, x1, x2 = (d.X.assign(Concentration_mg_mL=c) for c in [0., 1., 2.])
        f0, f1, f2 = (fit.logit(e.laws.term_values(fit.terms, x)) for x in [x0, x1, x2])
        q = (f2-2*f1+f0)/2
        lin = f1-f0-q
        slope = lin+2*q*d.X[e.laws.C].to_numpy()
        turns = np.divide(-lin, 2*q, out=np.full(len(q), np.nan), where=np.abs(q)>1e-12)
        slopes[name] = dict(decreasing_at_observation=int((slope<0).sum()), negative_curvature=int((q<0).sum()),
                            turning_0_300=int(((turns>0)&(turns<300)).sum()),
                            zero_concentration_probability_range=[float(expit(f0).min()), float(expit(f0).max())])
        err = d.meta[['experiment_id', 'polymer', 'solvent', 'abs_g_factor', 'family']].copy()
        err['model'] = name
        err['chemistry'] = d.chem
        err['concentration'] = d.X[e.laws.C]
        err['temperature_K'] = d.X.Temperature
        err['Mn'] = d.X.MW_Mn
        err['y'] = d.y
        err['p_loco'] = expit(pred['loco'])
        err['loss_loco'] = np.logaddexp(0, pred['loco'])-d.y*pred['loco']
        errors.append(err)
    pd.DataFrame(table).to_csv(e.HERE / 'development_comparison.csv', index=False)
    pd.concat(errors).to_csv(e.HERE / 'development_predictions.csv', index=False)
    e.save(e.HERE / 'development_paired_comparisons.json', comparisons)
    e.save(e.HERE / 'development_concentration_diagnostics.json', slopes)
    identical = d.X.copy()
    identical['chemistry'] = d.chem
    identical['y'] = d.y
    dup = identical.groupby([c for c in identical if c != 'y'], dropna=False).y.agg(['size', 'sum'])
    mixed = (dup['sum']>0)&(dup['sum']<dup['size'])
    e.save(e.HERE / 'development_data_audit.json', dict(n=len(d.y), positive=int(d.y.sum()),
           unique_descriptor_conditions=len(dup), mixed_label_conditions=int(mixed.sum()),
           observations_in_mixed_label_conditions=int(dup.loc[mixed, 'size'].sum()),
           minority_labels_with_identical_inputs=int(np.minimum(dup['sum'], dup['size']-dup['sum']).sum()),
           chemistry_count=len(np.unique(d.chem)), polymer_solvent_groups=len(np.unique(d.groups)),
           ranges={c: [float(d.X[c].min()), float(d.X[c].max())] for c in d.X}))
    print(pd.DataFrame(table).to_string(index=False))
    print(json.dumps(comparisons, indent=2))
    print(json.dumps(slopes, indent=2))


if __name__ == '__main__':
    main()
