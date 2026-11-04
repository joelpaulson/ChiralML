"""Reproduce fixed numerical results without redrawing or selecting models.

Recomputes all development folds and final all-data coefficients. Holdout metrics
are verified from its saved single-evaluation predictions, with no new candidate
evaluation. Generates verification.json. Run after the completed experiment.
"""
import json
import re
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
import experiment as e
import final_equation as printed


def main():
    report = {}
    cfg = e.read(e.HERE / 'split.json')
    selection = e.read(e.HERE / 'selection.json')
    assert e.sha(e.HERE / 'registry.json') == selection['registry_sha256']
    assert e.sha(e.HERE / 'folds.json') == selection['folds_sha256']
    for path, digest in cfg['source_hashes'].items():
        assert e.sha(e.ROOT / path) == digest
    assert e.sha(e.SOURCE) == cfg['workbook_sha256']
    features = pd.read_excel(e.SOURCE, sheet_name='Exp_ML_features')
    meta = pd.read_excel(e.SOURCE, sheet_name='Exp_ML_metadata', usecols=['experiment_id', 'polymer', 'solvent'])
    chem = e.data.chemistry(features, meta)
    hold = np.sort(np.random.default_rng(cfg['seed']).permutation(np.unique(chem))[:12])
    manifest = pd.read_csv(e.HERE / 'split_manifest.csv')
    assert np.array_equal(hold, cfg['holdout_groups'])
    assert np.array_equal(chem, manifest.chemistry)
    assert np.array_equal(np.isin(chem, hold), manifest.partition.eq('holdout'))
    assert manifest.experiment_id.is_unique and features.experiment_id.equals(meta.experiment_id)
    assert (manifest.groupby('polymer').partition.nunique() == 1).all()
    report['split_and_source_hashes'] = 'pass'
    d = e.load_partition()
    folds = e.folds(d)
    registry = e.read(e.HERE / 'registry.json')
    assert len(registry) == 4
    report['cv_max_logit_difference'] = {}
    for name, obj in registry.items():
        F = e.laws.term_values(obj['terms'], d.X)
        saved = np.load(e.HERE / 'results' / name / 'predictions.npz')
        diffs = []
        for fs, previous, groups in [(folds['loco'], saved['loco'], d.chem),
                      *[(fs, old, d.groups) for fs, old in zip(folds['pair'], saved['pair'])]]:
            eta = np.full(len(d.y), np.nan)
            for fold in fs:
                tr, te = np.array(fold['train']), np.array(fold['test'])
                assert not set(groups[tr]) & set(groups[te])
                fit = e.fitted(F[tr], d.y[tr], obj['terms'])
                eta[te] = fit.logit(F[te])
            assert np.isfinite(eta).all()
            diffs.append(float(np.max(np.abs(eta-previous))))
        report['cv_max_logit_difference'][name] = max(diffs)
        assert max(diffs) < 1e-10
        old = e.read(e.HERE / 'results' / name / 'metrics.json')
        assert all(abs(e.metrics(d.y, saved['loco'])[k]-old['loco'][k]) < 1e-12 for k in old['loco'])
        print(f'{name}: all 147 development fits reproduced; maximum logit difference {max(diffs):.3g}', flush=True)
    d = e.load_partition('all')
    obj = e.read(e.HERE / 'final_B32_fit.json')
    F = e.laws.term_values(obj['terms'], d.X)
    fit = e.fitted(F, d.y, obj['terms'])
    original = e.unpack(obj)
    report['all_data_refit_max_logit_difference'] = float(np.max(np.abs(fit.logit(F)-original.logit(F))))
    assert report['all_data_refit_max_logit_difference'] < 1e-10
    report['printed_max_logit_error'] = float(np.max(np.abs(printed.logit(d.X)-original.logit(F))))
    report['printed_max_probability_error'] = float(np.max(np.abs(printed.predict(d.X)-expit(original.logit(F)))))
    assert report['printed_max_logit_error'] < 1e-10
    assert report['printed_max_probability_error'] < 1e-12
    doc = (e.HERE / 'FINAL_EQUATION.md').read_text()
    text = re.search(r'```python\nf = \((.*?)\n\)\n```', doc, flags=re.S).group(1)
    lookup = dict(c='Concentration_mg_mL', M='MW_Mn', S='SC_C_atoms', B='BB_C_atoms', R='SC_n_branch_points',
                  U='SC_n_units', H='Solvent_NumHAcceptors', F='Backbone_fastsolv', Q='BB_conjugated_units', E='BB_hetero_atoms')
    # This is a local, authored expression with only arithmetic; evaluates the
    # actual printed transcription rather than recalling serialized coefficients.
    doc_eta = eval('('+text+')', {'__builtins__': {}}, {k: d.X[v].to_numpy() for k, v in lookup.items()})
    report['markdown_equation_max_probability_error'] = float(np.max(np.abs(expit(doc_eta)-expit(original.logit(F)))))
    assert report['markdown_equation_max_probability_error'] < 1e-12
    # Independent numerical solver, same objective, initialized at zero.
    Z = np.column_stack([np.ones(len(d.y)), (F-fit.mean)/fit.sd])
    penalty = np.r_[0., np.full(F.shape[1], e.laws.RIDGE)]
    def objective(theta):
        eta = (Z*theta).sum(axis=1)
        value = np.sum(np.logaddexp(0, eta)-d.y*eta)+.5*np.dot(penalty, theta**2)
        grad = (Z*(expit(eta)-d.y)[:, None]).sum(axis=0)+penalty*theta
        return value, grad
    alt = minimize(objective, np.zeros(Z.shape[1]), jac=True, method='BFGS', options={'gtol': 1e-7, 'maxiter': 1000})
    report['independent_optimizer_max_probability_difference'] = float(np.max(np.abs(expit((Z*alt.x).sum(1))-expit(original.logit(F)))))
    report['independent_optimizer_max_gradient'] = float(np.max(np.abs(alt.jac)))
    report['independent_optimizer_success'] = bool(alt.success)
    assert alt.success and report['independent_optimizer_max_gradient'] < 1e-7
    assert report['independent_optimizer_max_probability_difference'] < 1e-7
    # Audit the single frozen holdout result from saved predictions, not a new
    # trial and not a test of any rejected candidate.
    saved = np.load(e.HERE / 'holdout_predictions.npz')
    expected = e.read(e.HERE / 'holdout_metrics.json')['models']['B32']
    for k, v in e.metrics(saved['y'], saved['B32']).items():
        assert abs(v-expected[k]) < 1e-12
    report['holdout_saved_metric_recalculation'] = 'pass; no new candidate predictions'
    report['verified_utc'] = e.stamp()
    e.save(e.HERE / 'verification.json', report)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
