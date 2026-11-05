"""Fixed-model diagnostic. Writes only here; never refits or searches equations.

Run from repo root: PYTHONDONTWRITEBYTECODE=1 .venv/bin/python experiments/feature_gap_diagnostic/run_diagnostic.py
Use --summarize to recompute tables from the saved RF predictions without fitting.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from scipy.special import expit
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score, log_loss

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ASTRA = HERE.parent / 'astra_equation_repair'
sys.path.insert(0, str(ROOT / 'symbolic_regression'))
from chiral_sr import data, laws, protocol

MODELS = ['B32', 'RF15', 'RF41']
HOLDOUT_SEED = 20261104  # one seed, fixed before RF evaluation; no seed selection


def read(p):
    return json.loads(p.read_text())


def save(p, obj):
    p.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def source_hashes():
    paths = [p for folder in ['symbolic_regression', 'data_analysis', 'experimental_data',
                             'feature_generation', 'experiments/astra_equation_repair']
             for p in (ROOT / folder).rglob('*')
             if p.is_file() and '__pycache__' not in p.parts and '.mplcache' not in p.parts]
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)}


def inputs():
    d = data.load(ROOT)
    manifest = pd.read_csv(ASTRA / 'split_manifest.csv')
    assert np.array_equal(manifest.row, np.arange(len(d.y)))
    assert np.array_equal(manifest.chemistry, d.chem)
    assert np.array_equal(manifest.experiment_id, d.meta.experiment_id)
    assert sha(ROOT / 'data_analysis/data/analysis_data.xlsx') == read(ASTRA / 'split.json')['workbook_sha256']
    dev = np.flatnonzero(manifest.partition.eq('development'))
    hold = np.flatnonzero(manifest.partition.eq('holdout'))
    fs = read(ASTRA / 'folds.json')
    bdev = np.load(ASTRA / 'results/B32/predictions.npz')
    bhold = np.load(ASTRA / 'holdout_predictions.npz')
    for rows, b in [(dev, bdev), (hold, bhold)]:
        for key, v in [('rows', rows), ('y', d.y[rows]), ('chem', d.chem[rows])]:
            assert np.array_equal(b[key], v)
    for name, folds in [('loco', fs['loco']), *[(f'pair{r}', x) for r, x in enumerate(fs['pair'])]]:
        seen = []
        for f in folds:
            tr, te = dev[f['train']], dev[f['test']]
            assert not set(tr) & set(te)
            assert sorted(f['train'] + f['test']) == list(range(len(dev)))
            group = d.chem if name == 'loco' else d.groups
            assert not set(group[tr]) & set(group[te])
            seen.extend(f['test'])
        assert sorted(seen) == list(range(len(dev)))
    assert not set(d.chem[dev]) & set(d.chem[hold])
    # Confirm B32 saved standardization used precisely the matched training rows.
    terms = laws.LAWS['B32']['terms']
    F = laws.term_values(terms, d.X)
    checks = []
    for f, obj in zip(fs['loco'], read(ASTRA / 'results/B32/loco_fits.json')):
        checks.append((dev[f['train']], dev[f['test']], obj, bdev['loco'][f['test']]))
    checks.append((dev, hold, read(ASTRA / 'results/B32/development_fit.json'), bhold['B32']))
    max_error = 0.
    for tr, te, obj, eta in checks:
        assert obj['terms'] == terms
        np.testing.assert_allclose(obj['mean'], F[tr].mean(0), rtol=1e-12, atol=1e-12)
        sd = F[tr].std(0, ddof=1)
        np.testing.assert_allclose(obj['sd'], np.where(sd > 0, sd, 1), rtol=1e-12, atol=1e-12)
        fit = laws.Fit(terms, obj['intercept'], np.array(obj['effects']), np.array(obj['mean']), np.array(obj['sd']))
        error = np.max(np.abs(fit.logit(F[te]) - eta))
        max_error = max(max_error, float(error))
        assert error < 1e-10
    return d, manifest, dev, hold, fs, bdev, bhold, max_error


def path_contributions(model, X):
    """Exact tree-path probability decomposition, NOT SHAP or causal effects.

    Attribute child-minus-parent positive class probability to the split feature.
    Node probabilities incorporate training class weights. Follow sklearn float32
    comparisons and learned missing directions. Sum + mean root = predict_proba.
    """
    X = np.asarray(X, dtype=np.float32)
    contributions = np.zeros(X.shape)
    baseline = 0.
    for estimator in model.estimators_:
        t = estimator.tree_
        values = t.value[:, 0, 1] / t.value[:, 0, :].sum(axis=1)
        baseline += values[0]
        nodes = np.zeros(len(X), dtype=int)
        while True:
            rows = np.flatnonzero(t.children_left[nodes] != -1)
            if not len(rows):
                break
            parent = nodes[rows]
            features = t.feature[parent]
            x = X[rows, features]
            left = np.where(np.isnan(x), t.missing_go_to_left[parent].astype(bool), x <= t.threshold[parent])
            child = np.where(left, t.children_left[parent], t.children_right[parent])
            contributions[rows, features] += values[child] - values[parent]
            nodes[rows] = child
    return baseline / len(model.estimators_), contributions / len(model.estimators_)


def train(d, dev, hold, fs):
    for name, frame in [('RF15', d.X), ('RF41', d.X41)]:
        path = HERE / f'{name}_predictions.npz'
        if path.exists():
            raise RuntimeError(f'{path.name} exists; use --summarize. No silent refit.')
        X = frame.to_numpy(float)
        loco = np.full(len(dev), np.nan)
        contrib = np.full(X.shape, np.nan)
        roots = np.full(len(d.y), np.nan)
        pair = np.full((20, len(dev)), np.nan)
        support = []

        def fit_predict(tr, te, seed, explain=False):
            model = RandomForestClassifier(**protocol.RF_PARAMS, random_state=seed)
            model.fit(X[tr], d.y[tr])
            p = model.predict_proba(X[te])[:, 1]
            if explain:
                root, attr = path_contributions(model, X[te])
                np.testing.assert_allclose(root + attr.sum(1), p, rtol=0, atol=1e-12)
                contrib[te], roots[te] = attr, root
                for c in np.unique(d.chem[te]):
                    rows = te[d.chem[te] == c]
                    for j, feature in enumerate(frame.columns):
                        a, b = X[tr, j], X[rows, j]
                        support.append(dict(model=name, chemistry=int(c), feature=feature,
                                            train_min=float(np.nanmin(a)), train_max=float(np.nanmax(a)),
                                            test_min=float(np.nanmin(b)) if np.isfinite(b).any() else None,
                                            test_max=float(np.nanmax(b)) if np.isfinite(b).any() else None,
                                            train_missing=int(np.isnan(a).sum()), test_missing=int(np.isnan(b).sum()),
                                            outside_train_range=int(((b < np.nanmin(a)) | (b > np.nanmax(a))).sum()),
                                            train_unique=int(len(np.unique(a[np.isfinite(a)])))))
            return p

        for i, f in enumerate(fs['loco']):
            loco[f['test']] = fit_predict(dev[f['train']], dev[f['test']], 1000 + i, True)
            print(f'{name} LOCO {i + 1}/47', flush=True)
        hp = fit_predict(dev, hold, HOLDOUT_SEED, True)
        # Save primary results before the complementary, weaker transfer protocol.
        np.savez(HERE / f'{name}_primary.npz', loco=loco, holdout=hp, contributions=contrib, roots=roots)
        pd.DataFrame(support).to_csv(HERE / f'{name}_descriptor_support.csv', index=False)
        for r, folds in enumerate(fs['pair']):
            for i, f in enumerate(folds, 1):
                pair[r, f['test']] = fit_predict(dev[f['train']], dev[f['test']], (1000 + r) * 100 + i)
            print(f'{name} pair repeat {r + 1}/20', flush=True)
        assert np.isfinite(loco).all() and np.isfinite(pair).all() and np.isfinite(hp).all()
        np.savez(path, rows=dev, holdout_rows=hold, y=d.y[dev], holdout_y=d.y[hold],
                 loco=loco, pair=pair, holdout=hp, contributions=contrib, roots=roots)
        (HERE / f'{name}_primary.npz').unlink()


def probability_loss(y, p):
    q = np.clip(p, np.finfo(float).eps, 1 - np.finfo(float).eps)
    return -(y * np.log(q) + (1 - y) * np.log1p(-q))


def metrics(y, p, loss):
    both = len(np.unique(y)) == 2
    return dict(ap=float(average_precision_score(y, p)) if both else None,
                roc_auc=float(roc_auc_score(y, p)) if both else None, log_loss=float(np.mean(loss)))


def summarize(d, manifest, dev, hold, fs, bdev, bhold):
    saved = {m: np.load(HERE / f'{m}_predictions.npz') for m in MODELS[1:]}
    for z in saved.values():
        assert np.array_equal(z['rows'], dev) and np.array_equal(z['holdout_rows'], hold)
        assert np.array_equal(z['y'], d.y[dev]) and np.array_equal(z['holdout_y'], d.y[hold])
    eta = np.empty(len(d.y)); eta[dev] = bdev['loco']; eta[hold] = bhold['B32']
    P = {'B32': expit(eta)}
    L = {'B32': np.logaddexp(0, eta) - d.y * eta}
    for m, z in saved.items():
        P[m] = np.empty(len(d.y)); P[m][dev] = z['loco']; P[m][hold] = z['holdout']
        L[m] = probability_loss(d.y, P[m])
    errors = {m: (P[m] >= .5) != d.y for m in MODELS}
    out = d.meta.copy()
    out.insert(0, 'row', np.arange(len(d.y)))
    out['chemistry'], out['partition'], out['y'] = d.chem, manifest.partition, d.y
    out['MW_Mn_g_mol'] = d.X.MW_Mn
    for m in MODELS:
        out[f'{m}_p'], out[f'{m}_loss'], out[f'{m}_error'] = P[m], L[m], errors[m]
    out['all_three_error'] = np.logical_and.reduce(list(errors.values()))
    out.to_csv(HERE / 'predictions.csv', index=False)
    per = []
    for c in np.unique(d.chem):
        ix = np.flatnonzero(d.chem == c)
        item = dict(chemistry=int(c), partition=manifest.partition.iloc[ix[0]],
                    polymers='; '.join(sorted(d.meta.polymer.iloc[ix].unique())), n=len(ix), positives=int(d.y[ix].sum()))
        for m in MODELS:
            item.update({f'{m}_{k}': v for k, v in metrics(d.y[ix], P[m][ix], L[m][ix]).items()})
            item.update({f'{m}_errors': int(errors[m][ix].sum()),
                         f'{m}_false_negatives': int((errors[m][ix] & (d.y[ix] == 1)).sum()),
                         f'{m}_false_positives': int((errors[m][ix] & (d.y[ix] == 0)).sum())})
        for m in MODELS[1:]:
            item[f'{m}_loss_delta_vs_B32'] = float(np.mean(L[m][ix] - L['B32'][ix]))
            item[f'{m}_corrected_B32_errors'] = int((errors['B32'][ix] & ~errors[m][ix]).sum())
            item[f'{m}_new_errors_vs_B32'] = int((~errors['B32'][ix] & errors[m][ix]).sum())
        for a, b in itertools.combinations(MODELS, 2):
            item[f'{a}_{b}_shared_errors'] = int((errors[a][ix] & errors[b][ix]).sum())
        item['all_three_errors'] = int(out.all_three_error.iloc[ix].sum())
        per.append(item)
    per = pd.DataFrame(per)
    per.to_csv(HERE / 'per_chemistry_comparison.csv', index=False)
    scores, overlap, patterns, deltas = [], [], [], []
    for partition, rows in [('development_LOCO', dev), ('holdout', hold)]:
        for m in MODELS:
            scores.append(dict(protocol=partition, model=m, n=len(rows), positives=int(d.y[rows].sum()),
                               **metrics(d.y[rows], P[m][rows], L[m][rows]),
                               macro_log_loss=float(per.loc[per.chemistry.isin(d.chem[rows]), f'{m}_log_loss'].mean()),
                               errors=int(errors[m][rows].sum())))
        for a, b in itertools.combinations(MODELS, 2):
            ea, eb = errors[a][rows], errors[b][rows]
            overlap.append(dict(protocol=partition, model_a=a, model_b=b, both_wrong=int((ea & eb).sum()),
                                a_only_wrong=int((ea & ~eb).sum()), b_only_wrong=int((~ea & eb).sum()),
                                both_correct=int((~ea & ~eb).sum()), jaccard=float((ea & eb).sum() / (ea | eb).sum())))
        for pattern in itertools.product([False, True], repeat=3):
            selected = np.logical_and.reduce([errors[m][rows] == v for m, v in zip(MODELS, pattern)])
            patterns.append(dict(protocol=partition, **dict(zip(MODELS, pattern)), n=int(selected.sum())))
        # Paired chemistry bootstrap of fixed predictions. No refit or selection.
        W = protocol.cluster_weights(d.chem[rows], n_boot=2000, seed=20261105)
        for a, b in [('RF15', 'B32'), ('RF41', 'B32'), ('RF41', 'RF15')]:
            dl = L[a][rows] - L[b][rows]
            bs = (W * dl).sum(1) / W.sum(1)
            group_delta = pd.Series(dl).groupby(d.chem[rows]).mean()
            deltas.append(dict(protocol=partition, comparison=f'{a} minus {b}',
                               delta_log_loss=float(dl.mean()), ci95_low=float(np.quantile(bs, .025)),
                               ci95_high=float(np.quantile(bs, .975)), chemistries_lower_loss=int((group_delta < 0).sum()),
                               chemistries_total=len(group_delta), median_chemistry_delta=float(group_delta.median())))
    repeats = []
    for m in MODELS:
        for r in range(20):
            p = expit(bdev['pair'][r]) if m == 'B32' else saved[m]['pair'][r]
            loss = np.logaddexp(0, bdev['pair'][r]) - d.y[dev] * bdev['pair'][r] if m == 'B32' else probability_loss(d.y[dev], p)
            repeats.append(dict(model=m, repeat=r, **metrics(d.y[dev], p, loss)))
    repeats = pd.DataFrame(repeats)
    repeats.to_csv(HERE / 'pair_repeat_metrics.csv', index=False)
    for m in MODELS:
        sub = repeats[repeats.model == m]
        scores.append(dict(protocol='pair_CV_mean_20_repeats', model=m, n=len(dev), positives=int(d.y[dev].sum()),
                           **{k: float(sub[k].mean()) for k in ['ap', 'roc_auc', 'log_loss']},
                           **{f'{k}_sd': float(sub[k].std(ddof=1)) for k in ['ap', 'roc_auc', 'log_loss']}))
    for name, obj in [('metrics', scores), ('failure_overlap', overlap), ('failure_patterns', patterns), ('paired_loss_intervals', deltas)]:
        pd.DataFrame(obj).to_csv(HERE / f'{name}.csv', index=False)
    save(HERE / 'sensitivity_summary.json', dict(
        holdout_excluding_IDTBT={m: metrics(d.y[hold[d.chem[hold] != 57]],
            P[m][hold[d.chem[hold] != 57]], L[m][hold[d.chem[hold] != 57]]) for m in MODELS},
        pair_CV_mean_fraction_test_chemistry_in_training=float(np.mean([
            sum(np.isin(d.chem[dev[f['test']]], d.chem[dev[f['train']]]).sum() for f in folds) / len(dev)
            for folds in fs['pair']]))))
    # Attribute frozen predictions, including target-specific correction subsets.
    attrs = []
    for m, frame in [('RF15', d.X), ('RF41', d.X41)]:
        z = saved[m]
        np.testing.assert_allclose(z['roots'] + z['contributions'].sum(1), P[m], atol=1e-12, rtol=0)
        for c in np.unique(d.chem):
            subsets = [('all', d.chem == c),
                       ('corrected_B32_errors', (d.chem == c) & errors['B32'] & ~errors[m]),
                       ('persistent_B32_errors', (d.chem == c) & errors['B32'] & errors[m])]
            if m == 'RF41':
                subsets.append(('corrected_RF15_errors', (d.chem == c) & errors['RF15'] & ~errors[m]))
            for subset, mask in subsets:
                ix = np.flatnonzero(mask)
                if not len(ix):
                    continue
                for j, f in enumerate(frame):
                    attrs.append(dict(model=m, chemistry=int(c), subset=subset, n=len(ix), feature=f,
                                      mean_contribution=float(z['contributions'][ix, j].mean()),
                                      mean_abs_contribution=float(np.abs(z['contributions'][ix, j]).mean()),
                                      mean_root_probability=float(z['roots'][ix].mean())))
    pd.DataFrame(attrs).to_csv(HERE / 'descriptor_contributions.csv', index=False)
    audit(d, out)
    # Independently check exported row/chemistry aggregates and sklearn RF loss.
    disk_rows = pd.read_csv(HERE / 'predictions.csv')
    disk_chem = pd.read_csv(HERE / 'per_chemistry_comparison.csv')
    assert len(disk_rows) == 773 and len(disk_chem) == 59 and disk_rows.experiment_id.is_unique
    for partition in ['development', 'holdout']:
        g = disk_rows[disk_rows.partition.eq(partition)]
        c = disk_chem[disk_chem.partition.eq(partition)]
        for m in MODELS:
            np.testing.assert_allclose(np.average(c[f'{m}_log_loss'], weights=c.n),
                                       g[f'{m}_loss'].mean(), rtol=0, atol=1e-12)
            assert c[f'{m}_errors'].sum() == g[f'{m}_error'].sum()
            if m != 'B32':
                np.testing.assert_allclose(log_loss(g.y, g[f'{m}_p']), g[f'{m}_loss'].mean(), rtol=0, atol=1e-12)
    return pd.DataFrame(scores)


def audit(d, predictions):
    raw = pd.read_excel(ROOT / 'data_analysis/data/analysis_data.xlsx', sheet_name='Exp_ML_features')
    summary, conflicts = [], []
    condition = pd.DataFrame(dict(polymer=d.meta.polymer, solvent=d.meta.solvent,
                                  temperature_C=d.meta.temperature_C, concentration_mg_mL=d.meta.concentration_mg_mL))
    keys = {'recorded_condition': condition, 'X15': d.X, 'X41': d.X41,
            'all77': raw.drop(columns='experiment_id')}
    for partition in ['development', 'holdout', 'all']:
        rows = np.arange(len(d.y)) if partition == 'all' else np.flatnonzero(predictions.partition.eq(partition))
        for label, frame in keys.items():
            frame = frame.iloc[rows].copy()
            # Exact equality; missing values constitute a shared recorded value.
            ids = frame.groupby(list(frame), dropna=False, sort=False).ngroup().to_numpy()
            floor = mixed_n = mixed_rows = residual_phase_floor = 0
            for group in np.unique(ids):
                ix = rows[ids == group]; y = d.y[ix]
                minority = min(y.sum(), len(y) - y.sum())
                floor += int(minority)
                if not minority:
                    continue
                mixed_n += 1; mixed_rows += len(ix)
                phase_floor = 0
                for phase in d.meta.assembly_type.iloc[ix].unique():
                    yp = y[d.meta.assembly_type.iloc[ix].eq(phase)]
                    phase_floor += int(min(yp.sum(), len(yp) - yp.sum()))
                residual_phase_floor += phase_floor
                conflicts.append(dict(partition=partition, key=label, condition_id=int(group), n=len(ix), positives=int(y.sum()),
                                      minority_label_floor=int(minority), floor_with_phase=phase_floor,
                                      polymers='; '.join(sorted(d.meta.polymer.iloc[ix].unique())),
                                      solvents='; '.join(sorted(d.meta.solvent.iloc[ix].unique())),
                                      temperatures_C='; '.join(map(str, sorted(d.meta.temperature_C.iloc[ix].unique()))),
                                      concentrations_mg_mL='; '.join(map(str, sorted(d.meta.concentration_mg_mL.iloc[ix].unique()))),
                                      phases='; '.join(sorted(d.meta.assembly_type.iloc[ix].unique())),
                                      g_min=float(d.meta.abs_g_factor.iloc[ix].min()), g_max=float(d.meta.abs_g_factor.iloc[ix].max()),
                                      experiment_ids='; '.join(d.meta.experiment_id.iloc[ix]),
                                      distinct_source_records=int(d.meta.source_record.iloc[ix].nunique()),
                                      **{f'{m}_errors': int(predictions[f'{m}_error'].iloc[ix].sum()) for m in MODELS}))
            summary.append(dict(partition=partition, key=label, n=len(rows), distinct_conditions=len(np.unique(ids)),
                                conflicting_conditions=mixed_n, rows_in_conflicts=mixed_rows,
                                minority_label_floor=floor, minority_label_floor_with_phase=residual_phase_floor))
    pd.DataFrame(summary).to_csv(HERE / 'condition_audit.csv', index=False)
    pd.DataFrame(conflicts).to_csv(HERE / 'conflicting_conditions.csv', index=False)
    phases = []
    for partition in ['development', 'holdout', 'all']:
        sub = predictions if partition == 'all' else predictions[predictions.partition.eq(partition)]
        for phase, g in sub.groupby('assembly_type'):
            phases.append(dict(partition=partition, phase=phase, n=len(g), positives=int(g.y.sum()),
                               positive_rate=float(g.y.mean()), g_min=float(g.abs_g_factor.min()), g_max=float(g.abs_g_factor.max()),
                               **{f'{m}_errors': int(g[f'{m}_error'].sum()) for m in MODELS},
                               **{f'{m}_log_loss': float(g[f'{m}_loss'].mean()) for m in MODELS},
                               all_three_errors=int(g.all_three_error.sum())))
    pd.DataFrame(phases).to_csv(HERE / 'phase_audit.csv', index=False)
    save(HERE / 'label_audit.json', dict(n=len(d.y), positives=int(d.y.sum()),
        metadata_columns=d.meta.columns.tolist(), missing_labels=int(d.meta.abs_g_factor.isna().sum()),
        negative_absolute_g=int((d.meta.abs_g_factor < 0).sum()), zero_g=int(d.meta.abs_g_factor.eq(0).sum()),
        exactly_at_threshold=int(d.meta.abs_g_factor.eq(.001).sum()),
        between_half_and_twice_threshold=int(d.meta.abs_g_factor.between(.0005, .002).sum()),
        distinct_source_records=int(d.meta.source_record.nunique()),
        distinct_dataset_source_pairs=int(d.meta[['dataset_split', 'source_record']].drop_duplicates().shape[0]),
        missing_X15=d.X.isna().sum().to_dict(), missing_X41=d.X41.isna().sum().to_dict(),
        shared_features=sorted(set(d.X) & {c.removesuffix(data.N) for c in d.X41}),
        only_X15=sorted(set(d.X) - {c.removesuffix(data.N) for c in d.X41}),
        only_X41=sorted({c.removesuffix(data.N) for c in d.X41} - set(d.X))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--summarize', action='store_true')
    args = parser.parse_args()
    before = source_hashes()
    d, manifest, dev, hold, fs, bdev, bhold, berror = inputs()
    if not args.summarize:
        save(HERE / 'provenance.json', dict(source_hashes=before, rf_params=protocol.RF_PARAMS,
             seeds=dict(loco='1000 + index in sorted development chemistries', pair='(1000 + repeat) * 100 + fold (1..5)', holdout=HOLDOUT_SEED),
             versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, sklearn=sklearn.__version__),
             feature_sets=dict(RF15=list(d.X), RF41=list(d.X41)), preprocessing='RF: raw inputs, native NaN routing, no scaling/imputation; B32: saved train-fold mean/sample SD',
             scope='fixed B32 predictions; two fixed original RF configurations, no model or feature selection, no calibration or threshold tuning'))
        train(d, dev, hold, fs)
    else:
        assert before == read(HERE / 'provenance.json')['source_hashes']
    scores = summarize(d, manifest, dev, hold, fs, bdev, bhold)
    assert before == source_hashes(), 'A source or completed experiment file changed.'
    save(HERE / 'verification.json', dict(source_files_unchanged=len(before), matched_rows_and_folds=True,
         no_chemistry_overlap_in_LOCO_or_holdout=True, B32_LOCO_and_holdout_train_only_standardization_verified=True,
         max_B32_logit_reconstruction_error=berror, RF_path_probability_reconstruction_tolerance=1e-12,
         observations=773, development_chemistries=47, holdout_chemistries=12, fixed_RF_fits=296,
         exported_chemistry_row_reconciliation_passed=True, independent_sklearn_RF_loss_checks_passed=True,
         diagnostic_script_sha256=sha(Path(__file__)),
         output_sha256={p.name: sha(p) for p in sorted(HERE.iterdir()) if p.suffix in ['.csv', '.npz']}))
    print(scores.to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
