"""Staged, auditable equation repair. Run from any directory with the repo's .venv.

setup reads identifiers/features only. evaluate reads outcome cells only for the
development rows. final requires a frozen selection and is a one-shot operation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'symbolic_regression'))
os.environ.setdefault('MPLCONFIGDIR', str(HERE / '.mplcache'))
import numpy as np
import pandas as pd
import scipy
from scipy.special import expit
import sklearn
from sklearn.metrics import average_precision_score, roc_auc_score
import torch
from chiral_sr import data, laws, protocol

torch.set_num_threads(1)
SOURCE = ROOT / 'data_analysis/data/analysis_data.xlsx'
SEED = 20261104


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def read(path):
    return json.loads(Path(path).read_text())


def log(title, body):
    with (HERE / 'EXPERIMENT_LOG.md').open('a') as f:
        f.write(f'\n## {stamp()} — {title}\n\n{body}\n')


def setup():
    if (HERE / 'split.json').exists():
        raise RuntimeError('Split already exists; never redraw it.')
    features = pd.read_excel(SOURCE, sheet_name='Exp_ML_features')
    meta = pd.read_excel(SOURCE, sheet_name='Exp_ML_metadata',
                         usecols=['experiment_id', 'polymer', 'solvent'])
    assert features.experiment_id.equals(meta.experiment_id)
    chem = data.chemistry(features, meta)
    levels = np.unique(chem)
    hold = np.sort(np.random.default_rng(SEED).permutation(levels)[:int(np.ceil(.2 * len(levels)))])
    split = np.where(np.isin(chem, hold), 'holdout', 'development')
    manifest = meta.copy()
    manifest.insert(0, 'row', np.arange(len(meta)))
    manifest['chemistry'] = chem
    manifest['partition'] = split
    manifest.to_csv(HERE / 'split_manifest.csv', index=False)
    cols = [c for c in features.drop(columns='experiment_id')
            if c not in ('Temperature', 'Concentration_mg_mL', 'MW_Mn', 'MW_Mw', 'PDI')
            and (features[c].groupby(meta['polymer']).nunique(dropna=False) <= 1).all()]
    assert (manifest.groupby('polymer').chemistry.nunique() == 1).all()
    config = dict(seed=SEED, algorithm='numpy.default_rng(seed).permutation(sorted chemistry codes); first ceil(0.2*G)',
                  created_utc=stamp(), workbook_sha256=sha(SOURCE), group_columns=cols,
                  group_rule='chiral_sr.data.chemistry: polymer-constant descriptors except Mn, Mw, PDI and conditions; float, round(8), missing=-1e300, concatenate, first-occurrence factorize',
                  development_groups=levels[~np.isin(levels, hold)].tolist(), holdout_groups=hold.tolist(),
                  development_n=int((split == 'development').sum()), holdout_n=int((split == 'holdout').sum()),
                  source_hashes={p: sha(ROOT / p) for p in ['symbolic_regression/chiral_sr/data.py', 'symbolic_regression/chiral_sr/laws.py', 'symbolic_regression/chiral_sr/protocol.py', 'symbolic_regression/chiral_sr/symantic/fast/classification.py']},
                  versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                                scipy=scipy.__version__, sklearn=sklearn.__version__, torch=torch.__version__))
    save(HERE / 'split.json', config)
    save(HERE / 'registry.json', {'B32': dict(terms=laws.LAWS['B32']['terms'], registered_utc=stamp(), rationale='Original structure, coefficients refitted on development data.')})
    (HERE / 'EXPERIMENT_LOG.md').write_text('# Chronological equation-repair experiment\n\nDate: 2026-11-04 (America/Chicago); entries use UTC. No outcome data were inspected before the split.\n')
    log('Protocol and partition fixed before outcome inspection',
        f'Seed {SEED}; {len(levels)} chemistry groups, {len(hold)} held aside. Development n={config["development_n"]}; holdout n={config["holdout_n"]}. '
        'The manifest contains identifiers only. Group definitions and all defining columns are in split.json. '
        'All molecular-weight variants sharing the polymer descriptors stay together. No stratification or redraw.\n\n'
        'Primary: leave-one-chemistry-out (LOCO), pooled AP, ROC-AUC and stable log loss, plus equal-chemistry loss and mixed-class within-chemistry metrics. '
        'Complement: existing 20 repeats of 5-fold StratifiedGroupKFold, polymer-name || solvent, seeds 1000–1019. '
        'All models use identical saved folds, train-only term standardization (sample SD), summed logistic loss + 0.05 sum(beta_standardized^2), unpenalized intercept, SyMANTIC fit_glm, no chemistry offsets, no hyperparameter tuning. '
        'At most three new structures, each registered before evaluation. Prefer a five-term structure unless a sixth term brings consistent LOCO benefit. '
        'Assess transfer using AP together with log loss and chemistry-wise consistency; do not select on training scores or complementary CV alone. '
        'Bootstrap intervals resample chemistries of fixed OOF predictions, so omit selection and refitting uncertainty. '
        'Freeze selection before final labels are read. Final evaluates only selected and B32 once, then refits selected on all data. '
        'B32 was originally discovered on the full dataset: this is not a pristine prospective holdout.\n\n'
        'Original B32: f=b0+b1*(S-B)*c*R+b2*H*U*F+b3*c^2*F/Q+b4*(S-B)*U*B+b5*c*M*B/E. '
        'S=SC_C_atoms, B=BB_C_atoms, R=SC_n_branch_points, U=SC_n_units, H=Solvent_NumHAcceptors, F=Backbone_fastsolv, Q=BB_conjugated_units, E=BB_hetero_atoms, M=MW_Mn, c=Concentration_mg_mL. All count descriptors are normalized mol/kg. '
        'This is the original structure with numerical coefficients deliberately refitted within each training partition.')
    print(json.dumps({k: config[k] for k in ['seed', 'development_groups', 'holdout_groups', 'development_n', 'holdout_n']}, indent=2))


def load_partition(partition='development'):
    cfg = read(HERE / 'split.json')
    assert sha(SOURCE) == cfg['workbook_sha256'], 'Workbook changed.'
    manifest = pd.read_csv(HERE / 'split_manifest.csv')
    rows = manifest.row.to_numpy() if partition == 'all' else manifest.loc[manifest.partition == partition, 'row'].to_numpy()
    allowed = set(rows.tolist())
    # Skip held-out worksheet rows before parsing outcomes into a DataFrame.
    skip = lambda i: i > 0 and i - 1 not in allowed
    features = pd.read_excel(SOURCE, sheet_name='Exp_ML_features', skiprows=skip)
    meta = pd.read_excel(SOURCE, sheet_name='Exp_ML_metadata', skiprows=skip)
    assert features.experiment_id.equals(meta.experiment_id)
    assert meta.experiment_id.astype(str).tolist() == manifest.iloc[rows].experiment_id.astype(str).tolist()
    assert meta.abs_g_factor.notna().all()
    X = features[list(data.INPUTS)].astype(float)
    X.columns = [c.removesuffix(data.N) for c in X]
    norm = lambda s: s.astype('string').str.replace(r'\s+', ' ', regex=True).str.lower()
    groups = (norm(meta.polymer) + ' || ' + norm(meta.solvent)).to_numpy(str)
    return SimpleNamespace(X=X, meta=meta, y=(meta.abs_g_factor >= data.G_THRESHOLD).to_numpy(int),
                           chem=manifest.iloc[rows].chemistry.to_numpy(int), groups=groups, rows=rows)


def folds(d):
    path = HERE / 'folds.json'
    if path.exists():
        return read(path)
    rows = np.arange(len(d.y))
    out = {'loco': [{'train': rows[d.chem != c].tolist(), 'test': rows[d.chem == c].tolist(), 'chemistry': int(c)} for c in np.unique(d.chem)],
           'pair': [[{'train': tr.tolist(), 'test': te.tolist()} for tr, te in protocol.authors_folds(d, rows, r)] for r in range(20)]}
    for key, fs in [('loco', out['loco']), *[(f'pair{r}', fs) for r, fs in enumerate(out['pair'])]]:
        seen = []
        for f in fs:
            tr, te = np.array(f['train']), np.array(f['test'])
            groups = d.chem if key == 'loco' else d.groups
            assert not set(groups[tr]) & set(groups[te])
            seen.extend(te.tolist())
        assert sorted(seen) == rows.tolist()
    save(path, out)
    return out


def metrics(y, eta, weight=None):
    p = expit(eta)
    both = len(np.unique(y if weight is None else y[weight > 0])) == 2
    return dict(ap=float(average_precision_score(y, p, sample_weight=weight)) if both else None,
                auc=float(roc_auc_score(y, p, sample_weight=weight)) if both else None,
                log_loss=float(np.average(np.logaddexp(0, eta) - y * eta, weights=weight)),
                brier=float(np.average((p-y)**2, weights=weight)))


def fitted(F, y, terms):
    fit = laws.fit(F, y, np.arange(len(y)), terms)
    Z = (F-fit.mean)/fit.sd
    residual = expit(fit.logit(F))-y
    grad = np.r_[residual.sum(), Z.T @ residual + laws.RIDGE*fit.effects]
    assert np.max(np.abs(grad))/len(y) < 1e-7, f'Optimizer gradient: {grad}'
    return fit


def pack(fit):
    raw = fit.effects / fit.sd
    return dict(terms=fit.terms, intercept=fit.intercept, effects=fit.effects.tolist(),
                mean=fit.mean.tolist(), sd=fit.sd.tolist(), raw_intercept=float(fit.intercept-fit.mean@raw), raw_coefficients=raw.tolist())


def unpack(obj):
    return laws.Fit(obj['terms'], obj['intercept'], np.array(obj['effects']), np.array(obj['mean']), np.array(obj['sd']))


def register(name, terms, rationale):
    assert not (HERE / 'selection.json').exists(), 'Structures frozen.'
    registry = read(HERE / 'registry.json')
    assert name not in registry and len(registry) < 4
    assert len(terms) <= 6
    registry[name] = dict(terms=terms, rationale=rationale, registered_utc=stamp())
    save(HERE / 'registry.json', registry)
    log(f'PROPOSAL {name} — before any fit or evaluation', rationale + '\n\nExact terms (with free intercept and one fitted coefficient each):\n\n' + '\n'.join(f'- `{t}`' for t in terms))


def evaluate(name):
    assert not (HERE / 'selection.json').exists(), 'Development is frozen.'
    registry = read(HERE / 'registry.json')
    terms = registry[name]['terms']
    out = HERE / 'results' / name
    assert not (out / 'metrics.json').exists(), 'Already evaluated.'
    out.mkdir(parents=True, exist_ok=True)
    d = load_partition()
    fs = folds(d)
    F = laws.term_values(terms, d.X)
    def predict(fs):
        eta = np.full(len(d.y), np.nan)
        coeffs = []
        for fold in fs:
            tr, te = np.array(fold['train']), np.array(fold['test'])
            fit = fitted(F[tr], d.y[tr], terms)
            eta[te] = fit.logit(F[te])
            coeffs.append(pack(fit))
        assert np.isfinite(eta).all()
        return eta, coeffs
    loco, loco_fits = predict(fs['loco'])
    pair = np.stack([predict(repeat)[0] for repeat in fs['pair']])
    fit = fitted(F, d.y, terms)
    per = []
    for c in np.unique(d.chem):
        ix = d.chem == c
        per.append(dict(chemistry=int(c), n=int(ix.sum()), positives=int(d.y[ix].sum()),
                        polymers='; '.join(sorted(d.meta.loc[ix, 'polymer'].unique())), **metrics(d.y[ix], loco[ix])))
    per = pd.DataFrame(per)
    per.to_csv(out / 'per_chemistry.csv', index=False)
    m_pair = pd.DataFrame([metrics(d.y, p) for p in pair])
    m_pair.to_csv(out / 'pair_repeat_metrics.csv', index=False)
    summary = dict(model=name, n=len(d.y), positives=int(d.y.sum()), terms=len(terms),
                   loco=metrics(d.y, loco), chemistry_macro={k: float(per[k].mean()) for k in ['ap', 'auc', 'log_loss', 'brier']},
                   mixed_class_chemistries=int(per.ap.notna().sum()), pair_mean=m_pair.mean().to_dict(),
                   pair_sd=m_pair.std().to_dict(), training=metrics(d.y, fit.logit(F)), evaluated_utc=stamp())
    save(out / 'metrics.json', summary)
    save(out / 'development_fit.json', pack(fit))
    save(out / 'loco_fits.json', loco_fits)
    np.savez(out / 'predictions.npz', rows=d.rows, y=d.y, chem=d.chem, loco=loco, pair=pair, training=fit.logit(F))
    log(f'Feedback for {name}', 'Development only.\n\n```json\n' + json.dumps(summary, indent=2) + '\n```\n\nPer-chemistry results, all OOF predictions, and fitted coefficients saved under results/' + name + '/.')
    print(json.dumps(summary, indent=2))


def freeze(name, reason):
    assert not (HERE / 'selection.json').exists()
    registry = read(HERE / 'registry.json')
    assert all((HERE / 'results' / m / 'metrics.json').exists() for m in registry)
    selection = dict(selected=name, reason=reason, frozen_utc=stamp(), registry_sha256=sha(HERE / 'registry.json'),
                     terms=registry[name]['terms'], folds_sha256=sha(HERE / 'folds.json'),
                     development_results_sha256={m: sha(HERE / 'results' / m / 'metrics.json') for m in registry})
    save(HERE / 'selection.json', selection)
    log('STRUCTURE FROZEN before holdout labels or predictions', f'Selected {name}. {reason}\n\nNo further modifications or changes to selection are permitted. Registry and fold hashes stored in selection.json.')


def final():
    selection = read(HERE / 'selection.json')
    assert sha(HERE / 'registry.json') == selection['registry_sha256']
    assert not (HERE / 'final_started.json').exists(), 'One-shot holdout already opened.'
    save(HERE / 'final_started.json', dict(opened_utc=stamp(), selected=selection['selected']))
    log('Final holdout opened', 'Frozen model and B32 only; coefficients previously fitted on development. This is the first outcome read for held-aside rows.')
    d = load_partition('holdout')
    output, preds = {}, {}
    for name in dict.fromkeys(['B32', selection['selected']]):
        fit = unpack(read(HERE / 'results' / name / 'development_fit.json'))
        eta = fit.logit(laws.term_values(fit.terms, d.X))
        preds[name] = eta
        per = []
        for c in np.unique(d.chem):
            ix = d.chem == c
            per.append(dict(chemistry=int(c), n=int(ix.sum()), positives=int(d.y[ix].sum()), polymers='; '.join(sorted(d.meta.loc[ix, 'polymer'].unique())), **metrics(d.y[ix], eta[ix])))
        pd.DataFrame(per).to_csv(HERE / f'holdout_{name}_per_chemistry.csv', index=False)
        output[name] = metrics(d.y, eta)
    save(HERE / 'holdout_metrics.json', dict(n=len(d.y), positives=int(d.y.sum()), groups=len(np.unique(d.chem)), models=output))
    np.savez(HERE / 'holdout_predictions.npz', rows=d.rows, y=d.y, chem=d.chem, **preds)
    log('Final holdout results (selection unchanged)', '```json\n' + json.dumps(output, indent=2) + '\n```')
    all_d = load_partition('all')
    assert len(all_d.y) == 773
    # B32 refit supplies the requested algebraic comparator, not a new candidate.
    for name in dict.fromkeys(['B32', selection['selected']]):
        terms = read(HERE / 'registry.json')[name]['terms']
        fit = fitted(laws.term_values(terms, all_d.X), all_d.y, terms)
        save(HERE / f'final_{name}_fit.json', pack(fit))
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['setup', 'evaluate', 'final'])
    parser.add_argument('--model', default='B32')
    args = parser.parse_args()
    {'setup': setup, 'evaluate': lambda: evaluate(args.model), 'final': final}[args.stage]()
