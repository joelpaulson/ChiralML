"""Two concentration figures: development LOCO comparisons and frozen holdout.

Points are observed fractions above |g|>=1e-3 at identical descriptor conditions,
not g-factor magnitudes. Curves never fit the plotted chemistry's outcomes.
"""
import experiment as e
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import expit

COLORS = {'B32': '#123e66', 'A_overlap': '#df8b23', 'B_contrast': '#a869a4', 'C_selectivity': '#179381'}
LABELS = {'B32': 'B32 (selected)', 'A_overlap': 'A: overlap', 'B_contrast': 'B: solvation contrast', 'C_selectivity': 'C: added selectivity'}


def condition(d, polymer, solvent, temperature=298.15):
    ix = (d.meta.polymer.str.strip() == polymer) & (d.meta.solvent.str.strip() == solvent) & np.isclose(d.X.Temperature, temperature)
    assert ix.any(), (polymer, solvent, temperature)
    # All displayed rows must differ only in concentration, not another descriptor.
    assert (d.X.loc[ix].drop(columns=e.laws.C).nunique() == 1).all()
    return np.flatnonzero(ix)


def panel(ax, d, indices, fits, title):
    concentrations = d.X.iloc[indices][e.laws.C].to_numpy()
    grid = np.linspace(concentrations.min(), concentrations.max(), 201)
    X = pd.DataFrame(np.repeat(d.X.iloc[[indices[0]]].to_numpy(), len(grid), axis=0), columns=d.X.columns)
    X[e.laws.C] = grid
    for name, fit in fits.items():
        ax.plot(grid, expit(fit.logit(e.laws.term_values(fit.terms, X))), label=LABELS[name], color=COLORS[name], lw=2,
                linestyle='-' if name in ['B32', 'C_selectivity'] else '--')
    points = pd.DataFrame({'c': concentrations, 'y': d.y[indices]}).groupby('c').y.agg(['mean', 'size'])
    ax.scatter(points.index, points['mean'], s=30+12*points['size'], color='#20252b', edgecolor='white', zorder=6, label='Observed strong fraction')
    for c, row in points.iterrows():
        ax.annotate(f'n={int(row["size"])}', (c, row['mean']), xytext=(0, -16 if row['mean'] > .75 else 9),
                    textcoords='offset points', ha='center', fontsize=8)
    ax.set_title(title, fontsize=10, loc='left', pad=10)
    ax.set_ylim(-.09, 1.12)
    ax.set_yticks([0, .5, 1])
    ax.grid(alpha=.16)
    ax.set_xlabel('Concentration (mg/mL)')
    ax.spines[['top', 'right']].set_visible(False)


def main():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    folder = e.HERE / 'figures'
    folder.mkdir(exist_ok=True)
    d = e.load_partition()
    examples = [('N2200 - 9.1 k', 'Chlorobenzene'), ('N2200 - 9.1 k', 'Anisole'),
                ('PII-2T (C10-C10) - 23.8k', 'Toluene'), ('g42T-TT', 'Anisole'),
                ('AcDOT-[OE3]2', 'Chloroform'), ('AcDOT-[C10]2-co-[OE3]2', 'Chloroform')]
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.7), sharey=True)
    allfits = {name: e.read(e.HERE / 'results' / name / 'loco_fits.json') for name in e.read(e.HERE / 'registry.json')}
    folds = e.folds(d)['loco']
    selection = []
    for ax, (poly, solvent) in zip(axes.flat, examples):
        ix = condition(d, poly, solvent)
        c = int(d.chem[ix[0]])
        j = next(i for i, f in enumerate(folds) if f['chemistry'] == c)
        fits = {name: e.unpack(objs[j]) for name, objs in allfits.items()}
        panel(ax, d, ix, fits, f'{poly}\n{solvent}, 25 °C')
        selection.append(dict(partition='development', polymer=poly, solvent=solvent, temperature_K=298.15,
                              row_indices=d.rows[ix].tolist(), fit='LOCO excluding entire chemistry'))
    for ax in axes[:, 0]:
        ax.set_ylabel('P(strong) / observed strong fraction')
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=5, frameon=False, fontsize=9)
    fig.suptitle('Development: concentration response with the entire chemistry left out', fontsize=15, y=.985)
    fig.text(.5, .93, 'Diagnostic examples include targeted errors and adverse effects; curves span measured concentrations only.', ha='center', fontsize=10, color='#555555')
    fig.tight_layout(rect=[0, .05, 1, .91])
    fig.savefig(folder / 'development_concentration.png', dpi=180)
    fig.savefig(folder / 'development_concentration.pdf')
    plt.close(fig)

    d = e.load_partition('holdout')
    # Chosen to expose varied final behavior, including worst failure; explicitly descriptive.
    examples = [('PBTTT-C14', 'Chlorobenzene'), ('PTQ-6bO', 'o-Dichlorobenzene'), ('N2200-OEG', 'Anisole'),
                ('DPP2T-T2M2', 'Chlorobenzene'), ('IDTBT', 'Chlorobenzene'), ('F8BT', 'Chlorobenzene')]
    fit = e.unpack(e.read(e.HERE / 'results/B32/development_fit.json'))
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.7), sharey=True)
    for ax, (poly, solvent) in zip(axes.flat, examples):
        ix = condition(d, poly, solvent)
        panel(ax, d, ix, {'B32': fit}, f'{poly}\n{solvent}, 25 °C')
        selection.append(dict(partition='holdout', polymer=poly, solvent=solvent, temperature_K=298.15,
                              row_indices=d.rows[ix].tolist(), fit='frozen development fit'))
    for ax in axes[:, 0]:
        ax.set_ylabel('P(strong) / observed strong fraction')
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, frameon=False)
    fig.suptitle('Final chemistry holdout: the frozen B32 equation has uneven transfer', fontsize=15, y=.985)
    fig.text(.5, .93, 'Coefficients use development data only. Examples show successes and failures; no model changes followed.', ha='center', fontsize=10, color='#555555')
    fig.tight_layout(rect=[0, .05, 1, .91])
    fig.savefig(folder / 'holdout_concentration.png', dpi=180)
    fig.savefig(folder / 'holdout_concentration.pdf')
    plt.close(fig)
    e.save(folder / 'plotted_conditions.json', selection)


if __name__ == '__main__':
    main()
