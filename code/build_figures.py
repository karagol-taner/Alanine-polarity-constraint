"""Build the colour figures. Secondary encoding is used throughout so they
remain readable in greyscale.

    python3 code/build_figures.py

Writes Figure1, Figure2, Figure3 and Figure5 to figures/ as PNG and SVG.
Figure 4 is not generated here.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd, numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D

# Okabe-Ito pair, colourblind-safe
BLUE, VERM, NEUTRAL = '#0072B2', '#D55E00', '#EDEAE4'
# diverging: two hues with a neutral grey midpoint
DIV = LinearSegmentedColormap.from_list('bg_r', [BLUE, '#7FB9DC', NEUTRAL, '#E8A06A', VERM])
INK, MUTED, GRID = '#1A1A1A', '#5A5A5A', '#D8D6D1'

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 7.5,
    'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': MUTED,
    'xtick.color': MUTED, 'ytick.color': MUTED,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'savefig.dpi': 400, 'figure.dpi': 140, 'figure.facecolor': 'white'
})

import os, sys
ROOT = os.environ.get('PAPER_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def find(name, *dirs):
    """Locate an input under any of the candidate directories."""
    for d in dirs:
        p = os.path.join(ROOT, d, name)
        if os.path.exists(p):
            return p
    sys.exit('cannot find %s under any of %s (root %s)'
             % (name, ', '.join(dirs), ROOT))


X = find('TM_vs_nonTM_results.xlsx', 'results', 'tables')
OUT = os.path.join(ROOT, 'figures') + os.sep
ADJ = find('effects_all_adj.csv', 'results', 'consurf')
cells = pd.read_excel(X, sheet_name='Cells_TM_vs_nonTM')
pp = pd.read_excel(X, sheet_name='Per_protein_effects')
pairs = pd.read_excel(X, sheet_name='Pair_correlations')
AA = list('AVLIMFWCGPSTYNQHKRDE')   # nonpolar block, then polar block

# ============================== FIGURE 1 ==============================
fig, ax = plt.subplots(figsize=(4.9, 4.3))
M = np.full((20, 20), np.nan)
for _, r in cells.iterrows():
    if r.aa1 in AA and r.aa2 in AA:
        M[AA.index(r.aa1), AA.index(r.aa2)] = r.delta_TM_minus_nonTM
vmax = np.nanmax(np.abs(M))
im = ax.imshow(M, cmap=DIV, norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
# secondary encoding: the sign is also carried by a glyph
for i in range(20):
    for j in range(20):
        if not np.isnan(M[i, j]) and M[i, j] > 0:
            ax.plot(j, i, marker='+', ms=3.0, mew=0.75, color=INK)
ax.set_xticks(range(20)); ax.set_xticklabels(AA, fontsize=6.5)
ax.set_yticks(range(20)); ax.set_yticklabels(AA, fontsize=6.5)
ax.set_xlabel('substituted to', fontsize=8); ax.set_ylabel('substituted from', fontsize=8)
for s in ('top', 'right'):
    ax.spines[s].set_visible(True); ax.spines[s].set_linewidth(0.6); ax.spines[s].set_color(MUTED)
ax.axvline(9.5, color=INK, lw=0.9); ax.axhline(9.5, color=INK, lw=0.9)
ax.text(4.5, -1.55, 'nonpolar', ha='center', fontsize=7, color=MUTED)
ax.text(14.5, -1.55, 'polar', ha='center', fontsize=7, color=MUTED)
ax.text(-2.4, 4.5, 'nonpolar', va='center', rotation=90, fontsize=7, color=MUTED)
ax.text(-2.4, 14.5, 'polar', va='center', rotation=90, fontsize=7, color=MUTED)
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
cb.outline.set_linewidth(0.5); cb.outline.set_edgecolor(MUTED)
cb.set_label('median pathogenicity, membrane $-$ soluble', fontsize=7)
cb.ax.tick_params(labelsize=6.5)
ax.annotate('', xy=(19.6, 0), xytext=(20.9, 0),
            arrowprops=dict(arrowstyle='-', lw=0))
ax.legend(handles=[Line2D([], [], marker='+', ls='', color=INK, ms=4, mew=0.75,
                          label='costlier in the membrane (35 of 380 classes)')],
          loc='upper left', bbox_to_anchor=(0, -0.115), frameon=False, fontsize=7,
          handletextpad=0.3)
fig.tight_layout()
fig.savefig(OUT + 'Figure1.png', bbox_inches='tight')
fig.savefig(OUT + 'Figure1.svg', bbox_inches='tight')
plt.close(fig)

# ============================== FIGURE 2 ==============================
# the two residue groups marked in panel a; the claims are about the residue
# being replaced, so only that panel carries them
IONISABLE = set('DEHKR')
SMALLEST = set('AGS')

fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.5))
for ax, key, lab, panel in [(axes[0], 'aa1', 'substituted from', 'a'),
                            (axes[1], 'aa2', 'substituted to', 'b')]:
    g = cells.groupby(key)['delta_TM_minus_nonTM'].mean().sort_values()
    norm = TwoSlopeNorm(vmin=-0.28, vcenter=0, vmax=0.28)
    for k, (aa, v) in enumerate(g.items()):
        ax.barh(k, v, color=DIV(norm(v)), edgecolor=INK, lw=0.5, height=0.74, zorder=3)
    # group membership as a glyph on the tick label, so it survives greyscale
    if key == 'aa1':
        labels = [('\u25cf ' if aa in IONISABLE else '\u25a0 ' if aa in SMALLEST else '') + aa
                  for aa in g.index]
    else:
        labels = list(g.index)
    ax.set_yticks(range(len(g))); ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(0, color=INK, lw=0.9, zorder=4)
    ax.set_xlabel('mean difference, membrane $-$ soluble', fontsize=8)
    ax.set_title(lab, fontsize=8.5, color=INK, pad=4)
    ax.grid(axis='x', lw=0.4, color=GRID, zorder=0); ax.set_axisbelow(True)
    ax.text(-0.17, 1.045, panel, transform=ax.transAxes, fontsize=10, fontweight='bold')
axes[0].annotate('only A and P are\ncostlier in the membrane',
                 xy=(0.045, 18.6), xytext=(-0.20, 14.6), fontsize=6.8, color=INK,
                 arrowprops=dict(arrowstyle='->', lw=0.7, color=INK,
                                 connectionstyle='arc3,rad=-0.25'))
axes[0].legend(handles=[Line2D([], [], marker='o', ls='', color=INK, ms=3.2,
                              label='ionisable (D, E, H, K, R)'),
                        Line2D([], [], marker='s', ls='', color=INK, ms=3.2,
                              label='three smallest side chains (A, G, S)')],
               loc='upper left', bbox_to_anchor=(-0.02, -0.145), frameon=False,
               fontsize=6.6, handletextpad=0.35, labelspacing=0.35)
fig.tight_layout()
fig.savefig(OUT + 'Figure2.png', bbox_inches='tight')
fig.savefig(OUT + 'Figure2.svg', bbox_inches='tight')
plt.close(fig)

# ============================== FIGURE 3 ==============================
# Two panels, same layout twice, so the comparison survives greyscale.
adj = pd.read_csv(ADJ) if os.path.exists(ADJ) else None
sector = dict(zip(pp.acc, pp['membrane_sector_subunit'])) if 'acc' in pp.columns \
    else dict(zip(pp.iloc[:, 0], pp['membrane_sector_subunit']))

panels = [('a', 'unadjusted',
           {cl: pp.loc[pp['class'] == cl, ['acc', 'cliffs_delta_nonpolar_to_polar']]
                  .rename(columns={'cliffs_delta_nonpolar_to_polar': 'delta'})
            for cl in ('memb_complex', 'free_peptide')})]
if adj is not None:
    panels.append(('b', 'matched on within-protein conservation',
                   {cl: adj.loc[adj.cls == cl, ['acc', 'delta']]
                    for cl in ('memb_complex', 'free_peptide')}))

fig, axes = plt.subplots(len(panels), 1, figsize=(5.2, 2.15 * len(panels)), sharex=True)
axes = np.atleast_1d(axes)
STRATA = [('memb_complex', 'subunits of membrane\ncomplexes', 'o', BLUE),
          ('free_peptide', 'secreted peptides', 's', VERM)]

for ax, (tag, sub, data) in zip(axes, panels):
    rng = np.random.default_rng(0)
    for k, (cl, lab, mk, col) in enumerate(STRATA):
        s = data[cl]
        d = s.delta.values
        y = k + rng.uniform(-0.15, 0.15, len(d))
        sec = np.array([bool(sector.get(a, False)) for a in s.acc])
        ax.plot(d[~sec], y[~sec], mk, ms=4.6, mfc=col, mec='white', mew=0.7,
                ls='', alpha=0.95, zorder=3)
        ax.plot(d[sec], y[sec], mk, ms=7.2, mfc=col, mec=INK, mew=1.4, ls='', zorder=4)
        m = float(np.median(d))
        ax.plot([m, m], [k - 0.29, k + 0.29], '-', color=INK, lw=2.0, zorder=5)
        ax.text(m, k + 0.345, 'median %+.3f' % m, ha='center', fontsize=7, color=INK)
        ax.text(-0.205, k - 0.40, '%d of %d positive' % (int((d > 0).sum()), len(d)),
                ha='left', fontsize=7, color=MUTED)
    ax.axvline(0, color=INK, lw=0.9, ls=(0, (3, 2)), zorder=2)
    ax.set_ylim(-0.62, 1.58); ax.set_xlim(-0.215, 0.445)
    # panel letter bold
    ax.set_title(r'$\bf{%s}$  %s' % (tag, sub), fontsize=8.5, loc='left', color=INK, pad=6)
    ax.grid(axis='x', lw=0.4, color=GRID); ax.set_axisbelow(True)

for ax in axes:
    ax.set_yticks([0, 1]); ax.set_yticklabels([r[1] for r in STRATA], fontsize=7.5)
axes[-1].set_xlabel("per-protein effect size for nonpolar $\\rightarrow$ polar (Cliff's $\\delta$)",
                    fontsize=8)
axes[-1].legend(handles=[Line2D([], [], marker='o', ls='', mfc=BLUE, mec=INK, mew=1.4, ms=6.5,
                               label='membrane-sector subunit carrying no TM helix of its own')],
                loc='upper left', bbox_to_anchor=(0, -0.30), frameon=False, fontsize=7,
                handletextpad=0.4)
fig.tight_layout()
fig.savefig(OUT + 'Figure3.png', bbox_inches='tight')
fig.savefig(OUT + 'Figure3.svg', bbox_inches='tight')
plt.close(fig)

# ============================== FIGURE 5 ==============================
# Figure 5: raw against position-controlled pair correlation.
fig, ax = plt.subplots(figsize=(4.6, 3.6))
ax.axhspan(-0.75, 0, color=VERM, alpha=0.07, zorder=0)
ax.plot([-0.8, 1], [-0.8, 1], ls=(0, (3, 2)), color=MUTED, lw=0.8, zorder=1)
ax.axhline(0, color=INK, lw=0.8, zorder=2)
norm = TwoSlopeNorm(vmin=-0.75, vcenter=0, vmax=0.9)
ax.scatter(pairs.rho_raw, pairs.rho_position_controlled, s=22,
           c=[DIV(norm(v)) for v in pairs.rho_position_controlled],
           edgecolors=INK, linewidths=0.4, zorder=3)
inv = pairs[pairs.rho_position_controlled < -0.45]
for _, r in inv.head(5).iterrows():
    ax.annotate('%s/%s' % (r.target_A, r.target_B),
                (r.rho_raw, r.rho_position_controlled),
                textcoords='offset points', xytext=(7, -2), fontsize=6.8, color=INK)
ax.set_xlabel("Spearman $\\rho$, raw", fontsize=8)
ax.set_ylabel("Spearman $\\rho$, within-position mean removed", fontsize=8)
ax.text(0.035, 0.965,
        '190 target-residue pairs\nevery raw $\\rho$ > 0.63\n109 of 190 invert in sign',
        transform=ax.transAxes, va='top', fontsize=7.2, color=INK)
ax.text(0.97, 0.055, 'sign inverted', transform=ax.transAxes, ha='right',
        fontsize=7, color=VERM, style='italic')
ax.grid(lw=0.4, color=GRID); ax.set_axisbelow(True)
ax.set_ylim(-0.75, 1.0)
fig.tight_layout()
fig.savefig(OUT + 'Figure5.png', bbox_inches='tight')
fig.savefig(OUT + 'Figure5.svg', bbox_inches='tight')
plt.close(fig)

# Figure 4 is not regenerated by this script.
missing = [n for n in (1, 2, 3, 4, 5)
           if not os.path.exists(os.path.join(OUT, 'Figure%d.png' % n))]
print('Figures 1, 2, 3 and 5 written to', OUT)
print('Figure 4 (directional asymmetry) is NOT built here and was left untouched.')
if missing:
    raise SystemExit('MISSING after this run: %s.'
                     % ', '.join('Figure%d.png' % n for n in missing))
