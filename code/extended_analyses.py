#!/usr/bin/env python3
"""Analyses kept separate from the core pipeline so each can be rerun on its own:

  1. Physicochemical regression   is the alanine effect reducible to volume or hydropathy?
  2. Permutation test             is the stratum difference real, with labels shuffled?
  3. Bootstrap intervals          confidence intervals on each effect size
  4. Packing-register test        register composition against substitution cost
  5. Directional asymmetry        forward minus reverse cost, both environments

    python3 code/extended_analyses.py [--paper-root <repository root>]

Inputs, both included in this repository: data/nontm_substitutions.csv.gz and
results/TM_vs_nonTM_results.xlsx. Writes alanine_packing_register.csv and
asymmetry_190_pairs.csv to data/.
"""
import argparse, os, sys
import numpy as np
import pandas as pd
import numpy.linalg as la
from scipy.stats import t as tdist, mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classify import FREE_PEPTIDE, MEMB_COMPLEX, NONPOLAR, POLAR

# Kyte-Doolittle hydropathy; high is hydrophobic
KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4,
      'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8,
      'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
# Zamyatnin residue volumes, cubic angstroms
VOL = {'A': 88.6, 'R': 173.4, 'N': 114.1, 'D': 111.1, 'C': 108.5, 'Q': 143.8, 'E': 138.4,
       'G': 60.1, 'H': 153.2, 'I': 166.7, 'L': 166.7, 'K': 168.6, 'M': 162.9, 'F': 189.9,
       'P': 112.7, 'S': 89.0, 'T': 116.1, 'W': 227.8, 'Y': 193.6, 'V': 140.0}
CHARGED = set('DEKR')
ROOT = None


def find(name, *dirs):
    """Locate an input under any of the candidate directories."""
    for d in dirs:
        p = os.path.join(ROOT, d, name)
        if os.path.exists(p):
            return p
    sys.exit('cannot find %s under any of %s (root %s)'
             % (name, ', '.join(dirs), ROOT))
SMALL = set('AGS')


def ols(cols, names, y, label):
    X = np.column_stack([np.ones(len(y))] + cols)
    b, _, _, _ = la.lstsq(X, y, rcond=None)
    r = y - X @ b
    n, k = X.shape
    se = np.sqrt(np.diag((r @ r / (n - k)) * la.inv(X.T @ X)))
    p = 2 * (1 - tdist.cdf(np.abs(b / se), n - k))
    r2 = 1 - (r @ r) / ((y - y.mean()) @ (y - y.mean()))
    print(label)
    print("   %-30s %9s %8s %11s" % ("term", "coef", "SE", "p"))
    for nm, bb, ss, pp in zip(['intercept'] + names, b, se, p):
        print("   %-30s %+9.4f %8.4f %11.3g" % (nm, bb, ss, pp))
    print("   R2 = %.3f, n = %d\n" % (r2, n))
    return b, p, r2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--paper-root', default=os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
    a = ap.parse_args()

    global ROOT, MATRIX, RESULTS, OUT
    ROOT = a.paper_root
    MATRIX = find('nontm_substitutions.csv.gz', 'data')
    RESULTS = find('TM_vs_nonTM_results.xlsx', 'results', 'tables')
    OUT = os.path.join(ROOT, 'data')
    rng = np.random.default_rng(20260930)
    c = pd.read_excel(RESULTS, sheet_name='Cells_TM_vs_nonTM')
    pp = pd.read_excel(RESULTS, sheet_name='Per_protein_effects')

    # ---------- 1. physicochemical regression ----------
    print("=" * 74)
    print("1. IS THE ALANINE EFFECT REDUCIBLE TO VOLUME OR HYDROPHOBICITY?\n")
    c['dKD'] = c.aa2.map(KD) - c.aa1.map(KD)
    c['dVol'] = c.aa2.map(VOL) - c.aa1.map(VOL)
    c['srcA'] = (c.aa1 == 'A').astype(float)
    c['srcCh'] = c.aa1.isin(CHARGED).astype(float)
    y = c.delta_TM_minus_nonTM.values
    ols([c.dKD.values, c.dVol.values], ['change in hydropathy (KD)', 'change in volume'], y,
        "Model 1, physicochemistry alone")
    ols([c.dKD.values, c.dVol.values, c.srcA.values],
        ['change in hydropathy (KD)', 'change in volume', 'source is alanine'], y,
        "Model 2, adding an alanine-source indicator")
    b3, p3, _ = ols([c.dKD.values, c.dVol.values, c.srcA.values, c.srcCh.values],
                    ['change in hydropathy (KD)', 'change in volume', 'source is alanine',
                     'source is D/E/K/R'], y, "Model 3, adding a charged-source indicator")
    print("   alanine term after full adjustment: %+.4f, p = %.2g\n" % (b3[3], p3[3]))

    # ---------- 2. permutation ----------
    print("=" * 74)
    print("2. PERMUTATION TEST ON THE STRATUM DIFFERENCE\n")
    mem = pp[pp['class'] == 'memb_complex'].cliffs_delta_nonpolar_to_polar.values
    pep = pp[pp['class'] == 'free_peptide'].cliffs_delta_nonpolar_to_polar.values
    obs = np.median(mem) - np.median(pep)
    allv = np.concatenate([mem, pep])
    null = np.array([(lambda q: np.median(q[:len(mem)]) - np.median(q[len(mem):]))(rng.permutation(allv))
                     for _ in range(10000)])
    hits = int(np.sum(np.abs(null) >= abs(obs)))
    print("   observed %+.4f | null mean %+.4f sd %.4f | %d of 10,000 at least as extreme | p = %.4f\n"
          % (obs, null.mean(), null.std(), hits, (hits + 1) / 10001))

    # ---------- 3. bootstrap ----------
    print("=" * 74)
    print("3. BOOTSTRAP CONFIDENCE INTERVALS, 10,000 RESAMPLES\n")
    def boot(v, f=np.median, B=10000):
        s = np.array([f(rng.choice(v, len(v), replace=True)) for _ in range(B)])
        return f(v), np.percentile(s, 2.5), np.percentile(s, 97.5)
    for lab, v in [('membrane-complex subunits', mem), ('secreted peptides', pep)]:
        m, lo, hi = boot(v)
        print("   %-28s %+.4f  95%% CI [%+.4f, %+.4f]" % (lab, m, lo, hi))
    dd = np.array([np.median(rng.choice(mem, len(mem), True)) - np.median(rng.choice(pep, len(pep), True))
                   for _ in range(10000)])
    print("   %-28s %+.4f  95%% CI [%+.4f, %+.4f]" % ('difference between strata', obs,
          np.percentile(dd, 2.5), np.percentile(dd, 97.5)))
    m, lo, hi = boot(c.delta_TM_minus_nonTM.values)
    print("   %-28s %+.4f  95%% CI [%+.4f, %+.4f]\n" % ('median over 380 classes', m, lo, hi))

    # ---------- 4. packing register ----------
    print("=" * 74)
    print("4. PACKING-REGISTER TEST  (prediction set in advance; result came back OPPOSITE)\n")
    d = pd.read_csv(MATRIX).rename(columns={'acc': 'uniprot_acc',
                                            'entry': 'entry_name'})
    d = d.dropna(subset=['score', 'position']).copy()
    d['position'] = d['position'].astype(int)
    d['cls'] = np.where(d.uniprot_acc.isin(FREE_PEPTIDE), 'pep', 'mem')
    seqs = {}
    for acc, g in d.groupby('uniprot_acc'):
        wt = g.drop_duplicates('position').set_index('position')['aa1'].to_dict()
        seqs[acc] = ''.join(wt.get(i, 'X') for i in range(1, int(max(wt)) + 1))
    rows = []
    for acc, g in d.groupby('uniprot_acc'):
        s = seqs[acc]
        for pos, gg in g[(g.aa1 == 'A') & (g.aa2.isin(POLAR))].groupby('position'):
            i = int(pos) - 1
            nb = [s[i + o] for o in (-4, -3, 3, 4) if 0 <= i + o < len(s)]
            if len(nb) < 4:
                continue
            rows.append({'acc': acc, 'cls': g.cls.iloc[0], 'pos': int(pos),
                         'n_small': sum(1 for x in nb if x in SMALL), 'score': gg.score.median()})
    r = pd.DataFrame(rows)
    for cls, lab in [('mem', 'membrane-complex'), ('pep', 'secreted peptide')]:
        s = r[r.cls == cls]
        rho, p = spearmanr(s.n_small, s.score)
        lo, hi = s[s.n_small <= 1].score, s[s.n_small >= 3].score
        u, pm = mannwhitneyu(hi, lo, alternative='two-sided')
        print("   %-18s n=%4d  rho=%+.3f p=%.2g  |  low %.3f vs high %.3f  delta %+.3f p=%.3g"
              % (lab, len(s), rho, p, lo.median(), hi.median(),
                 2 * u / (len(hi) * len(lo)) - 1, pm))
    r.to_csv(os.path.join(OUT, 'alanine_packing_register.csv'), index=False)
    print("   -> prediction was HIGHER cost with more small residues; observed LOWER. See module docstring.\n")

    # ---------- 5. directional asymmetry ----------
    print("=" * 74)
    print("5. DIRECTIONAL ASYMMETRY\n")
    tm = {(x.aa1, x.aa2): x.tm_median for x in c.itertuples()}
    nt = {(x.aa1, x.aa2): x.nontm_median for x in c.itertuples()}
    AA = sorted({k[0] for k in tm})
    rows = []
    for i, x in enumerate(AA):
        for yy in AA[i + 1:]:
            if (x, yy) in tm and (yy, x) in tm:
                rows.append({'X': x, 'Y': yy,
                             'tm_asym': tm[(x, yy)] - tm[(yy, x)],
                             'nt_asym': nt[(x, yy)] - nt[(yy, x)]})
    aa = pd.DataFrame(rows)
    aa['abs_tm'] = aa.tm_asym.abs()
    aa['abs_nt'] = aa.nt_asym.abs()
    from scipy.stats import wilcoxon
    _, pw = wilcoxon(aa.abs_tm, aa.abs_nt)
    print("   mean |asymmetry| membrane %.4f vs soluble %.4f, larger in membrane for %d of %d, p = %.2e"
          % (aa.abs_tm.mean(), aa.abs_nt.mean(), int((aa.abs_tm > aa.abs_nt).sum()), len(aa), pw))
    print("   agreement between environments: Spearman rho = %.3f"
          % spearmanr(aa.tm_asym, aa.nt_asym)[0])
    aa.to_csv(os.path.join(OUT, 'asymmetry_190_pairs.csv'), index=False)
    print("\nwritten to %s: alanine_packing_register.csv, asymmetry_190_pairs.csv" % OUT)


if __name__ == '__main__':
    main()
