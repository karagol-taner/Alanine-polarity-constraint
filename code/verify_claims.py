#!/usr/bin/env python3
"""Recompute the reported values from the control-set substitution table and
assert each one.

Inputs, both included in this repository: data/nontm_substitutions.csv.gz and
results/TM_vs_nonTM_results.xlsx.

    python3 code/verify_claims.py [--paper-root <repository root>]

Runs 46 checks. Exits non-zero if any fails.
"""
import argparse, os, sys
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classify import FREE_PEPTIDE, MEMB_COMPLEX, NONPOLAR, POLAR

# Kyte-Doolittle hydropathy and Zamyatnin residue volumes, as used in
# code/extended_analyses.py
KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4,
      'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8,
      'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
VOL = {'A': 88.6, 'R': 173.4, 'N': 114.1, 'D': 111.1, 'C': 108.5, 'Q': 143.8, 'E': 138.4,
       'G': 60.1, 'H': 153.2, 'I': 166.7, 'L': 166.7, 'K': 168.6, 'M': 162.9, 'F': 189.9,
       'P': 112.7, 'S': 89.0, 'T': 116.1, 'W': 227.8, 'Y': 193.6, 'V': 140.0}

CHECKS = []
ROOT = None


def find(name, *dirs):
    """Locate an input under any of the candidate directories."""
    for d in dirs:
        p = os.path.join(ROOT, d, name)
        if os.path.exists(p):
            return p
    sys.exit('cannot find %s under any of %s (root %s)'
             % (name, ', '.join(dirs), ROOT))


def ck(name, got, want, tol=0.0):
    if isinstance(got, bool) or isinstance(want, bool):
        ok = got == want
    else:
        ok = abs(got - want) <= tol
    CHECKS.append((name, got, want, ok))
    g = ("%.4f" % got) if isinstance(got, float) else str(got)
    w = ("%.4f" % want) if isinstance(want, float) else str(want)
    print("  %-52s got %-12s want %-12s %s" % (name, g, w, "OK" if ok else "*** MISMATCH ***"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--paper-root', default=os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
    a = ap.parse_args()

    global ROOT
    ROOT = a.paper_root
    matrix = find('nontm_substitutions.csv.gz', 'data')
    results = find('TM_vs_nonTM_results.xlsx', 'results', 'tables')

    d = pd.read_csv(matrix)
    d = d.rename(columns={'acc': 'uniprot_acc', 'entry': 'entry_name'})
    d = d.dropna(subset=['score', 'position'])
    d['position'] = d['position'].astype(int)
    d['cls'] = np.where(d.uniprot_acc.isin(FREE_PEPTIDE), 'pep',
                np.where(d.uniprot_acc.isin(MEMB_COMPLEX), 'mem', None))

    print("=== counts ===")
    ck("total substitutions", len(d), 306622)
    ck("proteins", d.uniprot_acc.nunique(), 66)
    ck("free_peptide substitutions", int((d.cls == 'pep').sum()), 68571)

    print("\n=== global membrane vs soluble comparison ===")
    cells = pd.read_excel(results, sheet_name='Cells_TM_vs_nonTM')
    tm = {(r.aa1, r.aa2): float(r.tm_median) for r in cells.itertuples()}
    nt = d.groupby(['aa1', 'aa2'])['score'].median()
    common = [k for k in tm if k in nt.index]
    tmv = np.array([tm[k] for k in common])
    ntv = np.array([nt[k] for k in common])
    ck("comparable classes", len(common), 380)
    ck("classes where membrane is lower", int((tmv < ntv).sum()), 345)
    _, p = wilcoxon(tmv, ntv)
    ck("paired Wilcoxon p (log10)", float(np.log10(p)), float(np.log10(6.9e-51)), 0.3)

    print("\n=== source-residue ordering ===")
    delta = pd.DataFrame({'aa1': [k[0] for k in common], 'aa2': [k[1] for k in common], 'd': tmv - ntv})
    src = delta.groupby('aa1')['d'].mean()
    for aa, want in [('H', -0.2578), ('R', -0.2350), ('D', -0.2176), ('A', 0.0531), ('P', 0.0177)]:
        ck("mean difference, source %s" % aa, float(src[aa]), want, 0.001)
    ck("source residues with a positive mean", int((src > 0).sum()), 2)
    top6 = delta.nlargest(6, 'd')
    ck("six most membrane-penalised are all alanine", int((top6.aa1 == 'A').sum()), 6)
    toP = delta[delta.aa2 == 'P']
    ck("X->P cheaper in membrane, of 19", int((toP.d < 0).sum()), 19)

    ION = set('DEKRH')
    io = src[[a for a in src.index if a in ION]]
    rest = src[[a for a in src.index if a not in ION]]
    ck("ionisable source residues, mean", float(io.mean()), -0.2008, 0.001)
    ck("other fifteen source residues, mean", float(rest.mean()), -0.0666, 0.001)
    ck("ionisable vs rest, Mann-Whitney p (log10)",
       float(np.log10(mannwhitneyu(io, rest).pvalue)), float(np.log10(5.16e-4)), 0.1)
    rank = src.rank()
    for aa, pos, want in [('G', 8, -0.1325), ('S', 14, -0.0534), ('A', 20, 0.0531)]:
        ck("smallest residues, %s mean" % aa, float(src[aa]), want, 0.001)
        ck("smallest residues, %s rank of 20" % aa, int(rank[aa]), pos)
    ck("ordering vs hydropathy, Spearman rho",
       float(spearmanr([KD[a] for a in src.index], src.values)[0]), 0.303, 0.002)
    ck("ordering vs hydropathy, not significant",
       float(spearmanr([KD[a] for a in src.index], src.values)[1]) > 0.05, True)

    print("\n=== the alanine row on its own ===")
    arow = delta[delta.aa1 == 'A']
    apol = arow[arow.aa2.isin(POLAR)]['d']
    anon = arow[~arow.aa2.isin(POLAR)]['d']
    ck("A->polar classes", len(apol), 10)
    ck("A->nonpolar classes", len(anon), 9)
    ck("A->polar mean", float(apol.mean()), 0.1103, 0.001)
    ck("A->nonpolar mean", float(anon.mean()), -0.0105, 0.001)
    ck("A row, polar vs nonpolar p (log10)",
       float(np.log10(mannwhitneyu(apol, anon).pvalue)), float(np.log10(1.67e-3)), 0.1)
    dvol = [VOL[x] - VOL['A'] for x in arow.aa2]
    ck("A row vs volume increase, Spearman rho",
       float(spearmanr(dvol, arow.d.values)[0]), -0.042, 0.002)
    for tg, want in (('W', 0.0142), ('D', 0.1071)):
        ck("A->%s difference" % tg, float(arow[arow.aa2 == tg].d.iloc[0]), want, 0.001)

    print("\n=== per-protein effect sizes ===")
    d['dir'] = np.where(d.aa1.isin(NONPOLAR) & d.aa2.isin(POLAR), 'np',
                np.where(d.aa1.isin(POLAR) & d.aa2.isin(NONPOLAR), 'pn', 'x'))
    rows = []
    for acc, g in d.groupby('uniprot_acc'):
        x = g[g['dir'] == 'np']['score']
        y = g[g['dir'] == 'pn']['score']
        if len(x) >= 20 and len(y) >= 20:
            u, _ = mannwhitneyu(x, y)
            rows.append({'cls': g.cls.iloc[0], 'dl': 2 * u / (len(x) * len(y)) - 1, 'len': g.position.max()})
    pp = pd.DataFrame(rows)
    for cl, npos, med in [('mem', 39, 0.2382), ('pep', 15, 0.0714)]:
        s = pp[pp.cls == cl]
        ck("%s proteins with a positive effect" % cl, int((s.dl > 0).sum()), npos)
        ck("%s median Cliff's delta" % cl, float(s.dl.median()), med, 0.001)
    sub = pp[(pp['len'] >= 100) & (pp['len'] <= 400)]
    ck("length-matched mem n", len(sub[sub.cls == 'mem']), 28)
    ck("length-matched pep n", len(sub[sub.cls == 'pep']), 21)

    print("\n=== pair correlations ===")
    d['pm'] = d.groupby(['entry_name', 'position'])['score'].transform('mean')
    d['r'] = d.score - d.pm
    raws, ctrl = [], []
    tg = sorted(d.aa2.unique())
    for i, x in enumerate(tg):
        for y in tg[i + 1:]:
            A = d[d.aa2 == x][['entry_name', 'position', 'score', 'r']]
            B = d[d.aa2 == y][['entry_name', 'position', 'score', 'r']]
            j = A.merge(B, on=['entry_name', 'position'], suffixes=('_x', '_y'))
            if len(j) < 100:
                continue
            raws.append(spearmanr(j.score_x, j.score_y)[0])
            ctrl.append(spearmanr(j.r_x, j.r_y)[0])
    ck("pairs evaluated", len(raws), 190)
    ck("minimum raw rho", float(min(raws)), 0.634, 0.002)
    ck("pairs inverting in sign", int(sum(1 for v in ctrl if v < 0)), 109)

    print("\n=== directional asymmetry ===")
    asym_tm, asym_nt = [], []
    AAs = sorted({k[0] for k in tm})
    for i, x in enumerate(AAs):
        for y in AAs[i + 1:]:
            if (x, y) in tm and (y, x) in tm and (x, y) in nt.index and (y, x) in nt.index:
                asym_tm.append(tm[(x, y)] - tm[(y, x)])
                asym_nt.append(nt[(x, y)] - nt[(y, x)])
    asym_tm, asym_nt = np.array(asym_tm), np.array(asym_nt)
    ck("asymmetry pairs", len(asym_tm), 190)
    ck("mean |asymmetry|, membrane", float(np.abs(asym_tm).mean()), 0.2696, 0.001)
    ck("mean |asymmetry|, soluble", float(np.abs(asym_nt).mean()), 0.2397, 0.001)
    ck("asymmetry agreement between environments (rho)",
       float(spearmanr(asym_tm, asym_nt)[0]), 0.921, 0.002)

    bad = [c for c in CHECKS if not c[3]]
    print("\n" + "=" * 72)
    print("%d checks, %d passed, %d FAILED" % (len(CHECKS), len(CHECKS) - len(bad), len(bad)))
    for c in bad:
        print("   FAILED:", c[0], "got", c[1], "want", c[2])
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
