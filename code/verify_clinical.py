#!/usr/bin/env python3
"""Recompute the clinical-variant values from source and assert each one.

Inputs: data/humsavar.txt and data/nontm_substitutions.csv.gz.

    python3 code/verify_clinical.py [--paper-root <repository root>]

Runs 25 checks. Each prints PASS or FAIL; exits non-zero if any fails.
"""
import argparse
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd

AP = argparse.ArgumentParser()
AP.add_argument('--paper-root', default=os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))))
AP.add_argument('--humsavar', default=None,
                help='path to humsavar.txt; defaults to <paper-root>/data/humsavar.txt')
A = AP.parse_args()

ROOT = A.paper_root
HUMSAVAR = A.humsavar or os.path.join(ROOT, 'data', 'humsavar.txt')

checks = []


def check(name, got, want, tol=None):
    if tol is None:
        ok = got == want
    else:
        ok = abs(float(got) - float(want)) <= tol
    checks.append((name, ok, got, want))
    print('%-4s %-62s got=%-22s want=%s'
          % ('PASS' if ok else 'FAIL', name, got, want))
    return ok


# ---------------------------------------------------------------- load sources
if not os.path.exists(HUMSAVAR):
    sys.exit('missing humsavar: %s' % HUMSAVAR)

AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
       'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
       'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}

ROW = re.compile(
    r'^(\S+)\s+([A-Z0-9]+)\s+(VAR_\d+)\s+p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})\s+'
    r'(LB/B|LP/P|US)\b')

rows = []
with open(HUMSAVAR, encoding='utf-8', errors='replace') as fh:
    for line in fh:
        m = ROW.match(line)
        if not m:
            continue
        gene, acc, _vid, a1, pos, a2, sig = m.groups()
        a1, a2 = AA3.get(a1.upper()), AA3.get(a2.upper())
        if not a1 or not a2 or a1 == a2:
            continue
        rows.append((gene, acc, int(pos), a1, a2, sig))

hv = pd.DataFrame(rows, columns=['gene', 'acc', 'pos', 'aa1', 'aa2', 'sig'])
hv = hv.drop_duplicates(subset=['acc', 'pos', 'aa1', 'aa2'])
print('humsavar parsed: %d unique missense variants, %d accessions'
      % (len(hv), hv.acc.nunique()))

mat_path = os.path.join(ROOT, 'data', 'nontm_substitutions.csv.gz')
if not os.path.exists(mat_path):
    sys.exit('missing control-set matrix: %s' % mat_path)
with gzip.open(mat_path, 'rt') as fh:
    mat = pd.read_csv(fh)
cols = {c.lower(): c for c in mat.columns}


def col(*cands):
    for c in cands:
        if c in cols:
            return cols[c]
    sys.exit('matrix lacks any of %s; has %s' % (cands, list(mat.columns)))


C_ACC = col('uniprot_id', 'acc', 'accession', 'uniprot')
C_POS = col('protein_variant_pos', 'pos', 'position')
C_A1 = col('aa1', 'ref', 'reference_aa', 'wt')
C_A2 = col('aa2', 'alt', 'alternate_aa', 'mut')
C_SCORE = col('am_pathogenicity', 'score', 'alphamissense', 'pathogenicity')
C_STRAT = col('stratum', 'class', 'group', 'category')

mat = mat.rename(columns={C_ACC: 'acc', C_POS: 'pos', C_A1: 'aa1',
                          C_A2: 'aa2', C_SCORE: 'score', C_STRAT: 'stratum'})
mat['pos'] = mat['pos'].astype(int)
print('control-set matrix: %d rows, %d proteins, strata=%s'
      % (len(mat), mat.acc.nunique(), sorted(mat.stratum.unique())))

# ------------------------------------------------------------------- the join
j = hv.merge(mat[['acc', 'pos', 'aa1', 'aa2', 'score', 'stratum']],
             on=['acc', 'pos', 'aa1', 'aa2'], how='inner')
j = j.drop_duplicates(subset=['acc', 'pos', 'aa1', 'aa2'])

n_all = len(j)
n_prot = j.acc.nunique()
lab = j[j.sig.isin(['LP/P', 'LB/B'])].copy()
lab['y'] = (lab.sig == 'LP/P').astype(int)

print()
print('joined: %d variants in %d control-set proteins' % (n_all, n_prot))
print('labelled: %d (%d pathogenic, %d benign); unclassified: %d'
      % (len(lab), int(lab.y.sum()), int((1 - lab.y).sum()),
         int((j.sig == 'US').sum())))
print()


# ------------------------------------------------------------------------ AUC
def auc(y, s):
    y = np.asarray(y)
    s = np.asarray(s, float)
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float('nan')
    r = pd.Series(np.concatenate([pos, neg])).rank().values
    n1, n0 = len(pos), len(neg)
    u = r[:n1].sum() - n1 * (n1 + 1) / 2.0
    return u / (n1 * n0)


def boot_auc(y, s, n=20000, seed=0):
    rng = np.random.default_rng(seed)
    y = np.asarray(y)
    s = np.asarray(s, float)
    out = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(set(y[i].tolist())) < 2:
            continue
        out.append(auc(y[i], s[i]))
    out = np.array(out)
    return np.percentile(out, 2.5), np.percentile(out, 97.5)


def boot_diff(lab, n=20000, seed=0):
    rng = np.random.default_rng(seed)
    g = {k: v.reset_index(drop=True) for k, v in lab.groupby('stratum')}
    a, b = g['free_peptide'], g['memb_complex']
    out = []
    for _ in range(n):
        ia = rng.integers(0, len(a), len(a))
        ib = rng.integers(0, len(b), len(b))
        sa, sb = a.iloc[ia], b.iloc[ib]
        if sa.y.nunique() < 2 or sb.y.nunique() < 2:
            continue
        out.append(auc(sa.y, sa.score) - auc(sb.y, sb.score))
    out = np.array(out)
    return out.mean(), np.percentile(out, 2.5), np.percentile(out, 97.5)


a_all = auc(lab.y, lab.score)
lo_all, hi_all = boot_auc(lab.y, lab.score)

from scipy.stats import mannwhitneyu
u_stat, p_all = mannwhitneyu(lab.score[lab.y == 1], lab.score[lab.y == 0],
                             alternative='greater')

sub = lab[lab.stratum == 'memb_complex']
pep = lab[lab.stratum == 'free_peptide']
a_sub, a_pep = auc(sub.y, sub.score), auc(pep.y, pep.score)
lo_sub, hi_sub = boot_auc(sub.y, sub.score, seed=1)
lo_pep, hi_pep = boot_auc(pep.y, pep.score, seed=2)
d_mean, d_lo, d_hi = boot_diff(lab)

print('AUC all      %.3f [%.3f, %.3f]  n=%d (%d/%d)  MWU p=%.2e'
      % (a_all, lo_all, hi_all, len(lab), int(lab.y.sum()),
         int((1 - lab.y).sum()), p_all))
print('AUC subunits %.3f [%.3f, %.3f]  n=%d (%d/%d)'
      % (a_sub, lo_sub, hi_sub, len(sub), int(sub.y.sum()),
         int((1 - sub.y).sum())))
print('AUC peptides %.3f [%.3f, %.3f]  n=%d (%d/%d)'
      % (a_pep, lo_pep, hi_pep, len(pep), int(pep.y.sum()),
         int((1 - pep.y).sum())))
print('difference (peptides - subunits) %+.3f [%+.3f, %+.3f]'
      % (a_pep - a_sub, d_lo, d_hi))
print()

# alanine-to-polar specifically
POLAR = set('STNQHYCDEKR')
ala = lab[(lab.aa1 == 'A') & (lab.aa2.isin(POLAR))]
print('alanine-to-polar with a clinical label: %d (%d path, %d benign)'
      % (len(ala), int(ala.y.sum()), int((1 - ala.y).sum())))
for _, r in ala.iterrows():
    print('   %-8s %s%s%s  %-5s %-13s %.4f'
          % (r.acc, r.aa1, r.pos, r.aa2, r.sig, r.stratum, r.score))
print()

# ---------------------------------------------------------------------- checks
print('=' * 92)
check('join: total variants matched in control set', n_all, 216)
check('join: distinct proteins carrying a variant', n_prot, 43)
check('labelled variant count (LP/P or LB/B)', len(lab), 199)
check('pathogenic count', int(lab.y.sum()), 143)
check('benign count', int((1 - lab.y).sum()), 56)
check('unclassified (US) count', int((j.sig == 'US').sum()), 17)
check('every matched variant has an AlphaMissense score',
      int(j.score.notna().sum()), n_all)
# humsavar lists a variant once per associated disease, so counting lines rather
# than distinct substitutions would double-count three of them.
check('humsavar deduplicated to physical substitutions',
      int(hv.duplicated(subset=['acc', 'pos', 'aa1', 'aa2']).sum()), 0)
check('control-set matrix has no duplicate substitutions',
      int(mat.duplicated(subset=['acc', 'pos', 'aa1', 'aa2']).sum()), 0)
check('AUC overall', round(a_all, 3), 0.917, tol=0.0015)
check('AUC overall CI low', round(lo_all, 3), 0.865, tol=0.012)
check('AUC overall CI high', round(hi_all, 3), 0.962, tol=0.012)
check('Mann-Whitney p < 1e-15', p_all < 1e-15, True)
check('AUC subunits', round(a_sub, 3), 0.904, tol=0.0015)
check('AUC peptides', round(a_pep, 3), 0.941, tol=0.0015)
check('subunit n (path/benign)',
      (int(sub.y.sum()), int((1 - sub.y).sum())), (96, 26))
check('peptide n (path/benign)',
      (int(pep.y.sum()), int((1 - pep.y).sum())), (47, 30))
check('subunit + peptide = total labelled', len(sub) + len(pep), len(lab))
check('AUC difference point estimate', round(a_pep - a_sub, 3), 0.037,
      tol=0.0015)
check('AUC difference CI includes zero', bool(d_lo < 0 < d_hi), True)
check('AUC difference CI low', round(d_lo, 3), -0.050, tol=0.015)
check('AUC difference CI high', round(d_hi, 3), 0.129, tol=0.015)
check('alanine-to-polar labelled variants (underpowered)', len(ala), 7)
check('both strata exceed AUC 0.85', bool(a_sub > 0.85 and a_pep > 0.85), True)
check('overall AUC CI excludes 0.5 (no-skill)', bool(lo_all > 0.5), True)

print('=' * 92)
nfail = sum(1 for _, ok, _, _ in checks if not ok)
print('%d checks, %d passed, %d failed' % (len(checks), len(checks) - nfail,
                                           nfail))
if nfail:
    sys.exit('VERIFICATION FAILED')
print('All clinical-variant values reproduce from source.')
