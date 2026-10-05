#!/usr/bin/env python3
"""Recompute the conservation values from the Colab run outputs and assert
each one.

Inputs: conservation_qc.csv, effects_all_raw.csv and effects_all_adj.csv.

    python3 code/verify_conservation.py [--paper-root <repository root>]

Runs 36 checks. Exits non-zero if any fails.
"""
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

ROOT = os.environ.get('PAPER_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# results/ in this repository
C = next((d for d in (os.path.join(ROOT, 'consurf'), os.path.join(ROOT, 'results'))
          if os.path.exists(os.path.join(d, 'conservation_qc.csv'))), None)
if C is None:
    sys.exit("could not find conservation_qc.csv under %s in consurf/ or results/" % ROOT)
CHECKS = []


def ck(name, got, want, tol=0.0):
    ok = abs(got - want) <= tol if isinstance(want, (int, float)) else got == want
    CHECKS.append((name, got, want, ok))
    fmt = (lambda v: "%.4f" % v) if isinstance(want, float) else str
    print("  %-56s got %-11s want %-11s %s"
          % (name, fmt(got), fmt(want), "OK" if ok else "*** MISMATCH ***"))


def strat(d, cls):
    return d[d.cls == cls].delta.values


def main():
    qc = pd.read_csv(os.path.join(C, 'conservation_qc.csv'))
    raw = pd.read_csv(os.path.join(C, 'effects_all_raw.csv'))
    adj = pd.read_csv(os.path.join(C, 'effects_all_adj.csv'))

    print("=== coverage and depth ===")
    ck("proteins scored", len(qc), 64)
    ck("peptides scored", int((qc.stratum == 'free_peptide').sum()), 23)
    ck("subunits scored", int((qc.stratum == 'memb_complex').sum()), 41)
    ck("below the 50-homologue threshold", int(qc.below_floor.sum()), 14)
    ck("  of which peptides", int(qc[qc.stratum == 'free_peptide'].below_floor.sum()), 5)
    ck("  of which subunits", int(qc[qc.stratum == 'memb_complex'].below_floor.sum()), 9)
    for s, n in (('free_peptide', 23), ('memb_complex', 41)):
        sub = qc[qc.stratum == s]
        ck("below-floor share, %s (percent)" % s,
           round(100 * sub.below_floor.sum() / len(sub)), 22)
    ck("alignments from UniRef50 clusters", int((qc.source == 'uniref50').sum()), 45)
    ck("alignments from Pfam family searches", int((qc.source != 'uniref50').sum()), 19)
    ck("above the threshold", int((~qc.below_floor).sum()), 50)

    print("\n=== unadjusted, on the 64 scored proteins ===")
    m, p = strat(raw, 'memb_complex'), strat(raw, 'free_peptide')
    ck("membrane-complex median", float(np.median(m)), 0.238, 0.001)
    ck("membrane-complex positive", int((m > 0).sum()), 39)
    ck("peptide median", float(np.median(p)), 0.031, 0.001)
    ck("peptide positive", int((p > 0).sum()), 13)
    ck("peptide Wilcoxon p", float(wilcoxon(p, np.zeros(len(p)))[1]), 0.29, 0.01)
    ck("difference between strata", float(np.median(m) - np.median(p)), 0.208, 0.001)

    print("\n=== conservation-adjusted ===")
    m, p = strat(adj, 'memb_complex'), strat(adj, 'free_peptide')
    ck("membrane-complex median", float(np.median(m)), 0.235, 0.001)
    ck("membrane-complex positive", int((m > 0).sum()), 39)
    ck("membrane-complex Wilcoxon p (log10)",
       float(np.log10(wilcoxon(m, np.zeros(len(m)))[1])), float(np.log10(9.1e-12)), 0.15)
    ck("peptide median", float(np.median(p)), 0.107, 0.001)
    ck("peptide positive", int((p > 0).sum()), 18)
    ck("peptide Wilcoxon p", float(wilcoxon(p, np.zeros(len(p)))[1]), 8.5e-4, 1e-4)
    ck("difference between strata", float(np.median(m) - np.median(p)), 0.128, 0.001)
    u, pv = mannwhitneyu(m, p, alternative='two-sided')
    ck("between-strata Mann-Whitney p (log10)", float(np.log10(pv)), float(np.log10(1.4e-5)), 0.15)
    ck("between-strata Cliff's delta", 2 * u / (len(m) * len(p)) - 1, 0.659, 0.002)

    print("\n=== the two sensitivity cuts ===")
    sub = adj[~adj.below_floor]
    a, b = strat(sub, 'memb_complex'), strat(sub, 'free_peptide')
    ck("above-threshold n", len(sub), 50)
    ck("above-threshold difference", float(np.median(a) - np.median(b)), 0.138, 0.001)
    ck("above-threshold p (log10)",
       float(np.log10(mannwhitneyu(a, b, alternative='two-sided')[1])),
       float(np.log10(2.3e-4)), 0.15)
    uni = adj[adj.acc.isin(qc[qc.source == 'uniref50'].uniprot)]
    a, b = strat(uni, 'memb_complex'), strat(uni, 'free_peptide')
    ck("uniref50-only n", len(uni), 45)
    ck("uniref50-only difference", float(np.median(a) - np.median(b)), 0.127, 0.001)
    ck("uniref50-only p (log10)",
       float(np.log10(mannwhitneyu(a, b, alternative='two-sided')[1])),
       float(np.log10(2.7e-4)), 0.15)

    print("\n=== claims about individual proteins ===")
    ck("CCL20 absent from the scored set", 'P78556' in set(qc.uniprot), False)
    ck("INSL5 absent from the scored set", 'Q9Y5Q6' in set(qc.uniprot), False)
    sector = ['O75964', 'P61421', 'P56385', 'Q9Y3B6']
    s = adj[adj.acc.isin(sector)].delta
    ck("membrane-sector subunits, adjusted minimum", float(s.min()), 0.186, 0.001)
    ck("membrane-sector subunits, adjusted maximum", float(s.max()), 0.250, 0.001)

    bad = [c for c in CHECKS if not c[3]]
    print("\n" + "=" * 74)
    print("%d checks, %d passed, %d FAILED" % (len(CHECKS), len(CHECKS) - len(bad), len(bad)))
    for c in bad:
        print("   FAILED:", c[0], "got", c[1], "want", c[2])
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
