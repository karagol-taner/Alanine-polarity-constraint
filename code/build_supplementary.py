#!/usr/bin/env python3
"""Build Supplementary Table S1.

    python3 code/build_supplementary.py

Writes supplementary/Supplementary_Table_S1.xlsx with four sheets: the 66 control
proteins, the 5 excluded candidates, the clinical variants, and notes.
"""
import os, sys
import pandas as pd
import numpy as np

ROOT = os.environ.get('PAPER_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, 'supplementary')
os.makedirs(OUT, exist_ok=True)

def find(name, *dirs):
    """Locate an input under any of the candidate directories."""
    for d in dirs:
        p = os.path.join(ROOT, d, name)
        if os.path.exists(p):
            return p
    sys.exit('cannot find %s under any of %s (root %s)'
             % (name, ', '.join(dirs), ROOT))


sub = pd.read_csv(find('nontm_substitutions.csv.gz', 'data'))
tm = pd.read_csv(find('uniprot_transmem_check.tsv', 'data'), sep='\t')
qc = pd.read_csv(find('conservation_qc.csv', 'consurf', 'results'))
pp = pd.read_excel(find('TM_vs_nonTM_results.xlsx', 'tables', 'results'),
                   sheet_name='Per_protein_effects')
adj = pd.read_csv(find('effects_all_adj.csv', 'consurf', 'results'))

STRAT = {'free_peptide': 'freely secreted peptide',
         'memb_complex': 'soluble subunit of a membrane complex'}

# ---- sheet 1: the 66 retained proteins -------------------------------------
g = sub.groupby('acc').agg(entry=('entry', 'first'), gene=('gene', 'first'),
                           stratum=('class', 'first'), sector=('sector', 'first'),
                           substitutions=('score', 'size'),
                           positions=('position', 'nunique')).reset_index()
g['length'] = g.acc.map(dict(zip(tm.Entry, tm.Length)))
g['protein_name'] = g.acc.map(dict(zip(tm.Entry, tm['Protein names'])))
g['stratum'] = g.stratum.map(STRAT)
g['membrane_sector_subunit'] = g.acc.isin(pp.loc[pp['membrane_sector_subunit'].astype(bool), 'acc'])
g['cliffs_delta_unadjusted'] = g.acc.map(dict(zip(pp.acc, pp.cliffs_delta_nonpolar_to_polar)))

q = qc.set_index('uniprot')
for col, new in (('homologs', 'conservation_homologues'),
                 ('source', 'conservation_alignment_source'),
                 ('below_floor', 'below_50_homologue_threshold'),
                 ('mean_identity_to_query', 'mean_identity_to_query'),
                 ('tree_length', 'alignment_tree_length')):
    g[new] = g.acc.map(q[col])
g['cliffs_delta_conservation_adjusted'] = g.acc.map(dict(zip(adj.acc, adj.delta)))
g['conservation_scored'] = g.acc.isin(qc.uniprot)

g = g[['acc', 'entry', 'gene', 'protein_name', 'stratum', 'membrane_sector_subunit',
       'length', 'positions', 'substitutions', 'cliffs_delta_unadjusted',
       'conservation_scored', 'conservation_homologues', 'conservation_alignment_source',
       'below_50_homologue_threshold', 'mean_identity_to_query', 'alignment_tree_length',
       'cliffs_delta_conservation_adjusted']].sort_values(['stratum', 'acc'])
g.columns = ['UniProt', 'Entry name', 'Gene', 'Protein name', 'Stratum',
             'Membrane-sector subunit', 'Length (aa)', 'Positions scored', 'Substitutions',
             "Cliff's delta, unadjusted", 'Conservation computed', 'Homologues used',
             'Alignment source', 'Below 50-homologue threshold', 'Mean identity to query',
             'Tree length', "Cliff's delta, conservation-adjusted"]

# ---- sheet 2: the excluded candidates --------------------------------------
excl = tm[~tm.Entry.isin(g['UniProt'])].copy()
excl = excl[['Entry', 'Gene Names (primary)', 'Protein names', 'Length',
             'Transmembrane', 'Subcellular location [CC]']]
excl.columns = ['UniProt', 'Gene', 'Protein name', 'Length (aa)',
                'TRANSMEM feature (reason for exclusion)', 'Subcellular location']

# ---- sheet 3: clinical variants intersecting the control set ----------------
# parsed from humsavar.txt
import re as _re

AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
       'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
       'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}
_ROW = _re.compile(r'^(\S+)\s+([A-Z0-9]+)\s+(VAR_\d+)\s+'
                   r'p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})\s+(LB/B|LP/P|US)\b')
SIGLAB = {'LP/P': 'pathogenic or likely pathogenic',
          'LB/B': 'benign or likely benign',
          'US': 'uncertain significance'}

hv_path = find('humsavar.txt', 'data')
_rows = []
for _line in open(hv_path, encoding='utf-8', errors='replace'):
    _m = _ROW.match(_line)
    if not _m:
        continue
    _g, _acc, _vid, _a1, _pos, _a2, _sig = _m.groups()
    _a1, _a2 = AA3.get(_a1.upper()), AA3.get(_a2.upper())
    if not _a1 or not _a2 or _a1 == _a2:
        continue
    _rows.append((_g, _acc, int(_pos), _a1, _a2, _sig, _vid))
hv = pd.DataFrame(_rows, columns=['gene', 'acc', 'position', 'aa1', 'aa2',
                                  'sig', 'variant_id'])
# one row per substitution: humsavar lists a variant once per disease
hv = hv.drop_duplicates(subset=['acc', 'position', 'aa1', 'aa2'])

cv = hv.merge(sub[['acc', 'position', 'aa1', 'aa2', 'score', 'class']],
              on=['acc', 'position', 'aa1', 'aa2'], how='inner')
cv = cv.drop_duplicates(subset=['acc', 'position', 'aa1', 'aa2'])
cv['substitution'] = cv.aa1 + cv.position.astype(str) + cv.aa2
cv['classification'] = cv.sig.map(SIGLAB)
cv['stratum'] = cv['class'].map(STRAT)
cv = cv[['acc', 'gene', 'variant_id', 'substitution', 'aa1', 'position', 'aa2',
         'classification', 'stratum', 'score']].sort_values(['gene', 'position'])
cv.columns = ['UniProt', 'Gene', 'UniProt variant ID', 'Substitution',
              'Reference residue', 'Position', 'Target residue',
              'Clinical classification', 'Stratum', 'AlphaMissense score']

# ---- sheet 4: notes ---------------------------------------------------------
notes = pd.DataFrame({'Supplementary Table S1': [
    'Sheet "Control set": the 66 human proteins retained as the non-transmembrane control, '
    'each verified against UniProt TRANSMEM annotation rather than accepted by prior categorisation.',
    'Sheet "Excluded candidates": candidates carrying at least one annotated TRANSMEM feature, '
    'with the feature that excluded them.',
    'Sheet "Clinical variants": the catalogued human missense variants from the UniProt humsavar '
    'index (release 2026_03) that fall at scored positions in the control set, with their curated '
    'classification and predicted score. humsavar lists a variant once for each disease associated '
    'with it; entries here are reduced to distinct substitutions, so three substitutions that appear '
    'on two lines each in the source appear once here. Variants of uncertain significance are '
    'listed but were excluded from the discrimination analysis.',
    '',
    "Cliff's delta is the per-protein effect size for nonpolar-to-polar substitution against "
    "polar-to-nonpolar substitution in the same protein, computed independently within each protein.",
    'The conservation columns come from the per-residue run described in Methods. Two proteins, '
    'CCL20 (P78556) and INSL5 (Q9Y5Q6), yielded too few homologues inside the identity window to '
    'score and are marked Conservation computed = False.',
    'Alignment source "uniref50" means the UniRef50 cluster containing the query sufficed; '
    '"pfam:PFxxxxx:sp" means a Pfam family search restricted to reviewed entries was needed, '
    'and ":all" that unreviewed entries were also required.',
    'Below 50-homologue threshold marks proteins beneath the depth at which ConSurf warns against '
    'interpretation. These are retained and flagged rather than dropped, and the analysis is '
    'reported with and without them.',
    'Mean identity to query and tree length are descriptive diagnostics of alignment divergence. '
    'They depend on which homologues were retrieved and are not comparable across proteins drawn '
    'from different sources.',
    'The humsavar classification is curated from literature reports and is stated by UniProt to be '
    'for research use rather than clinical or diagnostic use, and to be subject to revision.',
]})

p = os.path.join(OUT, 'Supplementary_Table_S1.xlsx')
with pd.ExcelWriter(p, engine='openpyxl') as w:
    notes.to_excel(w, sheet_name='Notes', index=False)
    g.to_excel(w, sheet_name='Control set', index=False)
    excl.to_excel(w, sheet_name='Excluded candidates', index=False)
    cv.to_excel(w, sheet_name='Clinical variants', index=False)
    for name, df in (('Notes', notes), ('Control set', g),
                     ('Excluded candidates', excl), ('Clinical variants', cv)):
        ws = w.sheets[name]
        for i, col in enumerate(df.columns, 1):
            width = max(len(str(col)), *(len(str(v)) for v in df[col].head(80))) + 2
            ws.column_dimensions[ws.cell(1, i).column_letter].width = min(width, 60)
        ws.freeze_panes = 'A2'

print("wrote %s" % p)
print("  control set      %d rows, %d columns" % g.shape)
print("  excluded         %d rows" % len(excl))
print("  conservation scored %d of %d" % (int(g['Conservation computed'].sum()), len(g)))
print("  clinical variants %d rows (%d classified, %d uncertain)"
      % (len(cv), int((cv['Clinical classification'] != 'uncertain significance').sum()),
         int((cv['Clinical classification'] == 'uncertain significance').sum())))
