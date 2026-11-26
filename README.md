# Non-transmembrane control for the alpha-helical membrane substitution landscape

Data and code for Karagöl T and Karagöl A. Membrane proteins tolerate substitution more freely than soluble proteins except at alanine positions.

## What the analysis found

Against a verified 66-protein non-transmembrane control, the membrane environment is **more**
permissive to substitution than the soluble environment in 345 of 380 substitution classes
(90.8 percent, paired Wilcoxon p = 6.9e-51). The permissiveness reverses at only two source
residues, alanine and proline, and the six most membrane-penalised classes in the whole matrix
are all alanine to a polar residue. The nonpolar-to-polar penalty reported previously is
therefore not membrane-specific as a general claim, since it holds in the control too, but is
localised to alanine.

The penalty scales with proximity to a membrane rather than requiring insertion into one:
median Cliff's delta +0.238 across 41 soluble subunits of membrane complexes against +0.071
across 25 freely secreted peptides. Per-residue conservation computed for the control set by the
method underlying ConSurf leaves that separation intact at +0.128 (p = 1.4e-5) and raises the
peptide stratum to +0.107 (p = 8.5e-4), so the effect is graded across both strata rather than
confined to one.

Because every value compared is a prediction, the cross-stratum comparison was checked against
observed clinical variants. Of the 216 catalogued human missense variants in the UniProt humsavar
index that fall at scored positions in the control set, 199 carry a benign or pathogenic
classification; AlphaMissense separates them at AUC 0.917 [0.865, 0.962], and equally well in both
strata (0.904 for membrane-complex subunits against 0.941 for secreted peptides, difference +0.037
[-0.050, +0.129], which includes zero). The separation between strata is therefore not an artefact
of the predictor being better calibrated on one of them.

## Layout

```
code/               analysis and figure scripts; the Colab conservation notebook
data/               inputs and derived data
results/            effect tables, conservation QC, the summary of the controlled reanalysis
figures/            Figures 1 to 5 as PNG at 400 dpi and as editable SVG
supplementary/      Supplementary Table S1
```

### Every file in this deposit

The table is complete: if a file is here it is listed, and if it is listed it is here.

| file | what it is |
|---|---|
| `code/classify.py` | the control set itself, 66 accessions assigned to stratum, plus the polarity scheme; imported by `extended_analyses.py` |
| `code/extended_analyses.py` | the five supporting analyses: physicochemical regression, permutation test, length matching, packing register, directional asymmetry |
| `code/build_figures.py` | regenerates Figures 1, 2, 3 and 5 as PNG and SVG; Figure 4 is drawn by hand and is not generated |
| `code/ConSurf_control_set.ipynb` | the Colab notebook that computed per-residue conservation, resume-safe, writes to Drive |
| `data/nontm_substitutions.csv.gz` | 306,622 AlphaMissense substitutions over the 66 control proteins, deduplicated, the input to every analysis here |
| `data/uniprot_transmem_check.tsv` | the UniProt TRANSMEM query for all 71 candidates, including the 5 excluded |
| `data/homolog_survey.csv` | homologue retrieval measured per protein before the conservation run |
| `data/consurf_worklist.csv` | the PDB entry, chain, method and resolution chosen for each of the 66 proteins; the notebook reads this |
| `data/humsavar.txt` | UniProt humsavar release 2026_03, the curated human variant index used for the clinical check |
| `data/alanine_packing_register.csv` | small-residue count in the helical register around each alanine position; read by `extended_analyses.py` |
| `data/asymmetry_190_pairs.csv` | forward minus reverse cost for all 190 residue pairs in both environments; read by `extended_analyses.py` |
| `data/per_protein_effects.csv` | per-protein effect size for all 66 proteins, as a plain CSV; the same values appear with gene names in the `Per_protein_effects` sheet of the results workbook |
| `results/TM_vs_nonTM_results.xlsx` | the primary comparison: summary, all 380 substitution classes, per-protein effects, pair correlations |
| `results/RESULTS_conservation_controlled.txt` | the conservation-controlled reanalysis as it was reported, including the two proteins that failed homologue retrieval |
| `results/conservation_qc.csv` | per-protein conservation result: homologues used, alignment source, depth flag |
| `results/effects_all_raw.csv`, `results/effects_all_adj.csv` | per-protein effect sizes over the 64 proteins that completed conservation, unadjusted and adjusted |
| `results/clinical_variants_control_set.csv` | the 216 catalogued variants intersecting the control set, with classification and score |
| `figures/Figure1.png`, `figures/Figure1.svg` | distribution of the per-class difference across all 380 classes |
| `figures/Figure2.png`, `figures/Figure2.svg` | mean difference by source residue, and the alanine row in full |
| `figures/Figure3.png`, `figures/Figure3.svg` | per-protein effect sizes by stratum, unadjusted and conservation-adjusted |
| `figures/Figure4.png`, `figures/Figure4.svg` | the control set design, drawn by hand rather than generated |
| `figures/Figure5.png`, `figures/Figure5.svg` | the clinical-variant validation |
| `supplementary/Supplementary_Table_S1.xlsx` | the 66 proteins, the 5 excluded candidates, the 216 clinical variants, and notes |
| `README.md` | this file: what the analysis found, the layout, and the complete file list |
| `CITATION.cff` | how to cite this deposit |
| `LICENSE` | MIT, covering the code |
| `LICENSE-DATA` | CC BY 4.0, covering the derived data |

## Reproducing

```bash
python3 code/extended_analyses.py            # the five supporting analyses
python3 code/build_figures.py                # regenerate Figures 1, 2, 3 and 5
```

Every script runs on the files in this repository with no arguments. Figure 4 is not regenerated
by `build_figures.py`.

Two large outputs are available on request rather than included here: the per-residue
conservation master table (`All_Controls_Conservation_Master.tsv`, one row per scored residue)
and the 66 per-protein homologue alignments.

## Data provenance and licensing

Substitution scores are AlphaMissense predictions from Cheng et al. (2023),
*Science* 381(6664):eadg7492, [doi:10.1126/science.adg7492](https://doi.org/10.1126/science.adg7492),
obtained through the AlphaMissense data resource of Tordai et al. (2024), *Scientific Data* 11:495,
[doi:10.1038/s41597-024-03327-8](https://doi.org/10.1038/s41597-024-03327-8), at
[alphamissense.hegelab.org](https://alphamissense.hegelab.org). AlphaMissense predictions are
released by Google DeepMind under **CC BY 4.0**, which permits redistribution with attribution;
the derived tables here are redistributed on those terms.

Control-set assembly and all TRANSMEM feature queries were made against **UniProt release 2026_03**
(published 2 September 2026). Clinical variant classifications are from the UniProt humsavar
index of the same release. UniProt states that humsavar classifications are curated from literature
reports, are intended for research rather than clinical or diagnostic use, and are subject to
revision.

Protein annotations are from UniProt, released under CC BY 4.0. Homologue sets are drawn from
UniRef50 and Pfam. Conservation was computed with
[rate4site](https://www.tau.ac.il/~itaymay/cp/rate4site.html) and alignments with
[MAFFT](https://mafft.cbrc.jp/alignment/software/).

The code in `code/` is released under the MIT licence (`LICENSE`). The derived data in `data/`,
`results/` and `supplementary/` are released under CC BY 4.0 (`LICENSE-DATA`), consistent with
the terms of the sources they are derived from.

## Citing

See `CITATION.cff`. Until the paper appears, cite this repository.

## Contact

Taner Karagöl and Alper Karagöl. Issues and questions via the issue tracker.
