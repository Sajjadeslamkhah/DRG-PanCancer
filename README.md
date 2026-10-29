# DRG-PanCancer

Analysis code for the article

> **Disulfidptosis-Related Genes Comprise Two Functionally Distinct, Transcriptionally Uncoupled Subsets with Complementary Prognostic Associations: An Integrated Pan-Cancer Analysis of TCGA Data**
> Zeynep Kus, Sajjad Eslamkhah, Hakan Gulmez, Cuneyd Yavas, Yunus Emre Arvas, Ilknur Yilmaz
> *International Journal of Molecular Sciences* (under review)

The code reproduces every table and figure of the article from public data, starting from the raw files.

## Summary

Disulfidptosis is a regulated cell death caused by NADPH depletion and disulfide stress on the actin cytoskeleton in SLC7A11-high cells. Scores built from 24 disulfidptosis-related genes (DRGs) are widely used as prognostic markers. Analysing 10,242 primary tumours from 33 TCGA cancer types, we show that the DRG set consists of two functionally distinct and transcriptionally uncoupled subsets: a core subset of CRISPR-screen genes (SLC7A11, SLC3A2, RPN1, NCKAP1 and mitochondrial genes) and an actin-cytoskeleton subset (FLNA, FLNB, MYH9, TLN1, ACTB and others). The two subsets have different metabolic expression profiles and carry complementary, cancer-specific prognostic information, so composite DRG scores should be decomposed before biological interpretation.

## Repository contents

| Path | Content |
|---|---|
| `notebooks/01_DRG_pancancer_pipeline.ipynb` | Google Colab notebook that runs the full analysis and draws Figures 1–6 |
| `notebooks/02_DRG_revision_analyses.ipynb` | Colab notebook for the additional analyses performed during peer review |
| `src/drg_pancancer_pipeline.py` | Main pipeline: data download, sample selection, ssGSEA scoring, tumour–normal comparisons, survival analyses, metabolic programme analysis, LASSO-Cox validation, genomic alterations |
| `src/drg_v3_revision.py` | Revision analyses: null models with expression-matched random gene sets, metabolic specificity (oxidative PPP, NRF2 adjustment), tumour-microenvironment adjustment, glioma IDH/1p19q adjustment, proportional-hazards checks, LASSO benchmarks |
| `src/make_main_figures.py` | Assembles manuscript Figures 1–6 from the result tables |
| `src/make_tables.py` | Builds manuscript Tables 1 and 2 (CSV) from the result tables |

## Data

All input data are public and are downloaded automatically from the [UCSC Xena TCGA Pan-Cancer Atlas hub](https://xenabrowser.net) and the [GDC PanCanAtlas publication page](https://gdc.cancer.gov/about-data/publications/pancanatlas):

| Data | File |
|---|---|
| Batch-corrected RNA-seq (EB++) | `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz` |
| Survival endpoints (TCGA-CDR: OS, DSS, DFI, PFI) | `Survival_SupplementalTable_S1_20171025_xena_sp` |
| Sample phenotypes | `TCGA_phenotype_denseDataOnlyDownload.tsv.gz` |
| Somatic mutations (MC3, non-silent, gene level) | `mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz` |
| Gene-level copy number (SNP6, log2 ratio) | `broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena.gz` |
| Tumour purity (ABSOLUTE) | `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` |
| Glioma molecular subtypes (revision analyses) | `TCGASubtype.20170308.tsv.gz` |

About 1 GB is downloaded on the first run; later runs reuse a local cache.

## How to run

### Google Colab (recommended)

1. Open `notebooks/01_DRG_pancancer_pipeline.ipynb` in Google Colab.
2. Select **Runtime → Run all**. Google Drive is mounted and all data, caches, tables, figures and logs are written to `MyDrive/DRG_PanCancer/`.
3. After the main pipeline has finished, open `notebooks/02_DRG_revision_analyses.ipynb` and run it in the same way.

Approximate run times on a standard Colab runtime: main pipeline about 60 min on the first run (most of it the 100-split LASSO validation); revision analyses 25–40 min.

### Local machine

```bash
git clone https://github.com/Sajjadeslamkhah/drg-pancancer.git
cd drg-pancancer
pip install -r requirements.txt

python src/drg_pancancer_pipeline.py --base ./DRG_PanCancer          # main analysis
python src/make_main_figures.py --results ./DRG_PanCancer/results    # Figures 1-6
python src/make_tables.py --results ./DRG_PanCancer/results           # Tables 1-2
python src/drg_v3_revision.py --base ./DRG_PanCancer                 # revision analyses
```

At least 12 GB of RAM is recommended. A quick smoke test on synthetic data (about 2 min, no download) is available with
`python src/drg_pancancer_pipeline.py --base ./demo --demo`.

## Outputs

| Output folder | Content |
|---|---|
| `results/tables/` | All result tables (CSV), numbered by analysis module (`00_` … `11_`) |
| `results/figures/` | Per-module figures (PNG 300 dpi and PDF) |
| `results/figures_main/` | Manuscript Figures 1–6 (PNG/PDF/TIFF, 600 dpi) |
| `results/REPORT.md` | Summary of the main results |
| `results_v3/` | Tables, Supplementary Figure S9 and report of the revision analyses |
| `logs/` | Run logs |

Correspondence between code and manuscript:

| Manuscript item | Produced by |
|---|---|
| Figure 1 (expression, co-expression) | modules 01 and 10 of the pipeline; `figure1()` |
| Figure 2 (survival) | modules 03 and 04; `figure2()` |
| Figure 3 (core vs cytoskeleton subsets) | module 11; `figure3()` |
| Figure 4 (metabolic programmes) | module 06; `figure4()` |
| Figure 5, Table 2 (multivariable Cox, LASSO-Cox validation) | modules 07 and 08; `figure5()` |
| Figure 6 (clinical subgroups, mutations, copy number) | modules 05 and 09; `figure6()` |
| Tables 1 and 2 | `src/make_tables.py` |
| Supplementary Figure S9 and revision analyses | `src/drg_v3_revision.py` |

## Reproducibility

- All random procedures use the seed `2026`.
- The pipeline was run with Python 3.13, pandas 2.2, NumPy 2.1, SciPy, statsmodels, matplotlib and seaborn (`requirements.txt`).
- The closed-form ssGSEA implementation (τ = 0.25) was verified against a brute-force running-sum computation during development, and the revision module checks its fast ranking engine against the pipeline's ssGSEA at start-up.
- Each completed run replaces the previous outputs; an interrupted run leaves earlier outputs unchanged.

## Citation

If you use this code, please cite the article above (see `CITATION.cff`). The citation will be updated with the DOI upon publication.

## License

MIT License (see `LICENSE`).

## Contact

Sajjad Eslamkhah — Biruni Advanced Technology Application and Research Center (B@MER), Biruni University, Istanbul, Türkiye — mtaghizad@biruni.edu.tr
