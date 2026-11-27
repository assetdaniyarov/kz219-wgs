# Input data

This folder is excluded from version control (see `.gitignore`); only this file is tracked.
Individual-level genotypes and participant metadata (sex, self-reported tribe and zhuz) are
not distributed with the code. Raw reads are deposited in NCBI SRA under BioProject
[PRJNA1224597](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1224597); access to the
participant metadata is described in the Data availability statement of the article.

Expected layout, with the environment variable that overrides each location:

```
data/
  dragen/                                   KZ_DRAGEN_DIR   notebook 01
    <sample>/<sample>.mapping_metrics.csv
    <sample>/<sample>.wgs_coverage_metrics.csv
    <sample>/<sample>.ploidy_estimation_metrics.csv
    <sample>/<sample>.vc_metrics.csv
    <sample>/<sample>.fragment_length_hist.csv
    <sample>/<sample>.fastqc_metrics.csv
    <sample>/<sample>.vc_hethom_ratio_metrics.csv
    <sample>/<sample>.roh_metrics.csv

  popgen/                                   KZ_PROJECT      notebook 02
    results/
      groupB_kz_hgdp/                       PCA, ADMIXTURE, F_ST, ROH, AADR projection
      groupA_kz_1000g_sgdp/                 F_ST with 1000 Genomes and SGDP
      kz_tribes/                            within-cohort PCA, PC association, kinship
      kz_only/                              within-cohort ROH
      meta/                                 kz_219_tribes.tsv, panel_qc.tsv

  kz_219_hg19_dragen.annotated.header.txt   KZ_CALLSET      notebook 03
  kz_219_tribes.tsv                         KZ_SEX_FILE     notebook 03 (IID, Gender)

  chrM/                                     KZ_CHRM_DIR     notebook 05
    kz_219_mtdna_haplogroups.tsv            haplogrep3 --extend-report
    kz_219_chrM_metadata_final.xlsx         maternal tribe and zhuz

  chrY/                                     KZ_CHRY_DIR     notebook 06
    kz_219_ylt.hapresult.hg                 Y-LineageTracker v1.3.0, classify -b 37
    kz_219_ylt.lineageresult.txt
    metada_with_yhaplogroups.csv            paternal tribe and zhuz
```

Notebook 04 reads `work/03_af/kz_219_hg19_AF_from_GT_known.txt` (`KZ_KNOWN_FILE`), which
notebook 03 writes; run 03 first.
