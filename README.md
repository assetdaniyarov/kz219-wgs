# KZ-219: whole-genome sequencing of 219 Kazakh individuals

![cohort](https://img.shields.io/badge/cohort-KZ--219-1b9e77)
![platform](https://img.shields.io/badge/platform-DNBSEQ--G400-7570b3)
![pipeline](https://img.shields.io/badge/DRAGEN-v3.10.16-d95f02)
![reference](https://img.shields.io/badge/reference-GRCh37%2Fhg19-e7298a)
![annotation](https://img.shields.io/badge/ANNOVAR-v2025Mar02-666666)
![bioproject](https://img.shields.io/badge/BioProject-PRJNA1224597-333333)
![license](https://img.shields.io/badge/license-MIT-blue)
<!-- after the Zenodo release, replace XXXXXXX with the record number:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

Analysis code for the article

> Daniyarov A. *et al.* High-coverage whole genomes of 219 Kazakhs fill a Central Asian gap
> in human reference panels. *Communications Biology* (2026).
> <!-- doi: 10.1038/s42003-026-XXXXX-X -->

<img width="6248" height="3680" alt="figure1" src="https://github.com/user-attachments/assets/24bf9a7d-b62b-44cc-8c7e-33a2476f38f6" />

The upstream scripts in [`pipeline/`](pipeline) record the commands, software versions,
parameters and random seeds that produced the notebook inputs (DRAGEN, ANNOVAR, PLINK,
ADMIXTURE, samtools/BWA/mutserve/HaploGrep 3, Y-LineageTracker). The six Jupyter notebooks
take those outputs and produce the figures, tables and Supplementary Data of the article.

| | |
|---|---|
| **Cohort** | 219 healthy Kazakh adults (145 XX, 74 XY) |
| **Sequencing** | DNBSEQ-G400, MGIEasy PCR-Free, 100 bp paired-end |
| **Alignment and variant calling** | Illumina DRAGEN Bio-IT Platform v3.10.16; joint genotyping DRAGEN v4.3.6 |
| **Reference** | GRCh37 / hg19 |
| **Annotation** | ANNOVAR v2025Mar02 |
| **Raw data** | NCBI SRA, BioProject [PRJNA1224597](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1224597) |

## Notebooks and manuscript items

The notebooks were executed before the final renumbering of the manuscript, so output file
names and the item numbers printed inside the notebooks (set in the `MAIN_FIG` and `SUPP_NUM`
dictionaries, or in the file names in notebooks 01 and 02) follow the earlier numbering. The
table gives the numbering of the published article.

| Notebook | Topic | Main figures and tables | Supplementary figures | Supplementary tables | Supplementary Data |
|---|---|---|---|---|---|
| [`01_wgs_qc_metrics`](01_wgs_qc_metrics.ipynb) | sequencing, alignment and callset QC | | 2 | 1 | 1 |
| [`02_population_genetics`](02_population_genetics.ipynb) | PCA, ADMIXTURE, F<sub>ST</sub>, ROH, AADR projection | Fig. 2 | 3, 4, 5, 9 | 2, 3, 18, 19, 20 | 2 |
| [`03_allele_frequencies`](03_allele_frequencies.ipynb) | cohort allele frequencies, variants absent from reference panels | Table 4 | 7 | 9, 10, 21-28 | 4 |
| [`04_biomedical_pgx`](04_biomedical_pgx.ipynb) | pharmacogenomic and clinically relevant variants | Figs. 3, 4; Tables 1-3 | 6 | 4-8, 29 | 3 |
| [`05_chrM_analysis`](05_chrM_analysis.ipynb) | mitochondrial haplogroups | Fig. 5 | 8 | 11-13 | 5 |
| [`06_chrY_analysis`](06_chrY_analysis.ipynb) | Y-chromosome haplogroups | Fig. 6 | | 14-17 | 6 |

Fig. 1 and Supplementary Fig. 1 (study design and workflow diagrams) are not produced by
these notebooks. Main-text Tables 1-3 are drawn from Supplementary Tables 4 and 5 and
Supplementary Data 3.

<details>
<summary>Item-by-item correspondence between notebook outputs and the article</summary>

| Notebook | Label in notebook output | Item in the article |
|---|---|---|
| 01 | Supplementary Fig. 2, Table 1, Data 1 | Supplementary Fig. 2, Table 1, Data 1 |
| 02 | `Fig2_PCA_ADMIXTURE` | Fig. 2 |
| 02 | `Fig3_ROH` | Supplementary Fig. 5 |
| 02 | `Supplementary_Fig_3_tribe_PCA` | Supplementary Fig. 3 |
| 02 | `Supplementary_Fig_4_AADR_projection` | Supplementary Fig. 9 |
| 02 | `Supplementary_Fig_5_Fst_heatmaps` | Supplementary Fig. 4 |
| 02 | Supplementary Tables 2 / 3 / 4 / 5 / 6 | Supplementary Tables 18 / 2 / 3 / 19 / 20 |
| 02 | Supplementary Data 2 | Supplementary Data 2 |
| 03 | Supplementary Fig. 6 | Supplementary Fig. 7 |
| 03 | Supplementary Table 7 (callset composition) | Table 4 |
| 03 | Supplementary Tables 8 / 9 / 10 / 11 / 12 | Supplementary Tables 21 / 9 / 22 / 23 / 24 |
| 03 | Supplementary Tables 13 / 14 / 15 / 16 / 17 | Supplementary Tables 10 / 25 / 26 / 27 / 28 |
| 03 | Supplementary Data 3 | Supplementary Data 4 |
| 04 | Fig. 4 / Fig. 5 | Fig. 3 / Fig. 4 |
| 04 | Supplementary Fig. 7 | Supplementary Fig. 6 |
| 04 | Supplementary Tables 18 / 19 / 20 / 21 / 22 / 23 | Supplementary Tables 4 / 5 / 6 / 8 / 29 / 7 |
| 04 | Supplementary Data 4 | Supplementary Data 3 |
| 05 | Fig. 6; Supplementary Fig. 8 | Fig. 5; Supplementary Fig. 8 |
| 05 | Supplementary Tables 24 / 25 / 26 | Supplementary Tables 11 / 12 / 13 |
| 05 | Supplementary Data 5 | Supplementary Data 5 |
| 06 | Fig. 7 | Fig. 6 |
| 06 | Supplementary Tables 27 / 28 / 29 / 30 | Supplementary Tables 14 / 15 / 16 / 17 |
| 06 | Supplementary Data 6 | Supplementary Data 6 |

</details>

## Repository layout

```
.
├── 01_wgs_qc_metrics.ipynb
├── 02_population_genetics.ipynb
├── 03_allele_frequencies.ipynb
├── 04_biomedical_pgx.ipynb
├── 05_chrM_analysis.ipynb
├── 06_chrY_analysis.ipynb
├── pipeline/
│   ├── 01_dragen_per_sample.sh        DRAGEN v3.10.16 alignment and gVCF calling
│   ├── 02_dragen_joint_genotyping.sh  DRAGEN v4.3.6 iterative gVCF genotyper
│   ├── 03_annovar_nf.sh               ANNOVAR v2025Mar02 via annovar-nf
│   ├── 04_mtdna_haplogroups.py        chrM extraction, rCRS remapping, mutserve, HaploGrep 3
│   ├── 05_chrY_ylt.sh                 Y-LineageTracker v1.3.0
│   └── popgen/                        PLINK / ADMIXTURE / F_ST / ROH / AADR / tribe analyses
├── figstyle.py            plotting helpers for notebook 02
├── kz219_style.py         plotting helpers for notebook 03
├── run_all.sh             runs the notebooks in order
├── requirements.txt
├── environment.yml
├── data/                  inputs (not tracked, see data/README.md)
├── work/                  large intermediates from notebook 03 (not tracked)
└── results/               figures and tables (not tracked)
```

## Upstream processing

Paths are set through environment variables; nothing in the scripts depends on the host on
which the analysis was run.

| Step | Script | Tool and version | Output |
|---|---|---|---|
| 1 | `pipeline/01_dragen_per_sample.sh` | DRAGEN v3.10.16, hg19 hash table v8 | BAM, gVCF, CNV, SV, ROH and QC metrics per sample (input of notebook 01) |
| 2 | `pipeline/02_dragen_joint_genotyping.sh` | DRAGEN v4.3.6 iterative gVCF genotyper | multi-sample VCF of 219 genomes |
| 3 | `pipeline/03_annovar_nf.sh` | [annovar-nf](https://github.com/assetdaniyarov/annovar-nf), ANNOVAR v2025Mar02 | annotated callset (input of notebook 03) |
| 4 | `pipeline/04_mtdna_haplogroups.py` | samtools v1.21, BWA-MEM v0.7.12-r1039, mutserve v2.0.1, HaploGrep 3 v3.2.2 (phylotree-fu-rcrs 1.2) | `kz_219_mtdna_haplogroups.tsv` (input of notebook 05) |
| 5 | `pipeline/05_chrY_ylt.sh` | Y-LineageTracker v1.3.0, `classify -b 37` | `kz_219_ylt.*` (input of notebook 06) |
| 6 | `pipeline/popgen/run_all.sh` | PLINK v1.9.0-b.8, PLINK v2.00a5.12, ADMIXTURE v1.3.0 | `data/popgen/results` (input of notebook 02) |

### Population genetics (`pipeline/popgen`)

All parameters are defined once, in [`config.sh`](pipeline/popgen/config.sh), which every
stage sources; `00_env.sh` writes the versions of the binaries and the value of every
parameter to the log before the run.

```bash
cd pipeline/popgen
PROJECT=/path/to/popgen ./run_all.sh          # all stages
PROJECT=/path/to/popgen ./run_all.sh 04 05    # selected stages
```

| Stage | Content |
|---|---|
| `01_prepare.sh` | import and harmonisation of the cohort and reference panels (HGDP, 1000 Genomes, SGDP, AADR v62), genotype and sample QC |
| `02_merge.sh` | panel merging with explicit strand-flip and allele-mismatch classification |
| `03_relatedness.sh` | KING-robust kinship, `--king-cutoff 0.0884`, LD pruning (`--indep-pairwise 200kb 1 0.2`) |
| `04_pca.sh` | PCA, eigenvalue-based variance explained |
| `05_admixture.sh` | ADMIXTURE K = 2-10, ten seeds per K, capped training panel and projection (`-P`), uncapped joint fit at the selected K; `admixture_align.py` aligns replicates into modes |
| `06_fst.sh` | Hudson F<sub>ST</sub> (`--fst method=hudson`), block jackknife over 1,000-variant blocks |
| `07_roh.sh` | runs of homozygosity, PLINK 1.9 `--homozyg` |
| `08_project_aadr.sh` | projection of AADR individuals onto the reference axes |
| `09_kz_tribes.sh` | within-cohort PCA, PC association, kinship and F<sub>ST</sub> by tribe and juz |
| `kz233_relatedness.sh` | the 233 → 219 relatedness step of the cohort derivation |
| `popgen_tables.py` | summary tables read by notebook 02 |

### Random seeds

| Analysis | Seed |
|---|---|
| ADMIXTURE, capped training panel | `SEED_BASE + 100*K + rep`, `SEED_BASE = 20260101`, rep = 1-10 (`config.sh`, `05_admixture.sh`) |
| ADMIXTURE, uncapped joint fit | `SEED_BASE + 50000 + 100*K + rep` |
| ADMIXTURE, within-cohort | `SEED_BASE + 90000 + 100*K + rep` (`09_kz_tribes.sh`) |
| Notebooks (reservoir sampling of shared sites, bootstrap and permutation tests) | `numpy.random.default_rng(0)` |

The seed of every ADMIXTURE run is also written to `cv_raw.tsv` alongside its CV error and
log-likelihood.

## Installation

```bash
git clone https://github.com/assetdaniyarov/kz219-wgs.git
cd kz219-wgs

# pip
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# or conda
conda env create -f environment.yml && conda activate kz219
```

Figures use Arial; on Linux, Liberation Sans (metric-compatible) is used if Arial is not
installed.

## Running

Place the inputs as described in [`data/README.md`](data/README.md), then from the
repository root:

```bash
bash run_all.sh            # all six notebooks
bash run_all.sh 05 06      # selected notebooks
```

or open the notebooks in Jupyter and run them in order. Notebook 04 depends on the output of
notebook 03; the others are independent. All paths are relative to the repository root and
can be overridden with environment variables:

| Variable | Default | Notebook |
|---|---|---|
| `KZ_DRAGEN_DIR` | `data/dragen` | 01 |
| `KZ_PROJECT` | `data/popgen` | 02 |
| `KZ_CALLSET` | `data/kz_219_hg19_dragen.annotated.header.txt` | 03 |
| `KZ_SEX_FILE` | `data/kz_219_tribes.tsv` | 03 |
| `KZ_WORK_DIR` | `work/03_af` | 03 |
| `KZ_KNOWN_FILE` | `work/03_af/kz_219_hg19_AF_from_GT_known.txt` | 04 |
| `KZ_CHRM_DIR` | `data/chrM` | 05 |
| `KZ_CHRY_DIR` | `data/chrY` | 06 |
| `KZ_OUT_DIR` (`KZ_OUT` in 02) | `results/<notebook>` | all |

Notebook 03 makes one pass over the full joint callset (about 30 million rows), which takes
several hours and writes several GB to `work/03_af`; the scan is cached, so later runs finish
in seconds. Set `DRY_RUN_ROWS` in the notebook for a quick test.

## Outputs

Each notebook writes to its own folder in `results/`:

* figures as TIFF (CMYK, LZW) and PNG (RGB), 1200 dpi, Arial 8 pt, panel letters 10 pt;
* `Supplementary_Tables_0X.docx`, the supplementary tables in Arial without shading;
* `Supplementary_Data_X_*.xlsx`, the supplementary data, including source data for every
  plotted panel.

## Data availability

Raw sequencing reads are available from NCBI SRA under BioProject PRJNA1224597. Aggregate
allele frequencies and all summary tables are provided as Supplementary Data with the
article. Individual-level genotypes and participant metadata are not part of this repository;
see the Data availability statement of the article for access.

## Citation

See [`CITATION.cff`](CITATION.cff). Please cite the article and, for the code, the archived
release on Zenodo (DOI on the badge above).

## License

Code is released under the [MIT License](LICENSE).

## Contact

Asset Daniyarov ([asset.daniyarov@nu.edu.kz](mailto:asset.daniyarov@nu.edu.kz))  
Laboratory of Bioinformatics and Systems Biology, Center for Life Sciences,  
National Laboratory Astana, Nazarbayev University,  
53 Kabanbay Batyr Avenue, Astana 010000, Kazakhstan

Ulykbek Kairov ([ulykbek.kairov@nu.edu.kz](mailto:ulykbek.kairov@nu.edu.kz))  
Laboratory of Bioinformatics and Systems Biology, Center for Life Sciences,  
National Laboratory Astana, Nazarbayev University,  
53 Kabanbay Batyr Avenue, Astana 010000, Kazakhstan
