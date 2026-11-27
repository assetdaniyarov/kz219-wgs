# Changelog

## v1.0.0 (2026-11-27)

Release archived on Zenodo for the submission to *Communications Biology*.

* Upstream scripts added: population genetics (`pipeline/popgen`, PLINK / ADMIXTURE / F_ST /
  ROH / AADR projection / tribe analyses, with all parameters and seeds in `config.sh`),
  mitochondrial haplogroups (`pipeline/04_mtdna_haplogroups.py`) and Y-chromosome haplogroups
  (`pipeline/05_chrY_ylt.sh`), alongside the existing DRAGEN and ANNOVAR steps.
* Reference-panel metadata used by the popgen stages (`pipeline/popgen/meta`).
* README: correspondence between notebook outputs and the item numbering of the article;
  random seeds; versions of every upstream tool.
* `CITATION.cff` with the full author list and ORCIDs; `.zenodo.json` metadata for the
  Zenodo release; `.gitignore` and `.gitattributes`.
