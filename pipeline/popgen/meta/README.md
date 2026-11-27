# Reference-panel metadata

Hand-curated files read by the popgen stages from `${META_DIR}` (`${PROJECT}/results/meta`).
They describe the public reference panels only (HGDP, 1000 Genomes, SGDP, AADR v62) and
contain no participant data. Copy them into `${META_DIR}` before running `check_inputs.sh`;
`superpop_map.tsv` sits next to the scripts and is copied there by `setup_meta.sh`.

| File | Content |
|---|---|
| `aadr_keep_groups.txt` | AADR groups retained for the projection (Supplementary Fig. 9) |
| `aadr_selection.tsv` | AADR individuals selected from those groups |
| `duplicate_labels.tsv` | individuals present in more than one reference panel |
| `keep_nations_resolution.tsv` | resolution of the `KEEP_NATIONS` tokens in `config.sh` |
| `nation_aliases.tsv` | harmonisation of population labels across panels |
| `panel_exclude_iids.txt` | reference individuals excluded before analysis |
| `pseudohaploid_iids.txt` | pseudohaploid individuals, projected rather than used for the axes |

The cohort files `kz_219_metadata.tsv` and `kz_219_tribes.tsv`, also read from `${META_DIR}`,
are participant metadata and are not distributed (see `data/README.md`). The list of high-LD
regions is written by `02_merge.sh`.
