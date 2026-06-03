# eggNOG-mapper Functional Annotation Pipeline

## Overview

This workflow provides the code for functional annotation of *Lasioglossum malachurum* protein sequences using eggNOG-mapper v2.1.12.

The pipeline performs orthology-based annotation and assigns:
- Gene Ontology (GO) terms
- KEGG Orthology (KO) assignments
- COG functional categories
- Protein functional descriptions

---

## Input Data

Protein FASTA file:

LMAR_OGS_v2.1.1_pep.fasta

Source:
- Halictid Genome Browser
- https://beenomes.princeton.edu/
- Downloaded: April 2025

---

## Software Requirements

| Tool | Version |
|------|--------|
| eggNOG-mapper | v2.1.12 |
| DIAMOND | v2.0.11 |
| eggNOG database | v5.0.2 |
| Python | >=3.8 |

---

## Pipeline Overview

The workflow:
- Creates working directories
- Downloads eggNOG database (Hymenoptera, Tax ID: 7399)
- Builds DIAMOND database
- Runs eggNOG-mapper annotation
- Saves results in output folder

---

## Key Parameters

- Alignment method: DIAMOND (sensitive mode)
- E-value threshold: 0.001

---

## Reproducibility

To ensure reproducibility, DIAMOND version should be compatible with eggNOG-mapper:

- eggNOG-mapper v2.1.12
- DIAMOND v2.0.11
- eggNOG database v5.0.2

---

## Citations

[1] eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at metagenomic scale.
Cantalapiedra CP et al., 2021. Molecular Biology and Evolution
https://doi.org/10.1093/molbev/msab293

[2] eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource.
Huerta-Cepas J et al., 2019. Nucleic Acids Research
https://doi.org/10.1093/nar/gky1085

[3] Sensitive protein alignments at tree-of-life scale using DIAMOND.
Buchfink B et al., 2021. Nature Methods
https://doi.org/10.1038/s41592-021-01101-x



