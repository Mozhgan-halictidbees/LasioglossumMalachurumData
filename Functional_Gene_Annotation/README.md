# eggNOG-mapper Functional Annotation Pipeline (One-Click)

## Overview

This repository provides the code for functional annotation of *Lasioglossum malachurum* protein sequences using **eggNOG-mapper v2.1.12**.

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
| Python | ≥3.8 |

---

## One-Click Pipeline

Run the full workflow:

```bash
bash eggnog_pipeline.sh
What the pipeline does
Creates working directories
Downloads eggNOG database (Hymenoptera, Tax ID: 7399)
Builds DIAMOND database
Runs eggNOG-mapper annotation
Saves results in output folder
Output Files

All results are saved in:

output/
File	Description
.emapper.annotations	Functional annotation table
.emapper.hits	DIAMOND alignment results
.emapper.seed_orthologs	Orthology mapping results
Biological Method Summary

Orthology-based functional annotation of Lasioglossum malachurum was performed using eggNOG-mapper v2.1.12.

Protein sequences (LMAR_OGS_v2.1.1_pep.fasta), downloaded in April 2025, were annotated using eggNOG database v5.0.2.

Sequence similarity searches were performed using DIAMOND v2.0.11 in sensitive mode with an e-value threshold of 0.001.

Functional annotations include:

Gene Ontology (GO)
KEGG Orthology (KO)
COG functional categories
Gene descriptions
Reproducibility

To ensure reproducibility:

eggNOG-mapper v2.1.12
DIAMOND v2.0.11
eggNOG database v5.0.2
