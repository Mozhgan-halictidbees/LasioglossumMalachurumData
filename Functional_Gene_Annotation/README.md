# eggNOG-mapper Functional Annotation Pipeline 

## Overview

This workflow provides the code for functional annotation of *Lasioglossum malachurum* protein sequences using **eggNOG-mapper v2.1.12**.

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

## What the pipeline does
Creates working directories
Downloads eggNOG database (Hymenoptera, Tax ID: 7399)
Builds DIAMOND database
Runs eggNOG-mapper annotation
Saves results in output folder

---

##  Reproducibility

To ensure reproducibility:

eggNOG-mapper v2.1.12
DIAMOND v2.0.11
eggNOG database v5.0.2
