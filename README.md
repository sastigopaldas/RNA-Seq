
# RNA-seq Analysis Workflow

A reproducible RNA-seq analysis workflow implemented with **Snakemake**.  
This workflow processes raw RNA-seq data, performs quality control, aligns reads, quantifies gene expression, and identifies differentially expressed genes.

---

## Table of Contents

- [Features](#features)  
- [Requirements](#requirements)  
- [Installation](#installation)  
- [Workflow Usage](#workflow-usage)  
- [Configuration](#configuration)  
- [License](#license)  

---

## Features

- Fully automated RNA-seq workflow using **Snakemake**  
- Read quality control with **FastQC**  
- Alignment using **STAR**  
- Gene quantification with **featureCounts**  
- Differential expression analysis with **DESeq2**  
- Reproducible environment via **Conda**  

---

## Requirements

- Linux or MacOS system  
- Miniconda or Anaconda installed  
- 4 cores and at least 8 GB RAM (recommended)  

---

## Step 2: Create Conda Environment for Snakemake

### conda create -n snakemake7 python=3.10 snakemake


