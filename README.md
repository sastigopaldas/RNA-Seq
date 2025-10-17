
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

## Installation

### Step 1: Install Conda (if not installed)

```bash
# Download Miniconda (Linux example)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run installer
bash Miniconda3-latest-Linux-x86_64.sh

# Reload shell
source ~/.bashrc

# Verify installation
conda --version











A reproducible RNA-seq analysis workflow implemented with Snakemake.


**Step 1: If Conda is not installed on your system, install it first.**

**Step 2: Create a Conda environment for Snakemake:**

conda create -n snakemake7 python=3.10 snakemake
