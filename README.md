
# RNA-seq Analysis Workflow

A reproducible RNA-seq analysis workflow implemented with **Snakemake**.  
This workflow processes raw RNA-seq data, performs quality control, aligns reads, quantifies gene expression, and identifies differentially expressed genes.

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

- Linux  
- Miniconda or Anaconda installed  
- 4 cores and at least 8 GB RAM (recommended)  

---

## Step 2: Create Conda Environment for Snakemake

- conda create -n snakemake7 python=3.10 snakemake
---

## Step 3: Activate the Environment
- conda activate snakemake7
---

## Step 4: Clone the Repository
- git clone https://github.com/<your-username>/RNA_seq.git
- cd RNA_seq
---

## Step 5: Configure Samples
### Edit config/samples.tsv to include your RNA-seq samples:
#### sample	condition

- Sample1  control
- Sample2	control
- Sample3  control
- Sample4	control	
- Sample5  treated	
- Sample6	treated	
- Sample7  treated	
- Sample8	treated
---

## Step 6: Configure Paths
### Edit config/config.yaml to set paths for your data, reference genome, and results:
- paths:
-   data_dir: "data"
-   reference_dir: "reference"
-   results_dir: "results"
-   adapters: "reference/adapters.fa"
---

## Step 7: Activate Logging and Dry Run
- nakemake --use-conda --cores 4 -n
---

## Step 8: Run Snakemake Workflow
- snakemake --use-conda --cores 4
---
## Contact
- Developed by Sasti Gopal Das
- Email: sastigopaldas05@gmail.com
---


