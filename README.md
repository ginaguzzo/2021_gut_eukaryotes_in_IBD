# 2021_gut_eukaryotes_in_IBD

## Overview
This repository includes the code and output for detecting eukaryotes in gut microbiomes of individuals with inflammatory bowel disease (IBD) and individuals without IBD. 

The analysis is split into two sub-analyses of three shotgun metagenomic datasets and their respective metadata:

**1.	Cohort meta-analysis**
   -	Individuals with IBD from the 1000IBD cohort (Imhann et al. 2019) available here. *NOTE: the dataset is available per request following contact with the Data Access Committee.*
   -	Healthy individuals from the 500FG cohort (Schirmer et al. 2019) available here.
 
**2.	Faecal microbiota transplant (FMT) analysis**
   - FMT recipients and FMT donors from (Kong et al. 2020) available here.

## Dependencies

#### Conda packages
-	FastQC
-	MultiQC 
-	GNU Parallel
-	Fastp (v0.20.0)
-	Kneaddata (v0.7.2)
-	Pigz
-	Seqtk
-	Eukdetect (v1.2)
-	RiboTagger (v0.8.0)

#### R packages
- Phyloseq (v1.38.0)

## Workflow
The workflow is run in the following order:

**1. Bioinformatic_pipeline**
   - Shell commands for quality control, detecting eukaryotes in samples, and formatting output for R.

**2. R_analysis**
   - R code for results exploration and generating figures and tables.

See the respective README files for details on how to run each of these analyses.
