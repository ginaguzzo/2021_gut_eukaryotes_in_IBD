# 2021_gut_eukaryotes_in_IBD

## Overview
This repository includes the code and output for detecting eukaryotes in gut microbiomes of individuals with inflammatory bowel disease (IBD) and individuals without IBD. 

The analysis is split into two sub-analyses of three shotgun metagenomic datasets and their respective metadata:

**1.	Cohort meta-analysis**
   -	Individuals with IBD from the 1000IBD cohort ([Imhann et al. 2019](https://bmcgastroenterol.biomedcentral.com/articles/10.1186/s12876-018-0917-5)) available [here](https://ega-archive.org/datasets/EGAD00001004194). *NOTE: the dataset is available per request following contact with the Data Access Committee.*
   -	Healthy individuals from the 500FG cohort ([Schirmer et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5131922/)) available [here](https://www.ebi.ac.uk/ena/browser/view/PRJNA319574).
 
**2.	Faecal microbiota transplant (FMT) analysis**
   - FMT recipients and FMT donors from ([Kong et al. 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7725862/)) available [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA625520/).

## Dependencies

#### Conda packages
-	FastQC
-	MultiQC 
-	GNU Parallel
-	Fastp (v0.20.0)
-	KneadData (v0.7.2)
-	Pigz
-	Seqtk
-	EukDetect (v1.2)
-	RiboTagger (v0.8.0)

#### R packages
- R version >=4.0.2
- Phyloseq (v1.38.0)

## Workflow
The workflow is run in the following order:

[**1. Bioinformatic_pipeline**](https://github.com/ginaguzzo/2021_gut_eukaryotes_in_IBD/tree/main/Bioinformatic_pipeline)
   - Shell commands for quality control, detecting eukaryotes in samples, and formatting output for R.

**2. R_analysis**
   - R code for results exploration and generating figures and tables.

## Setting up
1.	Clone this repository onto your local machine.
```
git clone https://github.com/ginaguzzo/2021_gut_eukaryotes_in_IBD.git
cd 2021_gut_eukaryotes_in_IBD
```

2.	Download samples and move them to their respective directories.
```
mv EGA*gz Bioinformatic_pipeline/1000ibd_data/raw_files
mv SRR5*gz Bioinformatic_pipeline/500fg_data/raw_files
mv SRR1*gz Bioinformatic_pipeline/fmt_data/raw_files
```

3. Metadata files for the 1000IBD and 500FG cohorts should be moved to: 
```
R_analysis/cohort_data_analysis/metadata_files
```

4. Metadata for the FMT data should be moved to: 
```
R_analysis/fmt_data_analysis/metadata_files
```


## Running the pipeline
1. Change to the pipeline directory.
```
cd Bioinformatic_pipeline
```

2. Make a conda environment called `euks_in_ibd` with py2.7.
```
conda create -n euks_in_ibd python=2.7
```

3. Install the [required conda packages](https://github.com/ginaguzzo/2021_gut_eukaryotes_in_IBD#conda-packages) EXCEPT for eukdetect.

     *NOTE:* Clone the Kneaddata and RiboTagger repos to `Bioinformatic_pipeline/tools/` before installing
 
4.	Eukdetect should be installed in its own environment according to the instructions.

     *NOTE:* Clone the EukDetect repo to `Bioinformatic_pipeline/tools/`

5.	Edit the paths in the .yml files in `eukdetect_config_files` to include the absolute path from your home directory. 

6.	Copy config files to the Eukdetect installation directory.
```
cp eukdetect_config_files/*yml tools/EukDetect-master
```

7.	The pipeline scripts should be run in the order they are numbered, from 1-8.


## Running the R analysis
1. Install the [required R packages](https://github.com/ginaguzzo/2021_gut_eukaryotes_in_IBD#r-packages) in R studio. 

2. The scripts should be run in the order they are numbered.



*NOTE:* Due to data access constraints, all participant-associated metadata files have been removed. 

For reproducibility, example metadata files can be found in 
```
R_analysis/cohort_data_analysis/metadata_files
```
and 
```
R_analysis/fmt_data_analysis/metadata_files
```

These contain headers and mock samples of how the metadata was formatted.

