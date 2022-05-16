# Analysis of microbial eukaryotes in IBD gut microbiomes  

## Overview
This repository includes the code and output for detecting eukaryotes in gut microbiomes of individuals with (and without) inflammatory bowel disease (IBD). 

The analysis is split into two sub-analyses of three shotgun metagenomic datasets:

**1.	Cohort meta-analysis**
   -	Individuals with IBD from the 1000IBD cohort ([Imhann et al. 2019](https://bmcgastroenterol.biomedcentral.com/articles/10.1186/s12876-018-0917-5)).
   -	Healthy individuals from the 500FG cohort ([Schirmer et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5131922/)).
 
**2.	Faecal microbiota transplant (FMT) analysis**
   - FMT recipients and FMT donors from [Kong et al. 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7725862/).

## Downloading files

#### Metagenomic samples
- 1000IBD samples can be found [here](https://ega-archive.org/datasets/EGAD00001004194), and are accessible subject to approval by the Data Access Committee.
- 500FG samples can be downloaded [here](https://www.ebi.ac.uk/ena/browser/view/PRJNA319574).
- FMT study samples can be downloaded [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA625520/).

#### Sample metadata
Due to ethical requirements, all participant-associated metadata files have been removed from this repository. For reproducibility, example metadata files can be found in 
```
R_analysis/metadata_files
```


These contain headers and mock samples of how the metadata were formatted.


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
- phyloseq (v1.38.0)
- gridExtra (v2.3)
- patchwork (v1.1.1)
- cowplot (v1.1.1)
- pscl (v1.5.5)
- sandwich (v3.0-1)
- lmtest (v0.9-39)
- ggh4x (v0.2.1)


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

3. Metadata files for all three cohorts should be moved to: 
```
R_analysis/metadata_files
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

3. Create a directory to install conda packages.
```
mkdir tools
```

4. Install the [required conda packages](https://github.com/ginaguzzo/2021_gut_eukaryotes_in_IBD#conda-packages) EXCEPT for EukDetect.

     *NOTE:* Clone the KneadData and RiboTagger repos to `Bioinformatic_pipeline/tools/` before installing
 
5.	Eukdetect should be installed in its own environment according to the instructions.

     *NOTE:* Clone the EukDetect repo to `Bioinformatic_pipeline/tools/`

6.	Edit the paths in the .yml files in `eukdetect_config_files` to include the absolute path from your home directory. 

7.	Copy config files to the Eukdetect installation directory.
```
cp eukdetect_config_files/*yml tools/EukDetect-master
```

8.	The pipeline scripts should be run in the order they are numbered, from 1-8.


## Running the R analysis
1. Install the [required R packages](https://github.com/ginaguzzo/2021_gut_eukaryotes_in_IBD#r-packages) in R studio. 

2. The scripts are located in `R_analysis` and should be run in the order they are numbered.



