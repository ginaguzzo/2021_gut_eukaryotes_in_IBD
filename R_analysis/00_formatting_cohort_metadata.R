## Creating full metadata file for 500FG and 1000IBD ##

library(dplyr)
library(tidyr)

setwd("R_analysis/")


# Create 500FG metadata file ------------------------------------------------------------

## Metadata access from: https://hfgp.bbmri.nl/menu/500fgdata/dataexplorer?entity=FG500_public_data
## Click download, select 'Attribute labels', 'Entity labels' and download as CSV file.
## Save to R_analysis/metadata_files as 500FG_client_records.csv

## Import table
md1 <- read.csv("metadata_files/500FG_client_records.csv")

## Rename and sub-select columns
head(md1)

md1 <- md1 %>% 
  rename(ID = X500FG.identifier,
         sex = Wat.is.uw.geslacht.,
         age = Wat.is.uw.leeftijd.._._._.jaar,
         oral_contraceptive = Gebruikt.u.de.pil.als.anticonceptiemiddel.,
         smoking_status = Bent.u.een.roker.of.een.roker.geweest.) %>%
  select(ID, sex, age, smoking_status, oral_contraceptive)

## Translate data
unique(md1$smoking_status)

md1$smoking_status <- ifelse(md1$smoking_status == "nooit gerookt", "never",
                                ifelse(md1$smoking_status == "roker in het verleden", "past",
                                       ifelse(md1$smoking_status == "huidige roker", "current", NA)))



## Access sample identifiers on NCBI at: https://www.ebi.ac.uk/ena/browser/view/PRJNA319574?show=reads
## Under Column Selection, select 'run_accession' and 'experiment_title'.
## Download as TSV and save to R_analysis/metadata_files as 500FG_run_experiment.txt

## Import table
df_ids <- read.table(file = "metadata_files/500FG_run_experiment.txt", 
                     sep = "\t", header = TRUE)

## Remove unnecessary row data
df_ids$experiment_title <- gsub(x = df_ids$experiment_title,
                                pattern = "Illumina HiSeq 2000 sequencing; WGS of healthy adult gut metagenomes in the 500FG cohort: Subject  ", 
                                replacement = "") 

## Remove sample IDs that are not present in my file list
df_files <- read.table(file = "metadata_files/500FG_file_list.txt", 
                  sep = "\t", header = TRUE)

df_ids <- semi_join(df_ids, df_files, by = "run_accession")


## Rename columns
df_ids <- df_ids %>%
  rename(sampleid = run_accession,
         ID = experiment_title)


## Merge tables
md_500fg <- left_join(df_ids, md1, by = "ID")


## Add height and weight obtained from website administrator
md_bmi <- read.csv("metadata_files/500FG_height_weight.csv")

md_500fg <- inner_join(md_500fg, md_bmi, by = "ID", all.x=TRUE)


## Calculate BMI as separate column
md_500fg$bmi <- (md_500fg$weight/((md_500fg$height/100)^2))


## Save file
write.table(md_500fg, file = "metadata_files/500FG_md1.csv", sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)




# Create 1000IBD metadata file --------------------------------------------

## Import 1000IBD metadata and remove sample IDs that are not present in my file list
md2 <- read.csv("metadata_files/1000IBD_md1.csv")

df_files <- read.table(file = "metadata_files/1000IBD_file_list.txt", 
                       sep = "\t", header = TRUE)

md_1000ibd <- inner_join(df_files, md2, by = "sampleid")

## Save file
write.table(md_1000ibd, file = "metadata_files/1000IBD_md2.csv", sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)




# Merge metadata files ---------------------------------------------------

## Import tables
df_ctrl <- read.csv(file = "metadata_files/500FG_md1.csv")
df_ibd <- read.csv(file = "metadata_files/1000IBD_md2.csv")

## Format 500FG to match 1000IBD
df_ctrl$disease <- "control"

df_ctrl$sex <- ifelse(df_ctrl$sex == "male", "Male",
                             ifelse(df_ctrl$sex  == "female", "Female", ""))

## Add smoking_status to 1000IBD metadata
df_ibd$smoking_status <- ifelse(df_ibd$smoking_past == 0 & df_ibd$smoking_now == 0, "never",
                                ifelse(df_ibd$smoking_past == 1 & df_ibd$smoking_now == 0, "past",
                                       ifelse(df_ibd$smoking_now == 1, "current", NA)))
  

## Merge tables
df_total <- bind_rows(df_ibd, df_ctrl)
df_total = subset(df_total, select = -c(ID)) #Remove column with 500FG IDs


write.table(df_total, file = "metadata_files/cohort_metadata_combined.csv", 
            sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)
