## Merge raw FMT metadata files and simplify. ##

library(dplyr)
library(tidyr)
library(tibble)

## Import metadata
md1 <- read_tsv("metadata_files/raw_metadata_Sokoletal2020.txt")
md2 <- read_tsv ("metadata_files/raw_metadata_Kongetal2020.txt")

## Change column names to be merged
md1 <- md1 %>% dplyr::rename("sample_description" = "Description")
md2 <- md2 %>% dplyr::rename ("sampleid" = "run_accession", "sample_description" = "sample_title")

## Remove participant initials from sample_description in md1 to match md2
md1$sample_description <- gsub(x = md1$sample_description, pattern = "_[A-Z][A-Z]_", replacement = "_") 
md1$sample_description <- gsub(x = md1$sample_description, pattern = "_donor_", replacement = "_") 
md1$sample_description <- gsub(x = md1$sample_description, pattern = "_donor", replacement = "") 

## Simplify metadata
head(md1)
md1 <- dplyr::select(md1, -c(study_accession, sample_accession, secondary_sample_accession, fastq_file, `#SampleID`))
head(md1)

## As discovered in downstream analyses, participant 21 is missing info for two of their samples in the dataset
## Duplicate their info to the missing samples
which(md1$ID_sample == "1_21_PR")
md1 <- rbind(md1, md1[rep(105, 2),])
which(md1$ID_sample == "1_21_PR")

## Replace distinct values
which(colnames(md1) == "Timepoint")
md1[117, 9] = "D0" 
md1[118, 9] = "W-2"

which(colnames(md1) == "sample_description")
md1[117, 24] = "1_21_D0" 
md1[118, 24] = "1_21_W-2" 

## Remove duplicate samples from md2 that are from newer upload batch
md2 <- md2[!grepl("SRR14", md2$sampleid),]


## Merge 
md <- right_join(md1, md2, by = "sample_description")
md <- dplyr::relocate(md, "sampleid")

## Save
write.table(md, "metadata_files/fmt_metadata.csv", 
            sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)




