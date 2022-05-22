## RiboTagger results analysis

library(dplyr)
library(tidyr)
library(tibble)

setwd("R_analysis/")


# Create Table S1 ---------------------------------------------------------


## Import tables
df.ibd <- read.table(file = "results/1000ibd_ribotagger_all_v_region_hits.txt", 
                  sep = "\t", header = TRUE)
df.ctl <- read.table(file = "results/500fg_ribotagger_all_v_region_hits.txt", 
                  sep = "\t", header = TRUE)

## Create table of results
## Add column for disease state
df.ibd$disease <- 'IBD'
ncol(df.ibd)
df.ibd <- df.ibd[c(1, 19, 2:18)]

df.ctl$disease <- 'control'
ncol(df.ctl)
df.ctl <- df.ctl[c(1, 19, 2:18)]

## Combine results and save file
df.total <- bind_rows(df.ibd, df.ctl)

## Remove rows with plant IDs
df.total <- filter(df.total, c != "Embryophyta")

## Rename columns
tables1 <- df.total %>%
  rename("Sample ID" = sampleid, 
         "Disease status" = disease,
         "Confidence" = confidence,
         "Kingdom" = k,
         "Phylum" = p,
         "Class" = c,
         "Order" = o,
         "Family" = f,
         "Genus" = g,
         "18S V region" = v_region)

## Sub-select table
tables1 <- tables1 %>% 
  select(-c(tag, use, taxon_level, taxon_data, 
         long, long_total, long_this, support, s))

## Round decimal places in Confidence column
tables1$Confidence <- format(round(tables1$Confidence, 2), nsmall = 2)


write.table(tables1, file = "figures_and_tables/Table_S1_ribotagger_results.csv", 
            sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)



# Explore the data --------------------------------------------------------


## IBD samples ##
df.ibd1 <- unique(df.ibd[c("sampleid", "f",  "g")]) #Make new df of samples with unique results at family or genus level
df.ibd1 <- df.ibd1[!(df.ibd1$f=="" & df.ibd1$g==""),] #Remove samples with no family or genus listed
df.ibd1 <- df.ibd1[!(df.ibd1$f=="Liliopsida" | df.ibd1$g=="Sinapis"),] #Remove samples with plant taxa identified
#df.ibd1 <- distinct(df.ibd1, sampleid, f, .keep_all= TRUE) #Remove duplicate family-level hits in same samples

colSums(table(df.ibd1$f, df.ibd1$g, df.ibd1$f)) #Summarise taxa
length(unique(df.ibd1$sampleid)) #Total individuals with microbial eukaryotes
colSums(table(df.ibd1$sampleid, df.ibd1$g)) 

sum(duplicated(df.ibd1$sampleid)) #Number of duplicates
dups.ibd <-df.ibd1[duplicated(df.ibd1$sampleid)|duplicated(df.ibd1$sampleid, fromLast=TRUE),] #Find duplicate IDs
dups.ibd
table(dups.ibd$sampleid) #Number of duplicates per sample




## Control samples ##
df.ctl1 <- unique(df.ctl[c("sampleid", "f",  "g")]) #Make new df of samples with unique results at family or genus level
df.ctl1 <- df.ctl1[!(df.ctl1$f=="" & df.ctl1$g==""),] #Remove samples with no family or genus listed
df.ctl1 <- df.ctl1[!(df.ctl1$f=="Liliopsida" | df.ctl1$g=="Sinapis"),] #Remove samples with plant taxa identified
#df.ctl1 <- distinct(df.ctl1, sampleid, f, .keep_all= TRUE) #Remove duplicate family-level hits in same samples

colSums(table(df.ctl1$f, df.ctl1$g, df.ctl1$f)) #Summarise taxa
length(unique(df.ctl1$sampleid)) #Total individuals with microbial eukaryotes

sum(duplicated(df.ctl1$sampleid)) #Number of duplicates
dups.ctl <-df.ctl1[duplicated(df.ctl1$sampleid)|duplicated(df.ctl1$sampleid, fromLast=TRUE),] #Find duplicate IDs
dups.ctl
table(dups.ctl$sampleid) #Number of duplicates per sample

