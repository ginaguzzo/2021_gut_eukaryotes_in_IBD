## EukDetect results analysis on unrarefied samples

library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(readr)

setwd("R_analysis/")

# Data wrangling ----------------------------------------------------------

## Import tables and combine
df.1000ibd <- read_tsv(file = "results/1000ibd_eukdetect_all_hits_unrarefied.txt")
df.500fg <- read_tsv(file = "results/500fg_eukdetect_all_hits_unrarefied.txt")

## Combine
euks_df <- rbind(df.1000ibd, df.500fg)

## Remove extra columns
#euks_df <- euks_df %>%
#  select(sampleid, Taxid, Read_counts)

## Remove brackets from Candida glabrata
euks_df <- euks_df %>% 
  mutate(Name = replace(Name, Name == "[Candida] glabrata", "Candida glabrata"))

## Change column names
tables2 <- euks_df %>%
  rename("Sample_ID" = sampleid,
  "Taxa_ID" = Taxid,
  "Taxa" = Name)



## Save
write.table(tables2, file = "figures_and_tables/Table_S2_eukdetect_results_unrarefied.csv", 
            sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)


## Change column names
euks_df <- euks_df %>%
  rename("taxid" = "Taxid") %>%
  rename("read_count" = "Read_counts")

## Select smaller df
euks_df <- euks_df %>%
  dplyr::select(sampleid, taxid, read_count)

# Pivot table of hits to OTU table format
otu_euks <- pivot_wider(euks_df, names_from = sampleid, values_from=read_count, values_fill = 0)
otu_euks <- otu_euks %>% remove_rownames %>% column_to_rownames(var="taxid") # make first column row names


# Import taxa table
tax_df <- read.csv(file = "metadata_files/eukdetect_taxa_table.csv", sep = ",")
tax_df <-tax_df %>% remove_rownames %>% column_to_rownames(var="taxid")
# Import metadata
md <- read.csv(file = "metadata_files/cohort_metadata_combined.csv", sep = ",")
md <- md %>% remove_rownames %>% column_to_rownames(var="sampleid")
sampledata <- sample_data(data.frame(md))
sampledata

# Import OTU table to phyloseq
otu_ps <- otu_table(as.matrix(otu_euks), taxa_are_rows = TRUE)
otu_ps
# Import taxa table to phyloseq
tax_ps <- tax_table(as.matrix(tax_df))
tax_ps
# Import all to phyloseq object
physeq <- phyloseq(otu_ps, tax_ps, sampledata)
physeq
sample_names(physeq)[1:5] #check if ok
rank_names(physeq) #check if ok
taxa_names(physeq)

# Psmelt physeq into dataframe
physeq.df <- psmelt(physeq)
physeq.df

# Save physeq dataframe as table
write.table(physeq.df, "results/eukdetect_abundance_table_unrarefied.csv", 
            sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)

# Merge species hits to genus level
ps.genus <- tax_glom(physeq, taxrank="Genus", NArm=FALSE) # Merge by genera
ps.genus.df <- psmelt(ps.genus) # Make dataframe for genus-level hits only
write.table(ps.genus.df, "results/eukdetect_abundance_table_genus_unrarefied.csv", 
           sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)



# Exploring results -------------------------------------------------------

# Import metadata and hits
tax_df <- read.csv(file = "metadata_files/eukdetect_taxa_table.csv", sep = ",")
md_df <- read.csv(file = "metadata_files/cohort_metadata_combined.csv", sep = ",")



# Add taxa and patient metadata
df <- left_join(euks_df, tax_df, by = "taxid")
df <- left_join(df, md_df, by = "sampleid")

# Rename and add columns
df <- plyr::rename(df, c("read_count" = "Abundance"))
df$ibd_subtype <- df$diagnosis_last_record
df$ibd_subtype[df$disease == "control"] <- "Control"


## IBD samples ##
df.ibd <- df[(df$disease=="IBD"),] #Make new df for IBD patients only
df.ibd1 <- unique(df.ibd[c("sampleid", "Species")]) #Make new df of samples with unique results at the species level

length(unique(df.ibd1$sampleid)) #Total individuals with microbial eukaryotes
colSums(table(df.ibd1$sampleid, df.ibd1$Species)) #Species counts

sum(duplicated(df.ibd1$sampleid)) #Number of duplicates
dups.ibd <-df.ibd1[duplicated(df.ibd1$sampleid)|duplicated(df.ibd1$sampleid, fromLast=TRUE),] #Find duplicate IDs
dups.ibd
table(dups.ibd$sampleid) #Number of duplicates per sample



## Control samples ##
df.ctl <- df[(df$disease=="control"),] #Make new df for controls only
df.ctl1 <- unique(df.ctl[c("sampleid", "Species")]) #Make new df of samples with unique results at the species level

length(unique(df.ctl1$sampleid)) #Total individuals with microbial eukaryotes
colSums(table(df.ctl1$sampleid, df.ctl1$Species)) #Species counts

sum(duplicated(df.ctl1$sampleid)) #Number of duplicates
dups.ctl <-df.ctl1[duplicated(df.ctl1$sampleid)|duplicated(df.ctl1$sampleid, fromLast=TRUE),] #Find duplicate IDs
dups.ctl
table(dups.ctl$sampleid) #Number of duplicates per sample
