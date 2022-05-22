## Comparing results of EukDetect before and after rarefying
## all samples to 1.7 million sequences.

library(dplyr)
library(tidyr)
library(tibble)


setwd("R_analysis/")


# 2. After rarefying ------------------------------------------------------

# Import hits table
df_euks <- read.csv(file = "results/eukdetect_results_rarefied.csv", sep = ",")

# Remove extra columns
df_euks <- df_euks %>%
  dplyr::select(sampleid, Taxa_ID, Read_counts)

df_euks <- df_euks %>% relocate(Taxa_ID)

# Pivot to OTU table format
otu_euks <- pivot_wider(df_euks, names_from = sampleid, values_from=Read_counts, values_fill = 0)
otu_euks <- otu_euks %>% remove_rownames %>% column_to_rownames(var="Taxa_ID") # make first column row names

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


# Make new column for IBD subtype
physeq.df$ibd_subtype <- physeq.df$diagnosis_last_record
physeq.df$ibd_subtype[physeq.df$disease == "control"] <- "Control"

# Save physeq dataframe as table
write.table(physeq.df, "results/eukdetect_abundance_table_rarefied.csv", 
            sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)

# Merge species hits to genus level
ps.genus <- tax_glom(physeq, taxrank="Genus", NArm=FALSE) # Merge by genera
ps.genus.df <- psmelt(ps.genus) # Make dataframe for genus-level hits only

# Make new column for IBD subtype
ps.genus.df$ibd_subtype <- ps.genus.df$diagnosis_last_record
ps.genus.df$ibd_subtype[ps.genus.df$disease == "control"] <- "Control"

write.table(ps.genus.df, "results/eukdetect_abundance_table_genus_rarefied.csv", 
            sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)


# 1.1 Data wrangling ----------------------------------------------------------

# Before rarefying
df1 <- read.csv(file = "results/eukdetect_abundance_table_genus_unrarefied.csv", sep = ",")
df1 <- df1[!(df1$Abundance==0),] #Remove samples with 0 abundance
df1 <- unique(df1[c("Sample", "Genus", "disease")]) 
df1 <- df1[!(df1$Genus==""),] #Only keep samples with genus-level taxa identified
df1 <- df1 %>% 
  group_by(Genus, disease) %>%
  summarize(n())
names(df1)[1:3] <- c("Genus", "Disease_status", "Abundance_before_rarefying")

# After rarefying
df2 <- read.csv(file = "results/eukdetect_abundance_table_genus_rarefied.csv", sep = ",")
df2 <- df2[!(df2$Abundance==0),] #Remove samples with 0 abundance
df2 <- unique(df2[c("Sample", "Genus", "disease")]) 
df2 <- df2[!(df2$Genus==""),] #Only keep samples with genus-level taxa identified
df2 <- df2 %>% 
  group_by(Genus, disease) %>%
  summarize(n())
names(df2)[1:3] <- c("Genus", "Disease_status", "Abundance_after_rarefying")

df <- bind_rows(df1, df2)
df[is.na(df)] <- 0 
df <- df %>%
  group_by(Genus, Disease_status) %>%
  summarise(Abundance_before_rarefying=sum(Abundance_before_rarefying),
            Abundance_after_rarefying=sum(Abundance_after_rarefying))

df$Difference_in_abundance <- df$Abundance_before_rarefying - df$Abundance_after_rarefying # Make column for difference of before and after rarefying

write.table(df, file = "figures_and_tables/Table_S4_eukdetect_before_vs_after_rarefying.csv", 
            sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)




