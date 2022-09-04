## Logistic regression for EukDetect results ##

library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)

setwd("R_analysis/")

# 1. Data wrangling ----------------------------------------------------------

# Creating physeq table that includes all samples with 0 hits, as input for regression analysis

# Identify samples with 0 hits
hits <- read.csv(file = "results/eukdetect_results_rarefied.csv")
hits <- hits %>%
  dplyr::rename(taxid = Taxa_ID,
         read_count = Read_counts) %>%
  dplyr::select(sampleid, taxid, read_count)

all_samples <- read.csv(file = "metadata_files/cohort_metadata_combined.csv")
no_hits <- anti_join(all_samples, hits, by = "sampleid")

# Create table of 0 hits with arbitrary taxa ID 559292 in taxid column
df1 <- no_hits[1]
df1$taxid <- 559292
df1$read_count <- 0

# Append table of 0 hits to table with hits
df_euks <- bind_rows(hits, df1)

# Create total physeq dataframe that includes all 0 hit samples
otu_euks <- pivot_wider(df_euks, names_from = sampleid, values_from=read_count, values_fill = 0)
otu_euks <- otu_euks %>% remove_rownames %>% column_to_rownames(var="taxid") # make first column row names

# Import taxa table
tax_df <- read.csv(file = "metadata_files/eukdetect_taxa_table.csv", sep = ",")
tax_df <-tax_df %>% remove_rownames %>% column_to_rownames(var="taxid")

# Import metadata
md <- read.csv(file = "metadata_files/cohort_metadata_combined.csv", sep = ",")
md <- md %>% remove_rownames %>% column_to_rownames(var="sampleid")
sampledata <- sample_data(data.frame(md))

# Import OTU table to phyloseq
otu_ps <- otu_table(as.matrix(otu_euks), taxa_are_rows = TRUE)

# Import taxa table to phyloseq
tax_ps <- tax_table(as.matrix(tax_df))

# Import all to phyloseq object
physeq <- phyloseq(otu_ps, tax_ps, sampledata)
sample_names(physeq)[1:5] #check if ok
rank_names(physeq) #check if ok
taxa_names(physeq)

# Psmelt physeq into dataframe
physeq.df <- psmelt(physeq)
physeq.df

write.table(physeq.df, "results/eukdetect_abundance_table_rarefied_incl_no_hits.csv", 
            sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)

# Make genus-level dataframe
ps.genus <- tax_glom(physeq, taxrank="Genus", NArm=FALSE) # Merge by genera
ps.genus.df <- psmelt(ps.genus) # Make dataframe for genus-level only
write.table(ps.genus.df, "results/eukdetect_abundance_table_genus_rarefied_incl_no_hits.csv", 
            sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)



# 2. Create variables for regression -----------------------------------------

# Genus-only taxa dataframe 
physeq <- read.csv(file = "results/eukdetect_abundance_table_genus_rarefied_incl_no_hits.csv")

names(physeq)
attach(physeq)

# Make categorical variable of BMI
physeq$bmi_class
physeq[which(physeq$bmi<=18.5), 'bmi_class'] <- "Underweight"
physeq[which(physeq$bmi>18.5 & physeq$bmi<25), 'bmi_class'] <- "Normal"
physeq[which(physeq$bmi>=25 & physeq$bmi<30), 'bmi_class'] <- "Overweight"
physeq[which(physeq$bmi>=30), 'bmi_class'] <- "Obese"

# Make new column for IBD subtype
physeq$ibd_subtype <- physeq$diagnosis_last_record
physeq$ibd_subtype[physeq$disease == "control"] <- "no_IBD"

# Convert relevant columns to factors
physeq <- within(physeq, {
  disease <- factor(disease)
  ibd_subtype <- factor(ibd_subtype)
  sex <- factor(sex)
  smoking_status <- factor(smoking_status, levels = c("never", "past", "current"))
  bmi_class <- factor(bmi_class, levels = c("Normal", "Underweight", "Overweight", "Obese"))
})
summary(physeq)
attach(physeq)

# Saccharomyces subset dataframe
sacc <- subset(physeq, Genus == "Saccharomyces")
summary(sacc)
colSums(table(sacc$Abundance, sacc$disease)) #check if okay

# Blastocystis subset dataframe
blasto <- subset(physeq, Genus == "Blastocystis")
summary(blasto)
colSums(table(blasto$Abundance, blasto$disease)) #check if okay



# Regression --------------------------------------------
library(pscl)
## Load MASS for neg bin GLM
library(MASS)

# IBD vs no IBD - negative binomial regressions
# Saccharomyces
sacc.y <- formula(Abundance ~  disease + sex + age + bmi_class + smoking_status)
sacc.zb <- zeroinfl(sacc.y, data = sacc, dist = "negbin")
summary(sacc.zb)

# Blastocystis
blasto.y <- formula(Abundance ~  disease + sex + age + bmi_class + smoking_status)
blasto.zb <- zeroinfl(blasto.y, data = blasto, dist = "negbin")
summary(blasto.zb)



# Robust variance estimation ---------------------------------------------
## Load sandwich package for robust sandwich covariance matrix estimators
library(sandwich)
## Load lmtest package for coeftest
library(lmtest)


## IBD vs controls ##
# Saccharomyces
sink("results/regression_saccharomyces_all_ibd.txt")
print(coeftest(sacc.zb, vcov = sandwich))
sink()
# Blastocystis
sink("results/regression_blastocystis_all_ibd.txt")
print(coeftest(blasto.zb, vcov = sandwich))
sink()

## Final table formatted in Word.

