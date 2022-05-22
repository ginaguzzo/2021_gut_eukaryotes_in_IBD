## EukDetect FMT workflow ##
## This is for pre-rarefied data (i.e. absolute abundances)

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(readr)

setwd("R_analysis/")


# 1 Data wrangling ----------------------------------------------------------
# Import metadata and hits
euks_df <- read_tsv(file = "results/fmt_eukdetect_all_hits_unrarefied.txt")


## Change column names
tables5 <- euks_df %>%
  rename("Taxa_ID" = "Taxid") %>%
  rename("Taxa" = "Name")


## Save
write.table(tables5, file = "figures_and_tables/Table_S5_eukdetect_fmt_data_results_unrarefied.csv", 
            sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)

## Format table
euks_df <-
  rename(euks_df, "taxid" = "Taxid", "read_count" = "Read_counts")
tax_df <- read.csv(file = "metadata_files/eukdetect_fmt_taxa_table.csv", sep = ",")
md_df <- read.csv(file = "metadata_files/fmt_metadata.csv", sep = ",")


# Identify samples with 0 hits
no_hits <- anti_join(md_df, euks_df, by = "sampleid")

# Create table of 0 hits with arbitrary taxa ID 559292 (S. cerevisiae) in taxid column
df1 <- no_hits[1]
df1$taxid <- 559292
df1$read_count <- 0

# Append table of 0 hits to table with hits
euks_df <- bind_rows(euks_df, df1)

# Add taxa and patient metadata
df <- left_join(euks_df, tax_df, by = "taxid")
df <- left_join(df, md_df, by = "sampleid")


# Save dataframe as table
write.table(df, "results/fmt_eukdetect_abundance_table_unrarefied_incl_no_hits.csv", 
            sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)

# Convert column to binary
df <- plyr::rename(df, c("read_count" = "Abundance"))
df.bin <- df
df.bin$Abundance <- ifelse(df.bin$Abundance > 0, 1, 0)


# Group by pre or post-FMT
df.bin$timepoint_grouped <- ifelse(df.bin$Timepoint %in% c('W-2','D0'), 'pre',
                                   ifelse(df.bin$Timepoint %in% c('W2','W6','W10'), 'W2-10',
                                          ifelse(df.bin$Timepoint %in% c('W14','W18','W24'), 'W14-24', 'donor')))

dplyr::select(df.bin, "Timepoint", "timepoint_grouped") #check if okay

# Explore data
df.bin %>%
  group_by(ID_sample, FMT_SHAM) %>%
  group_by(FMT_SHAM) %>%
  summarise(n()) #Number of total samples in each treatment group
groups <- df.bin %>%
  group_by(ID_sample, FMT_SHAM) %>%
  summarise(n()) #Number of samples for each patient ID
groups %>%
  group_by(FMT_SHAM) %>%
  summarise(n()) #Number of individuals per group
df.bin %>%
  group_by(ID_sample, donor) %>%
  summarise() %>%
  print(n=30) #List of donors assigned to each patient

# Remove patients and donors with no hits, 
# only leaving 0's for time points in donors and samples that have hits at other time points
table(df.bin$Abundance, df.bin$ID_sample)
samples_to_remove <- c("1_19_BC", "1_22_LE", "1_5_JH", "1_7_LL", "3_1_MF", "7_2_LA", "7_34_KB", "7_55_PB", "7_59_BM")
df.bin <- filter(df.bin, !(ID_sample %in% samples_to_remove))

# Make new column with patient IDs that are simplified
df.bin$new_sampleid <- df.bin$ID_sample
df.bin$new_sampleid <- sub( "^1_","", df.bin$new_sampleid)
df.bin$new_sampleid <- sub( "_.{2}$","", df.bin$new_sampleid)
df.bin$new_sampleid <- sub( "^7_","", df.bin$new_sampleid)

# Reorder by patient number and FMT active vs sham/placebo
table(df.bin$FMT_SHAM, df.bin$new_sampleid) #Identify samples with FMT vs placebo
# FMT: "1", "4", "9", "15", "17", "21", Placebo: "2", "8", "10", "14", "18", "20", Donor: "47"
table(df.bin$Outcome, df.bin$new_sampleid) #Identify samples with success vs failure of FMT
# Success FMT: "1", "4", "15", "21", Failure FMT: "9", "17"
# Success placebo: "8", "18", Failure placebo: "2", "10", "14", "20"


# Update column with individual IDs for heatmap
df.bin$new_sampleid <- ifelse(df.bin$donor_patient=="donor", sub( "^","D", df.bin$new_sampleid),
                              ifelse(df.bin$donor_patient=="patient", sub( "^","P", df.bin$new_sampleid), NA))


df.bin$new_sampleid <- factor(df.bin$new_sampleid, order = TRUE, 
                              levels = c("P1", "P4", "P15", "P21", "P9", "P17",
                                         "P8", "P18", "P2", "P10", "P14", "P20",
                                         "D47")) #Reorder by clustering FMT vs placebo and success vs failure


# Reorder time points
df.bin$timepoint_grouped <- factor(df.bin$timepoint_grouped, order = TRUE, 
                                   levels = c("pre",
                                              "W2-10",
                                              "W14-24",
                                              "donor"))

# Make categorical variable for heatmap 
df.bin$hm_value <- ifelse(df.bin$timepoint_grouped %in% c("pre") & df.bin$Abundance==1, "Present pre-FMT",
                          ifelse(df.bin$timepoint_grouped %in% c("W2-10", "W14-24") & df.bin$Abundance==1, "Present post-FMT",
                                 ifelse(df.bin$timepoint_grouped %in% c("donor") & df.bin$Abundance==1, "Present in donor", "Absent")))

# Make nested facet categories
# FMT status
df.bin$facet_status <- df.bin$FMT_SHAM
df.bin$facet_status[df.bin$facet_status=="donor"] <- ""
df.bin$facet_status[df.bin$facet_status=="FMT"] <- "Active FMT group"
df.bin$facet_status[df.bin$facet_status=="SHAM"] <- "Sham FMT group"

# FMT outcome
df.bin$facet_outcome <- df.bin$Outcome
df.bin$facet_outcome <- ifelse(df.bin$facet_outcome %in% c('success_FMT','success_SHAM'), 'Success',
                               ifelse(df.bin$facet_outcome %in% c('Failure_FMT','failure_SHAM'), 'Failure', ''))
df.bin$facet_outcome[df.bin$FMT_SHAM=="donor"] <- ""

# Reorder levels based on outcome
df.bin$facet_outcome <- factor(df.bin$facet_outcome, order = TRUE, 
                               levels = c("Success", "Failure", "NA"))

df.bin$facet_status <- factor(df.bin$facet_status, order = TRUE, 
                              levels = c("Active FMT group", "Sham FMT group", ""))


# Custom y labels
colSums(table(df.bin %>% 
                group_by(Abundance, Species) %>%
                summarize(n()))) #Calculate number of samples in each group

my_y_titles <- rev(c(
  expression(paste(italic("Blastocystis"), ~"sp. subtype 1")),
  expression(paste(italic("Blastocystis"), ~"sp. subtype 2")),
  expression(paste(italic("Blastocystis"), ~"sp. subtype 4")),
  expression(paste(italic("Candida albicans"))),
  expression(paste(italic("Candida tropicalis"))),
  expression(paste(italic("Clavispora lusitaniae"))),
  expression(paste(italic("Debaryomyces hansenii"))),
  expression(paste(italic("Diutina catenulata"))),
  expression(paste(italic("Kluyveromyces lactis"))),
  expression(paste(italic("Kluyveromyces marxianus"))),
  expression(paste(italic("Penicillium nalgiovense"))),
  expression(paste(italic("Penicillium roqueforti"))),
  expression(paste(italic("Penicillium"), ~"sp.")),
  expression(paste(italic("Pichia fermentans"))),
  expression(paste(italic("Pichia kluyveri"))),
  expression(paste(italic("Pichia kudriavzevii"))),
  expression(paste(italic("Saccharomyces cerevisiae"))),
  expression(paste(italic("Torulaspora delbrueckii")))
)) 




# 2. Heatmap -----------------------------------------------------------------
library(ggh4x)


hm <- ggplot(df.bin, aes(timepoint_grouped, y = Species, fill = hm_value)) + 
  geom_tile(aes(width=0.9, height=0.9)) +
  theme_bw() +
  facet_nested(~ facet_status + facet_outcome + new_sampleid, scales = "free_x", space = "free_x") +
  scale_y_discrete(labels = my_y_titles,
                   limits = rev(levels(as.factor(df.bin$Species)))) +
  theme(
    plot.title = element_text(hjust = -0.12, vjust = -23, colour = "black", size = 12, face = "bold"),
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 12, colour = "black", face = "bold", vjust = 2),
    axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 10, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    panel.grid.major.x = element_blank(),
    legend.key = element_rect(color="black"),
    legend.text = element_text(colour = "black", size = 12),
    legend.title = element_text(face = "bold", colour = "black", size = 12),
    panel.border = element_rect(colour = "black"),
    strip.placement = "top",
    strip.background = element_rect(fill = c("#EEEEEE", "black"), color = "#FFFFFF"),
    strip.text = element_text(colour = "black", size = 12),
    panel.spacing = unit(0.5, "lines")) +
  scale_fill_manual(name = expression(bold(atop("Presence of", "eukaryotes"))), 
                    values=c("Present pre-FMT" = "#94DBFF", 
                             "Present post-FMT" = "#0073c2",
                             "Present in donor" = "#ff7f00",
                             "Absent"= NA), na.value = NA) +
  ggtitle("ID number →")


hm



# Modify grid grobs
library(grid)
library(gtable)

g.hm <- ggplot_gtable(ggplot_build(hm))
g.hm$widths
# Add spacing between select panels
g.hm$widths[27] <- 4*g.hm$widths[27]
g.hm$widths[51] <- 4*g.hm$widths[51]
# Remove extra donor facet strips
s1 <- which(g.hm$layout$name=="strip-t-3")
g.hm$grobs[[s1]]$grobs[[1]]$children[[2]]$children[[1]]$label <- ""
g.hm$grobs[[s1]]$grobs[[1]]$children[[1]]$gp$fill <- NA
s2 <- which(g.hm$layout$name=="strip-t-8")
g.hm$grobs[[s2]]$grobs[[1]]$children[[2]]$children[[1]]$label <- ""
g.hm$grobs[[s2]]$grobs[[1]]$children[[1]]$gp$fill <- NA


# Change facet strip colours
strip_both <- which(grepl('strip-', g.hm$layout$name))
fills <- c("black", "black", NA, "#acacac", "#acacac", "#acacac", "#acacac", NA,
           "#EEEEEE", "#EEEEEE", "#EEEEEE", "#EEEEEE", "#EEEEEE", "#EEEEEE", 
           "#EEEEEE", "#EEEEEE", "#EEEEEE", "#EEEEEE", "#EEEEEE", "#EEEEEE", "#EEEEEE")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g.hm$grobs[[i]]$grobs[[1]]$childrenOrder))
  g.hm$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}


# Change top facet trip font colour
s3 <- which(g.hm$layout$name=="strip-t-1")
g.hm$grobs[[s3]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- "white"
s4 <- which(g.hm$layout$name=="strip-t-2")
g.hm$grobs[[s4]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- "white"

grid.draw(g.hm)


ggsave(g.hm, filename = "figures_and_tables/Figure_3_fmt_heatmap_species.tiff",
       height = 8, width = 15, dpi = 300, units = "in", device='tiff')

