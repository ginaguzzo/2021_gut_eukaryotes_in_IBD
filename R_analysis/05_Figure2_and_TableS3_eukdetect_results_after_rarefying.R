# EukDetect results analysis after rarefying

library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(ggplot2)

setwd("R_analysis/")

# Make Table S3 ----------------------------------------------------------

## Import tables and combine
df.1000ibd <- read_tsv(file = "results/1000ibd_eukdetect_all_hits_rarefied.txt")
df.500fg <- read_tsv(file = "results/500fg_eukdetect_all_hits_rarefied.txt")

## Combine
euks_df <- rbind(df.1000ibd, df.500fg)

## Remove extra columns
#euks_df <- euks_df %>%
#  select(sampleid, Taxid, Read_counts)

## Remove brackets from Candida glabrata
euks_df <- euks_df %>% 
  mutate(Name = replace(Name, Name == "[Candida] glabrata", "Candida glabrata"))

## Change column names
tables3 <- euks_df %>%
  rename("Taxa_ID" = "Taxid") %>%
  rename("Taxa" = "Name")

## Save
write.table(tables3, file = "figures_and_tables/Table_S3_eukdetect_results_rarefied.csv", 
            sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)

# 1. Data wrangling ----------------------------------------------------------

## Change column names
euks_df <- euks_df %>%
  rename("taxid" = "Taxid") %>%
  rename("read_count" = "Read_counts")

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



# 1.1 Exploring results ---------------------------------------------------

## IBD samples ##
df.ibd <- df[(df$disease=="IBD"),] #Make new df for IBD patients only
df.ibd1 <- unique(df.ibd[c("sampleid", "Species", "ibd_subtype", "Genus")]) #Make new df of samples with unique results at the species level

length(unique(df.ibd1$sampleid)) #Total individuals with microbial eukaryotes
colSums(table(df.ibd1$sampleid, df.ibd1$Species)) #Species counts
table(df.ibd1$ibd_subtype, df.ibd1$Species) #Species counts per subtype
table(df.ibd1$ibd_subtype, df.ibd1$Genus) #Genera counts per subtype

sum(duplicated(df.ibd1$sampleid)) #Number of duplicates
dups.ibd <-df.ibd1[duplicated(df.ibd1$sampleid)|duplicated(df.ibd1$sampleid, fromLast=TRUE),] #Find duplicate IDs
dups.ibd
table(dups.ibd$sampleid) #Number of duplicates per sample

# Unique IDs with Blastocystis
df.ibd.blasto <- df[(df$disease=="IBD" & df$Genus=="Blastocystis" & df$ibd_subtype!="Control"),]
length(unique(df.ibd.blasto$sampleid, df.ibd.blasto$Species))
df.ibd.blasto[duplicated(df.ibd.blasto$sampleid)|duplicated(df.ibd.blasto$sampleid, fromLast=TRUE),]



## Control samples ##
df.ctl <- df[(df$disease=="control"),] #Make new df for controls only
df.ctl1 <- unique(df.ctl[c("sampleid", "Species")]) #Make new df of samples with unique results at the species level

length(unique(df.ctl1$sampleid)) #Total individuals with microbial eukaryotes
colSums(table(df.ctl1$sampleid, df.ctl1$Species)) #Species counts

sum(duplicated(df.ctl1$sampleid)) #Number of duplicates
dups.ctl <-df.ctl1[duplicated(df.ctl1$sampleid)|duplicated(df.ctl1$sampleid, fromLast=TRUE),] #Find duplicate IDs
dups.ctl
table(dups.ctl$sampleid) #Number of duplicates per sample

# Unique IDs with Blastocystis
df.ctl.blasto <- df[(df$disease=="control" & df$Genus=="Blastocystis"),]
length(unique(df.ctl.blasto$sampleid, df.ctl.blasto$Species))
df.ctl.blasto[duplicated(df.ctl.blasto$sampleid)|duplicated(df.ctl.blasto$sampleid, fromLast=TRUE),]


# 2. Plotting ----------------------------------------------------------------


# 2.1.1 Heatmap data prep ----------------------------------------------------------------

library(RColorBrewer)
library(scales)
library(grid)
library(gtable)

# List species
df %>% 
  group_by(Species) %>% 
  summarise(n()) 

# Order species alphabetically
df$Species <- factor(df$Species, order = TRUE, 
                            levels = c("Blastocystis hominis",
                                       "Blastocystis sp. subtype 1",
                                       "Blastocystis sp. subtype 2",
                                       "Blastocystis sp. subtype 3",
                                       "Blastocystis sp. subtype 4",
                                       "Candida albicans",
                                       "Candida glabrata",
                                       "Candida sake",
                                       "Clavispora lusitaniae",
                                       "Cyberlindnera jadinii NRRL Y-1542",
                                       "Debaryomyces hansenii",
                                       "Malassezia restricta",
                                       "Penicillium roqueforti",
                                       "Pichia kudriavzevii",
                                       "Saccharomyces cerevisiae",
                                       "Wickerhamomyces anomalus"))

# Order by disease status
df$ibd_subtype <- factor(df$ibd_subtype, order = TRUE, 
                                levels = c("CD",
                                           "UC",
                                           "IBDU",
                                           "Control"))


# Change 0 to NA
df.na <- df
df.na[df.na == 0] <- NA


## SPECIES-LEVEL HEATMAP ##
# Italicise species names for axis labels
my_y_titles <- rev(c(
  expression(paste(italic("Blastocystis hominis"))),
  expression(paste(italic("Blastocystis"), ~"sp. subtype 1")),
  expression(paste(italic("Blastocystis"), ~"sp. subtype 2")),
  expression(paste(italic("Blastocystis"), ~"sp. subtype 3")),
  expression(paste(italic("Blastocystis"), ~"sp. subtype 4")),
  expression(paste(italic("Candida albicans"))),
  expression(paste(italic("Candida glabrata"))),
  expression(paste(italic("Candida sake"))),
  expression(paste(italic("Clavispora lusitaniae"))),
  expression(paste(italic("Cyberlindnera jadinii"), ~"NRRL Y-1542")),
  expression(paste(italic("Debaryomyces hansenii"))),
  expression(paste(italic("Malassezia restricta"))),
  expression(paste(italic("Penicillium roqueforti"))),
  expression(paste(italic("Pichia kudriavzevii"))),
  expression(paste(italic("Saccharomyces cerevisiae"))),
  expression(paste(italic("Wickerhamomyces anomalus")))
))

# Add subtitles for facet grid labels
# Calculate total number of individuals with hits in each cohort group (disregards multiple hits/individual)
colSums(table(df %>% 
                group_by(sampleid, ibd_subtype) %>%
                summarize(n()))) #Calculate number of samples in each group
# List of labels
vnames <-list(
  "CD" = expression(atop(paste(bold("CD")), "(n=29)")),
  "UC" = expression(atop(paste(bold("UC")), "(n=20)")),
  "IBDU" = expression(atop(paste(bold("IBDU")), "(n=6)")),
  "Control" = expression(atop(paste(bold("Control subjects")), "(n=108)")))

# Function with labels
vlabeller <- function(variable,value){
  return(vnames[value])
}



# 2.1.2 Heatmap: species-level --------------------------------------------

g1 <- ggplot(df.na, aes(sampleid, Species, fill=Abundance)) + 
  geom_tile() +
  facet_grid(~ ibd_subtype, scales = "free_x", space = "free_x", switch = "x", labeller = vlabeller) +
  theme_bw() +
  scale_fill_gradient(name = "Abundance",
                      low = "#66CCFF",
                      high = "#000033",
                      na.value = NA) +
  scale_y_discrete(labels = my_y_titles,
                   limits = rev(levels(as.factor(df$Species)))) +
  theme(plot.title = element_text(hjust = -0.15, colour = "black", size = 12, face = "bold"),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.ticks.x = element_blank(),
        legend.text = element_text(colour = "black"),
        legend.title = element_text(face = "bold", colour = "black"),
        panel.border = element_rect(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"),
        strip.text = element_text(colour = "black", size = 12),
        panel.spacing = unit(0.5, "lines"),
        plot.margin = unit(c(0,0,1,0), "cm")) +
  ggtitle("Species")


g1

# Adjust label clipping by turning plot into grob
pg <- ggplotGrob(g1)
pg
pg$grobs[[18]]$layout$clip <- "off" #Turn off clipping for IBDU label

grid::grid.draw(pg)

# Save
#ggsave(pg, filename = "figures_and_tables/eukdetect_heatmap_species.tiff",
#       height = 8, width = 15, dpi = 300, units = "in", device='tiff')




# 2.2 Stacked bar plots: data prep ---------------------------------------------------------------

require(pals)

df %>% 
  group_by(Genus) %>% 
  summarise(n()) 

mycolours <- c(
  "Blastocystis" = "#0073c2",
  "Candida" = "#f2db80",
  "Clavispora" = "#868686",
  "Cyberlindnera" = "#f7c1b9",
  "Debaryomyces" = "#002d4b",
  "Malassezia" = "#92ccf4",
  "Nakaseomyces" = "#d08567",
  "Penicillium" = "#faf5e1",
  "Pichia" = "#722b20",
  "Saccharomyces" = "#EEEEEE",
  "Wickerhamomyces" = "#efbf00")


mycolours2 <- c(
  "Blastocystis" = "#0073c2",
  "Candida" = "#94DBFF",
  "Clavispora" = "#33a02c",
  "Cyberlindnera" = "#b2df8a",
  "Debaryomyces" = "#722b20",
  "Malassezia" = "#f9ca8e",
  "Nakaseomyces" = "#ff7f00",
  "Penicillium" = "#ffff99",
  "Pichia" = "#525252",
  "Saccharomyces" = "#acacac",
  "Wickerhamomyces" = "black")

mycolours2 <- c(
  "Blastocystis" = "#0073c2",
  "Candida" = "#94DBFF",
  "Clavispora" = "#21711C",
  "Cyberlindnera" = "#b2df8a",
  "Debaryomyces" = "#722b20",
  "Malassezia" = "#F9D7D2",
  "Nakaseomyces" = "#ff7f00",
  "Penicillium" = "#ffff99",
  "Pichia" = "#525252",
  "Saccharomyces" = "#acacac",
  "Wickerhamomyces" = "black")


# 2.2.1 Stacked bar plots: IBD vs control ----------------------------------
df.nz <- df[!(df$Abundance==0),] #Remove samples with 0 abundance
colSums(table(df.nz %>% 
                group_by(sampleid, disease) %>%
                summarize(n()))) #Calculate number of samples in each group

# Make custom x labels
my_x_titles = c(
  "IBD" = expression(atop(paste(bold("IBD")), "(n = 55)")),
  "control" = expression(atop(paste(bold("Control subjects")), "(n = 108)")))


# Stacked bar plot IBD vs control: proportional
b1.2 <- ggplot(df.nz, aes(x=disease, y=Abundance, fill = Genus)) +
  geom_bar(stat="identity", position="fill", color=NA, width=0.8) +
  ylab("Abundance") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", face = "bold", size = 12, vjust = 2),
        axis.text.x = element_text(colour = "black", size = 12, vjust = -1),
        axis.text.y = element_text(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        legend.title = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(colour = "black", face = "italic", size = 12),
        legend.position = "right",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = unit(c(0,0,2,0), "cm")) +
  scale_x_discrete(labels = my_x_titles) +
  scale_y_continuous(labels=percent, expand = c(0,0.015)) +
  scale_fill_manual(values = mycolours2)

b1.2


#ggsave(b1.2, filename = "figures_and_tables/eukdetect_stacked_bar_plot_genus_percent.tiff",
#       height = 8, width = 7, dpi = 300, units = "in", device='tiff')




# 2.2.2 Stacked bar plots: IBD subtype ------------------------------------------
df.nz.nc <- df.nz[!(df.nz$disease=="control"),] #Remove control samples
colSums(table(df.nz.nc %>% 
                group_by(sampleid, ibd_subtype) %>%
                summarize(n()))) #Calculate number of samples in each group

# Make custom x labels
my_x_titles2 = c(
  "CD" = expression(atop(bold(paste("CD")), "(n = 29)")),
  "UC" = expression(atop(paste(bold("UC")), "(n = 20)")),
  "IBDU" = expression(atop(paste(bold("IBDU")), "(n = 6)")))


# Stacked bar plot IBD subtype: proportional
b2.2 <- ggplot(df.nz.nc, aes(x=ibd_subtype, y=Abundance, fill = Genus)) +
  geom_bar(stat="identity", position="fill", color=NA, width=0.8) +
  ylab("Abundance") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", face = "bold", size = 12, vjust = 2),
        axis.text.x = element_text(colour = "black", size = 12, vjust = -1),
        axis.text.y = element_text(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,2,0), "cm")) +
  scale_x_discrete(labels = my_x_titles2) +
  scale_y_continuous(labels=percent, expand = c(0,0.015)) +
  scale_fill_manual(values = mycolours2)

b2.2


#ggsave(b2.2, filename = "figures_and_tables/eukdetect_stacked_bar_plot_genus_percent_ibd_subtype.tiff",
#       height = 8, width = 7, dpi = 300, units = "in", device='tiff')



# 2.3 Proportions numbers for results ---------------------------------------------------------------------

# Get the exact numbers of each taxon's proportions for the results section of the paper
df.nz %>% 
  group_by(disease, Genus) %>%
  summarise(total_ab = sum(Abundance, na.rm = TRUE)) %>%
  mutate(freq = total_ab / sum(total_ab)) #Taxa abundances in controls vs IBD

df.nz.nc %>% 
  group_by(ibd_subtype, Genus) %>%
  summarise(total_ab = sum(Abundance, na.rm = TRUE)) %>%
  mutate(freq = total_ab / sum(total_ab)) #Taxa abundances in IBD subtypes



# 2.4 Venn diagram --------------------------------------------------------

# Find shared and unique species for each group

# IBD and control only dataframes
ibd.df <- df[(df$disease == "IBD"),]
ibd.df <- ibd.df[!(ibd.df$Abundance==0),] 
ibd.df <- subset(ibd.df, select = c(Species))
ibd.df <- unique(ibd.df$Species)

ctl.df <- df[(df$disease == "control"),]
ctl.df <- ctl.df[!(ctl.df$Abundance==0),] 
ctl.df <- subset(ctl.df, select = c(Species))
ctl.df <- unique(ctl.df$Species)

# Unique species
setdiff(ibd.df, ctl.df) #Unique to IBD
setdiff(ctl.df, ibd.df) #Unique to controls

# Shared species
intersect(ibd.df, ctl.df)

# Venn diagram made manually in Adobe Illustrator.



# 3. Arrange plots for figure -------------------------------------------------------
library(patchwork)

patchwork <- (b1.2 | b2.2) /
  pg
  
patchwork <- patchwork + plot_annotation(tag_levels = 'A', ) &
  theme(plot.tag = element_text(size = 18, face = "bold"))
patchwork <- patchwork + plot_layout(heights = c(1.75, 2.75))
patchwork

ggsave(patchwork, filename = "figures_and_tables/Figure_2_eukdetect_results_rarefied_draft.tiff",
       height = 15, width = 15, dpi = 300, units = "in", device='tiff')

## Fig 2D Venn Diagram and lines connecting 2A to 2B were added in Adobe Illustrator.