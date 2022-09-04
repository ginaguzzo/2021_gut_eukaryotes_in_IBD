# Generating tables and plots to compare RiboTagger and EukDetect identifications of microbial 
# eukaryotes in gut microbiome samples

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)

setwd("R_analysis/")

# 1.1 Data wrangling ----------------------------------------------------------

# Import ribotagger and eukdetect tables of hits
# RiboTagger results
df.ribo <- read.csv(file = "figures_and_tables/Table_S1_ribotagger_results.csv", sep = ",")
df.ribo <- df.ribo %>%
  rename("sampleid" = Sample.ID,
         "Disease" = Disease.status)
df.ribo <- unique(df.ribo[c("sampleid", "Genus", "Disease")]) 
df.ribo <- df.ribo[!(df.ribo$Genus=="Sinapis"),] #Remove samples with plant taxa, Sinapis
df.ribo <- df.ribo[!(df.ribo$Genus==""),] #Only keep samples with genus-level taxa identified
df.ribo <- df.ribo %>% 
  group_by(Genus, Disease) %>%
  summarize(n())
names(df.ribo)[3] <- c("Ribotagger")


# EukDetect results
df.ed <- read.csv(file = "results/eukdetect_abundance_table_genus_unrarefied.csv", sep = ",")
df.ed <- df.ed[!(df.ed$Abundance==0),] #Remove samples with 0 abundance
df.ed <- unique(df.ed[c("Sample", "Genus", "disease")]) 
df.ed <- df.ed[!(df.ed$Genus==""),] #Only keep samples with genus-level taxa identified
df.ed <- df.ed %>% 
  group_by(Genus, disease) %>%
  summarize(n())
names(df.ed)[1:3] <- c("Genus", "Disease", "Eukdetect")


# Group RiboTagger and EukDetect hits together into IBD vs control dataframes
df.total <- bind_rows(df.ribo, df.ed)
df.total[is.na(df.total)] <- 0 
df.total <- df.total %>%
  group_by(Genus, Disease) %>%
  summarise(Ribotagger=sum(Ribotagger),
            Eukdetect=sum(Eukdetect))
#write.table(df.total, file = "figures_and_tables/ribotagger_vs_eukdetect_unrarefied.csv", 
#            sep=",", quote = FALSE, col.names=TRUE, row.names=FALSE)

# Make separate tables for IBD and controls
ibd.z <- df.total %>% 
  filter(Disease=="IBD")
ctl.z <- df.total %>% 
  filter(Disease=="control")

# Add missing genera for each dataframe
ibd.miss <- data.frame(Genus = c("Giardia", "Hanseniaspora", "Pichia"),
                       Disease = "IBD",
                       Ribotagger = 0,
                       Eukdetect = 0)
ctl.miss <- data.frame(Genus = c("Clavispora", "Galactomyces", "Meyerozyma", "Nakaseomyces", "Wickerhamomyces"),
                       Disease = "control",
                       Ribotagger = 0,
                       Eukdetect = 0)

ibd.z <- rbind(ibd.z, ibd.miss)
ctl.z <- rbind(ctl.z, ctl.miss)


# Reorder alphabetically by genus
ibd.z <- with(ibd.z, ibd.z[order(Genus),])
ctl.z <- with(ctl.z, ctl.z[order(Genus),])

# Change zeros to NAs for plots
ibd <- ibd.z
ibd[ibd == 0] <- NA
ctl <- ctl.z
ctl[ctl == 0] <- NA 


# 2.1 Fig 1C - Pyramid lollipop plot ---------------------------------------------------
library(grid)
library(gridExtra)


g1 <- ggplot(ibd) +
  theme_bw() +
  geom_segment(data = ibd.z, aes(x=Genus, xend=Genus, y=Ribotagger, yend=Eukdetect), 
               color="black", size = 0.5) +
  geom_point( aes(x=Genus, y=Ribotagger, color="Ribotagger"), size=7.5) +
  geom_point( aes(x=Genus, y=Eukdetect, color="Eukdetect"), size=7.5) +
  scale_color_manual(values = c("#0073c2", "black"),
                     guide  = guide_legend(), 
                     name   = "",
                     labels = c("EukDetect", "RiboTagger")) +
  geom_text(aes(x=Genus, y=Eukdetect, label=Eukdetect), color = "white", size = 3) +
  geom_text(aes(x=Genus, y=Ribotagger, label=Ribotagger), color = "white", size = 3) +
  geom_hline(yintercept = -1.95, colour = "black") +
  ggtitle(expression(atop(paste(bold("IBD")), "(n = 355)"))) +
  theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold",  size = 12),
        axis.title.x = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x = element_text(colour = "black"),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = "black"),
        panel.border = element_blank(),
        plot.margin = unit(c(1,2,1,1), "mm"),
        legend.position = "left",
        legend.text = element_text(colour = "black"),
        legend.text.align=0) +
  scale_y_reverse(limits = c(71, -2), expand = c(0, 0)) +
  scale_x_discrete(limits = rev(factor(ibd$Genus))) +
  coord_flip() +
  ylab("Number of hits") +
  xlab(NULL) 

g1


g2 <- ggplot(ctl) +
  theme_bw() +
  geom_segment(data = ctl.z, aes(x=Genus, xend=Genus, y=Ribotagger, yend=Eukdetect), 
               color="black", size = 0.5) +
  geom_point( aes(x=Genus, y=Ribotagger), color=("black"), size=7.5) +
  geom_point( aes(x=Genus, y=Eukdetect), color=("#0073c2"), size=7.5) +
  geom_text(aes(x=Genus, y=Eukdetect, label=Eukdetect), color = "white", size = 3) +
  geom_text(aes(x=Genus, y=Ribotagger, label=Ribotagger), color = "white", size = 3) +
  ggtitle(expression(atop(paste(bold("Control subjects")), "(n = 471)"))) +
  theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
        axis.title.x = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(1,2,1,1), "mm")
        ) +
  scale_y_continuous(limits = c(-3,151), expand = c(0, 0)) +
  scale_x_discrete(limits = rev(factor(ctl$Genus))) +
  coord_flip() +
  ylab("Number of hits") +
  xlab(NULL)
  

g2


g.mid <- ggplot(ibd, aes(x=1, y = reorder(Genus, desc(Genus)))) +
  geom_text(aes(label=Genus, fontface = "italic")) +
  geom_segment(aes(x=0.1,xend=0.2,yend=Genus)) +
  geom_segment(aes(x=1.1,xend=1.2,yend=Genus)) +
  ylab(NULL)+
  ggtitle(expression(atop(paste("", bold("Genus"))))) +
  scale_x_continuous(expand=c(-1,1), limits=c(0.94,1.065)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -8, colour = "black", size = 12, face = "bold"),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,5,-1), "mm"))

g.mid


# Arrange plots together
gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))

gg.full <- grid.arrange(gg1, gg.mid, gg2, ncol=3, widths=c(4/9,1/10,4/9)) 
gg.full 
#ggsave(gg.full, filename = "figures_and_tables/ribotagger_vs_eukdetect_lollipop_plot.tiff",
#       height = 6, width = 15, dpi = 300, units = "in", device='tiff')




# 3.1 Bar plots: IBD vs controls -------------------------------------------------------

# Create stacked bar plots comparing eukaryote hits to total group sizes in cohort


# 3.1.1 Preparing data for plots ----------------------------------------------------

# Calculate total number of individuals with hits
# RiboTagger
df.ribo <- read.csv(file = "figures_and_tables/Table_S1_ribotagger_results.csv", sep = ",")
df.ribo <- unique(df.ribo[c("sampleid", "f", "g", "disease")]) 
df.ribo <- df.ribo[!(df.ribo$f=="" & df.ribo$g==""),] #Remove samples with no family and no genus listed
df.ribo <- df.ribo[!(df.ribo$f=="Liliopsida" | df.ribo$g=="Sinapis"),] #Remove samples with plant taxa identified
df.ribo <- unique(df.ribo[c("sampleid", "disease")]) 
df.ribo <- df.ribo %>% 
  group_by(disease) %>%
  summarise(n())
names(df.ribo)[1:2] <- c("Cohort", "number")
df.ribo$euk_hits <- "euks"

df.r.nohits <- df.ribo
df.r.nohits$euk_hits <- "no_euks"
df.r.nohits$number[1] <- 471 - df.r.nohits$number[1]
df.r.nohits$number[2] <- 355 - df.r.nohits$number[2]

df.ribo <- bind_rows(df.ribo, df.r.nohits)
df.ribo$tool <- "RT"

# EukDetect
df.ed <- read.csv(file = "results/eukdetect_abundance_table_genus_unrarefied.csv", sep = ",")
df.ed <- df.ed[!(df.ed$Abundance==0),] #Remove samples with 0 abundance
df.ed <- unique(df.ed[c("Sample", "disease")]) 
df.ed <- df.ed %>% 
  group_by(disease) %>%
  summarize(n())
names(df.ed)[1:2] <- c("Cohort", "number")
df.ed$euk_hits <- "euks"

df.e.nohits <- df.ed
df.e.nohits$euk_hits <- "no_euks"
df.e.nohits$number[1] <- 471 - df.e.nohits$number[1]
df.e.nohits$number[2] <- 355 - df.e.nohits$number[2]

df.ed <- bind_rows(df.ed, df.e.nohits)
df.ed$tool <- "ED"

# Merge tables
df.total2 <- bind_rows(df.ribo, df.ed)
df.total2$Cohort[df.total2$Cohort=="control"] <- "Control"

# Reorder by cohort
df.total2$Cohort <- factor(df.total2$Cohort, order = TRUE, 
                                levels = c("Control", "IBD"))

# Reorder by no hits
df.total2$euk_hits <- factor(df.total2$euk_hits, order = TRUE, 
                            levels = c("no_euks", "euks"))

# Reorder by no hits
df.total2$tool <- factor(df.total2$tool, order = TRUE, 
                            levels = c("RT", "ED"))


# Make column for specific bar chart categories
df.total2$bar_fill <- df.total2$euk_hits
df.total2$bar_fill <- ifelse(df.total2$euk_hits=="euks" & df.total2$tool=="RT", "euks_RT",
                             ifelse(df.total2$euk_hits=="euks" & df.total2$tool=="ED", "euks_ED",
                                    ifelse(df.total2$euk_hits=="no_euks", "no_euks", "")))

df.total2$bar_fill <- factor(df.total2$bar_fill, order = TRUE, 
                         levels = c("no_euks", "euks_RT", "euks_ED"))

## Facet labels for plot
facet_labels1 = list(
  "Controls" = expression(atop(paste(bold("Control subjects")), "(n = 471)")),
  "IBD" = expression(atop(paste(bold("IBD")), "(n = 355)")))
vlabeller1 <- function(variable,value){
  return(facet_labels1[value])
}


# 3.1.2 Stacked bar plot: proportional -------------------------------------

df.total2 <- df.total2 %>% 
  mutate(prop = if_else(Cohort == "IBD", number/355, number/471))

b1.1 <- ggplot(df.total2, aes(x=tool, y=number, fill = bar_fill)) +
  geom_bar(stat="identity", position="fill", color="black", width=0.8) +
  facet_grid(~Cohort, switch = "x", labeller = vlabeller1) +
  scale_y_continuous(labels = scales::percent,  expand = c(0,0.015)) +
  geom_text(aes(label = scales::percent(prop), y=prop), position = position_stack(vjust = 0.5), colour = "white") +
  ylab("Percentage of cohort") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        legend.text = element_text(colour= "black" , size = 12),
        legend.title = element_text(colour = "black", face = "bold", size = 12),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(1, "cm"),
        strip.text = element_text(colour = "black", size = 12),
        strip.placement = "outside",
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        plot.margin = unit(c(0,1,2,0), "cm")) +
  scale_fill_manual(name = expression(bold(atop("Presence of", "eukaryotes"))),
                    values = c("#acacac", "black", "#0073c2"), 
                    labels = c("Absent", "Present - RiboTagger", "Present - EukDetect"))

b1.1


#ggsave(b1.1, filename = "figures_and_tables/ribotagger_vs_eukdetect_stacked_barplot_percent.tiff",
#       height = 8, width = 7, dpi = 300, units = "in", device='tiff')



# 3.2 Bar plots: IBD subtypes -------------------------------------

# Create stacked bar plots comparing eukaryote hits in IBD subtypes to total group sizes in cohort


# 3.2.1 Preparing data for plots ------------------------------------------------

## IBD subtypes ##
# RiboTagger results
df.ribo <- read.csv(file = "figures_and_tables/Table_S1_ribotagger_results.csv", sep = ",")
df.ribo <- unique(df.ribo[c("sampleid", "f", "g", "disease")]) 
df.ribo <- df.ribo[!(df.ribo$f=="" & df.ribo$g==""),] #Remove samples with no family and no genus listed
df.ribo <- df.ribo[!(df.ribo$disease=="control"),] #Remove control samples
df.ribo <- df.ribo[!(df.ribo$f=="Liliopsida" | df.ribo$g=="Sinapis"),] #Remove samples with plant taxa identified
df.ribo <- unique(df.ribo[c("sampleid", "disease")]) 

# EukDetect results
df.ed <- read.csv(file = "results/eukdetect_abundance_table_genus_unrarefied.csv", sep = ",")
df.ed <- df.ed[!(df.ed$Abundance==0),] #Remove samples with 0 abundance
df.ed <- df.ed[!(df.ed$disease=="control"),] #Remove control samples
df.ed <- unique(df.ed[c("Sample", "disease")]) 
names(df.ed)[names(df.ed) == 'Sample'] <- 'sampleid' #Rename column before joining

# Add disease subtype from metadata
md <- read.csv(file = "metadata_files/cohort_metadata_combined.csv", sep = ",")
md <- select(md, sampleid, diagnosis_last_record)

df.ribo <- inner_join(df.ribo, md, by="sampleid")
df.ed <- inner_join(df.ed, md, by="sampleid")


# Clean up tables for bar plot
# RiboTagger table
df.ribo <- df.ribo %>% 
  group_by(diagnosis_last_record) %>%
  summarise(n())
names(df.ribo)[1:2] <- c("IBD subtype", "number")
df.ribo$euk_hits <- "euks"

df.r.nohits <- df.ribo
df.r.nohits$euk_hits <- "no_euks"
df.r.nohits$number[1] <- 206 - df.r.nohits$number[1]
df.r.nohits$number[2] <- 23 - df.r.nohits$number[2]
df.r.nohits$number[3] <- 126 - df.r.nohits$number[3]

df.ribo <- bind_rows(df.ribo, df.r.nohits)
df.ribo$tool <- "RT"

# EukDetect table
df.ed <- df.ed %>% 
  group_by(diagnosis_last_record) %>%
  summarize(n())
names(df.ed)[1:2] <- c("IBD subtype", "number")
df.ed$euk_hits <- "euks"

df.e.nohits <- df.ed
df.e.nohits$euk_hits <- "no_euks"
df.e.nohits$number[1] <- 206 - df.e.nohits$number[1]
df.e.nohits$number[2] <- 23 - df.e.nohits$number[2]
df.e.nohits$number[3] <- 126 - df.e.nohits$number[3]

df.ed <- bind_rows(df.ed, df.e.nohits)
df.ed$tool <- "ED"

# Merge tables
df.total3 <- bind_rows(df.ribo, df.ed)

# Reorder by subtype
df.total3$`IBD subtype` <- factor(df.total3$`IBD subtype`, order = TRUE, 
                           levels = c("CD", "UC", "IBDU"))

# Reorder by no hits
df.total3$euk_hits <- factor(df.total3$euk_hits, order = TRUE, 
                             levels = c("no_euks", "euks"))

# Reorder by no hits
df.total3$tool <- factor(df.total3$tool, order = TRUE, 
                         levels = c("RT", "ED"))


# Make column for specific bar chart categories
df.total3$bar_fill <- df.total3$euk_hits
df.total3$bar_fill <- ifelse(df.total3$euk_hits=="euks" & df.total3$tool=="RT", "euks_RT",
                             ifelse(df.total3$euk_hits=="euks" & df.total3$tool=="ED", "euks_ED",
                                    ifelse(df.total3$euk_hits=="no_euks", "no_euks", "")))

df.total3$bar_fill <- factor(df.total3$bar_fill, order = TRUE, 
                             levels = c("no_euks", "euks_RT", "euks_ED"))


## Facet labels for plot
facet_labels2 = list("CD" = expression(atop(paste(bold("CD")), "(n = 206)")),
                     "UC" = expression(atop(paste(bold("UC")), "(n = 126)")),
                     "IBDU" = expression(atop(paste(bold("IBDU")), "(n = 23)")))

vlabeller2 <- function(variable,value){
  return(facet_labels2[value])
}


# 3.2.2 Stacked bar plot: proportional -------------------------------------

df.total3 <- df.total3 %>% 
  mutate(prop = if_else(`IBD subtype` == "CD", number/206,
                        if_else(`IBD subtype` == "UC", number/126, number/23)))

b2.1 <- ggplot(df.total3, aes(x=tool, y=number, fill = bar_fill)) +
  geom_bar(stat="identity", position="fill", color="black", width=0.9) +
  facet_grid(~`IBD subtype`, switch = "x", labeller = vlabeller2) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0.015)) +
  geom_text(aes(label = scales::percent(prop), y=prop), position = position_stack(vjust = 0.5), colour = "white") +
  ylab("Percentage of cohort") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(1, "cm"),
        strip.text = element_text(colour = "black", size = 12),
        strip.placement = "outside",
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        plot.margin = unit(c(0,1,2,0), "cm")) +
  scale_fill_manual(values = c("#acacac", "black", "#0073c2"), 
                    labels = c("Eukaryotes absent", "Eukaryotes present - RiboTagger", "Eukaryotes present - EukDetect"))



b2.1


#ggsave(b2.1, filename = "figures_and_tables/ribotagger_vs_eukdetect_ibd_subtypes_stacked_barplot_percent.tiff",
#       height = 8, width = 6, dpi = 300, units = "in", device='tiff')




# 4.1 Arrange plots for figure -----------------------------------------------------------
library(patchwork)

patchwork <- (b1.1 | b2.1) /
              gg.full
patchwork <- patchwork + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 18, face = "bold"))
patchwork <- patchwork + plot_layout(heights = c(1.5, 2.7))
patchwork


ggsave(patchwork, filename = "figures_and_tables/Figure_1_ribotagger_vs_eukdetect_draft.tiff",
       height = 14, width = 15, dpi = 300, units = "in", device='tiff')


## Lines and boxes connecting Fig 1A and 1B were added in Adobe Illustrator.