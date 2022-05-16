## Calculating library sizes of study group samples

library(dplyr)
library(tibble)

setwd("R_analysis/")

# BEFORE SEQUENCE QUALITY CONTROL i.e. raw reads --------------------------

# Import data
ibd <- read.table(file = "results/1000ibd_sequence_counts_before_qc.txt",
                  sep = "\t", header = TRUE)
ctl <- read.table(file = "results/500fg_sequence_counts_before_qc.txt",
                  sep = "\t", header = TRUE)


## Make column to sort reads
ibd$for_rev <- ifelse(grepl(".1$", ibd$Sample), 'forward', 'reverse')
ctl$for_rev <- ifelse(grepl("_1$", ctl$Sample), 'forward', 'reverse')

## Group by only forward reads
ibd <- ibd %>% 
  filter(for_rev=="forward")
ctl <- ctl %>% 
  filter(for_rev=="forward")

## Column for total reads
ibd$Total <- ibd$Unique.Reads + ibd$Duplicate.Reads
ctl$Total <- ctl$Unique.Reads + ctl$Duplicate.Reads

## Mean + SD 
ibd %>% 
  summarise(mean=mean(Total), sd=sd(Total))
ctl %>% 
  summarise(mean=mean(Total), sd=sd(Total))


# Min and max range
ibd %>% 
  summarise(min=min(Total), max=max(Total))
ctl %>% 
  summarise(min=min(Total), max=max(Total))





# AFTER SEQUENCE QUALITY CONTROL ------------------------------------------


# Import data
ibd <- read.table(file = "results/1000ibd_sequence_counts_after_qc.txt",
                  sep = "\t", header = TRUE)
ctl <- read.table(file = "results/500fg_sequence_counts_after_qc.txt",
                  sep = "\t", header = TRUE)
fmt <- read.table(file = "results/fmt_data_sequence_counts.txt",
                  sep = "\t", header = TRUE)


## Make column to sort reads
ibd$for_rev <- ifelse(grepl("_1$", ibd$Sample), 'forward', 'reverse')
ctl$for_rev <- ifelse(grepl("_1$", ctl$Sample), 'forward', 'reverse')

## Group by only forward reads
ibd <- ibd %>% 
  filter(for_rev=="forward")
ctl <- ctl %>% 
  filter(for_rev=="forward")

## FMT samples are single end so don't need sorting

## Column for total reads
ibd$Total <- ibd$Unique.Reads + ibd$Duplicate.Reads
ctl$Total <- ctl$Unique.Reads + ctl$Duplicate.Reads
fmt$Total <- fmt$Unique.Reads + fmt$Duplicate.Reads

## Mean + SD 
ibd %>% 
  summarise(mean=mean(Total), sd=sd(Total))
ctl %>% 
  summarise(mean=mean(Total), sd=sd(Total))
fmt %>% 
  summarise(mean=mean(Total), sd=sd(Total))

# Min and max range
ibd %>% 
  summarise(min=min(Total), max=max(Total))
ctl %>% 
  summarise(min=min(Total), max=max(Total))
fmt %>% 
  summarise(min=min(Total), max=max(Total))
