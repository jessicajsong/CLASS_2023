
---
  
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
options(stringsAsFactors = F)
library(GenomicRanges)
library(tidyverse)
# library(Gviz)
library(IRanges)
```

```

# setting file path to peak files

basepath <- "/scratch/Shares/rinnclass/CLASS_2023/kurt"
path <- "CLASS_2023/group_chip/kurt/results/bwa/mergedLibrary/macs/broadPeak"

# Load in peak files in a tab separated format (tsv)
ATF3_peak1 <- read_tsv(file.path(basepath, "CLASS_2023/group_chip/kurt/results/bwa/mergedLibrary/macs/broadPeak/ATF3_R1_peaks.broadPeak"), col_names=FALSE)
ATF3_peak2 <- read_tsv(file.path(basepath, "CLASS_2023/group_chip/kurt/results/bwa/mergedLibrary/macs/broadPeak/ATF3_R2_peaks.broadPeak"), col_names=FALSE)
ATF3_peak3 <- read_tsv(file.path(basepath, "CLASS_2023/group_chip/kurt/results/bwa/mergedLibrary/macs/broadPeak/ATF3_R1_peaks.broadPeak"), col_names=FALSE)
ATF3_peak4 <- read_tsv(file.path(basepath, "CLASS_2023/group_chip/kurt/results/bwa/mergedLibrary/macs/broadPeak/ATF3_R2_peaks.broadPeak"), col_names=FALSE)

# Rename the column headings
names(ATF3_peak1) <- c('chromosome', 'start', 'end', 'name', 'score', 'strand', 
                       'signalValue', 'pValue', 'qValue')
names(ATF3_peak2) <- c('chromosome', 'start', 'end', 'name', 'score', 'strand', 
                       'signalValue', 'pValue', 'qValue')
names(ATF3_peak3) <- c('chromosome', 'start', 'end', 'name', 'score', 'strand', 
                       'signalValue', 'pValue', 'qValue')
names(ATF3_peak4) <- c('chromosome', 'start', 'end', 'name', 'score', 'strand', 
                       'signalValue', 'pValue', 'qValue')

# To find overlaps
# first we read the peak files in as gRanges object with rtracklayer function.
peaks1 <- rtracklayer::import(file.path(basepath,path, "ATF3_R1_peaks.broadPeak"))
peaks2 <- rtracklayer::import(file.path(basepath,path, "ATF3_R2_peaks.broadPeak"))
peaks3 <- rtracklayer::import(file.path(basepath,path, "ATF3_R3_peaks.broadPeak"))
peaks4 <- rtracklayer::import(file.path(basepath,path, "ATF3_R4_peaks.broadPeak"))

## Overlaps between 1 and 2
ovf <- findOverlaps(peaks1, peaks2)
length(ovf)
# 3848
length(peaks1)
# 3848
length(peaks2)
# 34018

ovf <- findOverlaps(peaks2, peaks1)
length(ovf)
# 3224

## Overlaps between 1 and 3
ovf <- findOverlaps(peaks1, peaks3)
length(ovf)
# 3210
length(peaks1)
# 3848
length(peaks3)
# 49209

ovf <- findOverlaps(peaks3, peaks1)
length(ovf)
# 3210 


## Overlaps between 1 and 4
ovf <- findOverlaps(peaks1, peaks4)
length(ovf)
# 3169
length(peaks1)
# 3848
length(peaks4)
# 60518

ovf <- findOverlaps(peaks4, peaks1)
length(ovf)
# 3169

## Overlaps between 2 and 3
ovf <- findOverlaps(peaks2, peaks3)
length(ovf)
# 16417

ovf <- findOverlaps(peaks3, peaks2)
length(ovf)
# 16417

## Overlaps between 2 and 4
ovf <- findOverlaps(peaks2, peaks4)
length(ovf)
# 17650

ovf <- findOverlaps(peaks4, peaks2)
length(ovf)
# 17650

## Overlaps between 3 and 4
ovf <- findOverlaps(peaks3, peaks4)
length(ovf)
# 33314

ovf <- findOverlaps(peaks4, peaks3)
length(ovf)
# 33314

# Width of peaks
summary(width(peaks1))

# How many "consensus peaks do you have" - what percentage (roughly) of average peak number per file
## (Consensus = peaks in all replicates).
/scratch/Shares/rinnclass/CLASS_2023/kurt/group_chip/kurt/results/bwa/mergedLibrary/macs/broadPeak/consensus/ATF3/ATF3.consensus_peaks.bed
length(ATF3_consensus_peaks)
# 3212 ranges?

# Find how many peaks in the replicate overlap with the consensus peak
con_ovf <-countOverlaps(peaks4, ATF3_consensus_peaks)
sum(con_ovf)
table(con_ovf)

# What does the raw data (IGV) look like for a good consensus peak
##  Need to do on IGV
# What does the raw data look like for a non consensus peak
##  Need to do on IGV

# Load Gencode-v32: for genome features.
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/data/data/genomes/gencode.v32.annotation.gtf")

# How many of ATF3 consensus peaks overlap promoters?
# let's add 1Kb upstream and downstream from the TSS to define "promoters"
gencode_promoters <- promoters(gencode_gr[gencode_gr$type == "gene"], 
                               upstream = 1e3, 
                               downstream = 1e3)
ATF3_consensus_peaks_ov_promoters <- subsetByOverlaps(ATF3_consensus_peaks, gencode_promoters)
length(ATF3_consensus_peaks_ov_promoters)
#1868 ranges

# subset overlaps from peaks1 and promoters
promoter_overlaps_ATF3_peaks1 <- subsetByOverlaps(peaks1, gencode_promoters)
## 1967
promoter_overlaps_ATF3_peaks2 <- subsetByOverlaps(peaks2, gencode_promoters)
## 11659
promoter_overlaps_ATF3_peaks3 <- subsetByOverlaps(peaks3, gencode_promoters)
## 14329
promoter_overlaps_ATF3_peaks4 <- subsetByOverlaps(peaks4, gencode_promoters)
## 15796


```