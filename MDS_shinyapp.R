library(shiny)
library(RColorBrewer)
library(edgeR)
library(DESeq)
library(limma)

# Working directory
setwd("/Users/hnatri/Dropbox (ASU)/ICGC_HCC/")

# Defining colors
viralPalette <- brewer.pal(8, "Set1")
hbvColor <- viralPalette[1]
hcvColor <- viralPalette[2]
bothColor <- viralPalette[3]
neitherColor <- viralPalette[4]

sexTissuePalette <- brewer.pal(12, "Paired")
maleTumorColor <- sexTissuePalette[4]
maleAdjacentColor <- sexTissuePalette[3]
femaleTumorColor <- sexTissuePalette[6]
femaleAdjacentColor <- sexTissuePalette[5]

voomCounts <- readRDS("voomDGElist.Rdata")
