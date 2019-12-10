
library(dplyr)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(gridExtra)
library(data.table)

## Set global ggplot2-themes and coloring schemes
theme_set(theme_classic())


add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))

setwd("/Users/hru/Dropbox/citeseq/")

## Run all fun_* codes
for(code in list.files("src/R/scRNAseq/", "fun", full.names = T, recursive = T)){
  
  print(code)
  source(code)
  
}

## mutate() can stop working with scater, and thus unloading scater might be needed
# detach("package:scran", unload=TRUE)
# detach("package:scater", unload=TRUE)
