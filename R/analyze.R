
global_folder <- here::here()
setwd(global_folder)

## The output folder
output_dir <- "results/"
dir.create(output_dir, showWarnings = F)

## ====== Run the analyses ======

## Initialize

require(Seurat)
require(dplyr)
require(ggplot2)
require(data.table)
require(gridExtra)
require(RColorBrewer)

theme_set(theme_classic())

add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))




## Run all fun_* codes
for(code in list.files("src/R/", "fun", full.names = T, recursive = T)){
  message(code)
  source(code)
}

## Download and preprocess the data
source("src_manuscript/preprocess/run_preprocessRNA.R")
source("src_manuscript/preprocess/run_preprocessTCRab.R")
source("src_manuscript/R/preprocess/run_preprocessRNATCRab.R")

## Do quality control; as this is bead sorted for CD4+ cells, remove B-cell and monocytes
source("src_manuscript/R/qc/run_qc.R.R")
source("src_manuscript/R/qc/run_removeNonTcells.R.R")

## For the QC-passed cells, do clustering and analyze the clustering results
source("src_manuscript/R/clustering/run_hvg.R")
source("src_manuscript/R/clustering/run_plotClusters.R")

## Analyze the mutated clonotype
source("src_manuscript/R/mutated_clonotype/run_mutatedClonotype.R")

## Visualize the results
source("src_manuscript/R/run_plotManuscript.R")

message("Fin")
