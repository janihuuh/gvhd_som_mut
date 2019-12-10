
## Work with seurat 2.3.4
# devtools::install_version(package = 'Seurat', version = package_version('2.3.4'))
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RColorBrewer)

me=system("whoami", intern = T)
setwd(paste0("/Users/", me, "/Dropbox/gvhd_scrnaseq/"))

theme_set(theme_classic())

add_guide   <- guides(colour = guide_legend(override.aes = list(size=10)))

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))

## Run all fun_* codes
for(code in list.files("src_manuscript//R/", "fun", full.names = T, recursive = T)){
  
  message(code)
  source(code)
  
}


## mutate() can stop working with scater, and thus unloading scater might be needed
# detach("package:scran", unload=TRUE)
# detach("package:scater", unload=TRUE)
