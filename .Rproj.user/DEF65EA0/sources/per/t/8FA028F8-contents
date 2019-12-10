


####################

scale_factor = 10e3     ## Scaling factor
n_hvgs       = 3000     ## How many HVGs to calculate (minus the clonality and unwanted variation genes)
n_pc_stdev   = 1.5      ## The minimun stdev value for PC to be accepted
res          = c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)

####################

## Init output options
output_dir       <- paste0("results/hvg/")
dir.create(output_dir, showWarnings = F)

## Normalize
gvhd <- NormalizeData(gvhd, normalization.method = "LogNormalize", scale.factor = scale_factor)

## ==== Find HVGs with Seurat
# gvhd        <- FindVariableFeatures(gvhd, selection.method = "mvp", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
gvhd <- FindVariableFeatures(gvhd, selection.method = "vst", nfeatures = 1e3)

seurat_hvg  <-  VariableFeatures(gvhd)

## Remove clonality and unwanted genes
clonality_genes <- getClonalityGenes(gvhd)
unwanted_genes  <- getUnwantedGenes(gvhd)

seurat_hvg      <- seurat_hvg[!seurat_hvg %in% clonality_genes]
seurat_hvg      <- seurat_hvg[!seurat_hvg %in% unwanted_genes]
write.table(seurat_hvg, paste0(output_dir, "seurat_hvg.txt"), sep = "\t", quote = F, row.names = F)


## Visualize the HVG function
VariableFeatures(gvhd) <- seurat_hvg

p <- plotSeuratHVG(gvhd, top_n = 50)
ggsave(plot = p, paste0(output_dir, "seurat_hvg.png"), width = 8, height = 6)




## Scale data based on all genes 
gvhd <- ScaleData(gvhd, features = rownames(gvhd))

## Perform  PCA based on the HVGs
gvhd <- RunPCA(gvhd, features = seurat_hvg, npcs = 50)

## Choose nPC based on stdev
nPCs <- sum(gvhd[["pca"]]@stdev > n_pc_stdev)
message(paste("nPCs:", nPCs, "\n"))
write.table(nPCs, paste0(output_dir, "nPCs.txt"), sep = "\t", quote = F, row.names = F)

## RunUMAP-function does not work
set.seed(123)
gvhd <- RunUMAP(gvhd, dims = 1:nPCs)

DimPlot(gvhd, reduction = "umap", label = T, group.by = "timepoint")
DimPlot(gvhd, reduction = "umap", label = T)



## Clustering
gvhd$RNA_snn_res.0.3 <- NULL
res       <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
gvhd      <- FindNeighbors(gvhd, reduction = "pca", dims = 1:nPCs)
gvhd      <- FindClusters(object = gvhd, resolution = res)
clustering_columns <- grep("res", colnames(gvhd@meta.data), value = T)



## Plot clustering results
p <- NULL
i <- 1

for(clustering_column in clustering_columns){
  
  message(clustering_column)
  nColors <- gvhd@meta.data[,clustering_column] %>% levels %>% length
  
  p[[i]] <- DimPlot(gvhd, reduction = "umap", group.by = clustering_column, cols = getPalette(nColors), label = T) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(title = clustering_column)
  i <- i + 1
  
}

png("results/final_clusters/umap_per_clustering_column.png", width = 1024, height = 1024)
do.call(grid.arrange, c(p, ncol = 4))
dev.off()


q <- NULL
i <- 1

for(clustering_column in clustering_columns){
  
  message(clustering_column)
  q[[i]] <- gvhd@meta.data[,clustering_column] %>% levels %>% length
  i <- i + 1
  
}

data.frame(res,q) %>% ggplot(aes(res,q,label=q)) + theme_bw() + labs(x = "res parameter", y = "nClusters") + geom_text()
ggsave("results/final_clusters/scatter_res_per_cluster.png", width = 4, height = 3)

DimPlot(gvhd, reduction = "umap", group.by = "RNA_snn_res.0.4", cols = getPalette3(15), label = T)

#### Decide clustering on res.0.4 as it agrees best with the latent umap and plateus
Idents(gvhd) <- gvhd$RNA_snn_res.0.4


### Remove cluster 5 as it is monocyte doublets
idents.to.keep = c(0:4, 6)
gvhd <- subset(gvhd, idents = idents.to.keep)



### Do the other steps 
gvhd            <- FindVariableFeatures(gvhd, selection.method = "vst", nfeatures = 1e3)
seurat_hvg      <- VariableFeatures(gvhd)
clonality_genes <- getClonalityGenes(gvhd)
unwanted_genes  <- getUnwantedGenes(gvhd)
seurat_hvg      <- seurat_hvg[!seurat_hvg %in% clonality_genes]
seurat_hvg      <- seurat_hvg[!seurat_hvg %in% unwanted_genes]
write.table(seurat_hvg, paste0(output_dir, "seurat_hvg.txt"), sep = "\t", quote = F, row.names = F)

VariableFeatures(gvhd) <- seurat_hvg
p <- plotSeuratHVG(gvhd, top_n = 50)
ggsave(plot = p, paste0(output_dir, "seurat_hvg.png"), width = 8, height = 6)

## Scale data based on all genes 
gvhd <- ScaleData(gvhd, features = rownames(gvhd))
gvhd <- RunPCA(gvhd, features = seurat_hvg, npcs = 50)
nPCs <- sum(gvhd[["pca"]]@stdev > n_pc_stdev)
message(paste("nPCs:", nPCs, "\n"))
write.table(nPCs, paste0(output_dir, "nPCs.txt"), sep = "\t", quote = F, row.names = F)
gvhd <- RunUMAP(gvhd, dims = 1:nPCs)

DimPlot(gvhd, reduction = "umap", label = T, group.by = "timepoint")
DimPlot(gvhd, reduction = "umap", label = T)

saveRDS(gvhd, "results/gvhd_tcr_fin.rds")

