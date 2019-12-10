

####################

scale_factor = 10e3     ## Scaling factor
n_hvgs       = 3000     ## How many HVGs to calculate (minus the clonality and unwanted variation genes)
n_pc_stdev   = 1.5      ## The minimun stdev value for PC to be accepted
res          = c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)

####################


## Normalize
gvhd <- NormalizeData(gvhd, normalization.method = "LogNormalize", scale.factor = scale_factor)

## Find HVGs with Seurat
gvhd <- FindVariableFeatures(gvhd, selection.method = "vst", nfeatures = 1e3)
seurat_hvg  <- VariableFeatures(gvhd)

## Remove clonality and unwanted genes
clonality_genes <- getClonalityGenes(gvhd)
unwanted_genes  <- getUnwantedGenes(gvhd)

seurat_hvg      <- seurat_hvg[!seurat_hvg %in% clonality_genes]
seurat_hvg      <- seurat_hvg[!seurat_hvg %in% unwanted_genes]

## Scale data based on the HVGs
gvhd <- ScaleData(gvhd, features = rownames(gvhd))

## Perform  PCA based on the HVGs
gvhd <- RunPCA(gvhd, features = seurat_hvg, npcs = 50)

## Choose nPC based on stdev
nPCs <- sum(gvhd[["pca"]]@stdev > n_pc_stdev)
message(paste("nPCs:", nPCs, "\n"))

## RunUMAP
gvhd <- RunUMAP(gvhd, dims = 1:nPCs, learning.rate = 1)

## Cluster to find B-cells and monocytes
gvhd <- FindNeighbors(gvhd, reduction = "pca", dims = 1:nPCs)
gvhd <- FindClusters(object = gvhd, resolution = 0.3)

DimPlot(gvhd, reduction = "umap", label = T, group.by = "RNA_snn_res.0.3")

## Remove B-cells and monocytes (the outliers in the UMAP)
Idents(gvhd) <- gvhd$RNA_snn_res.0.3
idents.to.keep <- c(0:8)[!0:8 %in% c(4,8)]
bcell_monocyte <- subset(gvhd, idents = c(4,8))
fwrite(as.data.frame(colnames(bcell_monocyte)), "results/qc/bcell_monocytes_idents.txt", sep = "\t", quote = F, row.names = F)

gvhd <- subset(gvhd, idents = idents.to.keep)

