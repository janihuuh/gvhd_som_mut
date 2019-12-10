
## Download the citeseq data. Fetched from 10X-website on 3.10.2019


# Load in the RNA UMI matrix
cbmc.rna <- as.sparse(read.csv("data/citeseq/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", sep = ",", header = TRUE, row.names = 1))
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

cbmc.adt <- as.sparse(read.csv("data/citeseq/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",", header = TRUE, row.names = 1))
cbmc.adt <- cbmc.adt[setdiff(rownames(x = cbmc.adt), c("CCR5", "CCR7", "CD10")), ]

cite_seurat <- CreateSeuratObject(counts = cbmc.rna)


# We will define an ADT assay, and store raw counts for it
cite_seurat[["ADT"]] <- CreateAssayObject(counts = cbmc.adt)
cite_seurat <- cite_seurat %>% NormalizeData(assay = "ADT", normalization.method = "CLR") %>% ScaleData(assay = "ADT")
adt.data    <- GetAssayData(cite_seurat, slot = "data", assay = "ADT")
adt.dist    <- dist(t((adt.data)))

# Now, we rerun tSNE using our distance matrix defined only on ADT (protein) levels.
cite_seurat[["tsne_adt"]] <- RunTSNE(adt.dist, assay = "ADT", reduction.key = "adtTSNE_")
cite_seurat[["adt_snn"]]  <- FindNeighbors(adt.dist)$snn
cite_seurat               <- FindClusters(cite_seurat, resolution = 0.2, graph.name = "adt_snn")
new.cluster.ids <- c("CD4 T", "CD14+ Mono", "NK", "B", "CD8 T", "NK", "CD34+", "T/Mono doublets", "CD16+ Mono", "pDCs", "B")
names(new.cluster.ids) <- levels(cite_seurat)
cite_seurat <- RenameIdents(cite_seurat, new.cluster.ids)

saveRDS(cite_seurat, "results/cite_seurat.rds")



####################

scale_factor = 10e3     ## Scaling factor
n_hvgs       = 3000     ## How many HVGs to calculate (minus the clonality and unwanted variation genes)
n_pc_stdev   = 2        ## The minimun stdev value for PC to be accepted

####################

## Select only the CD4+ cells, based on protein expression
cd4_seurat <- subset(cite_seurat, idents = "CD4 T")


###############################################
## Recluster it based on gene-expression
###############################################




## Normalize
cd4_seurat <- NormalizeData(cd4_seurat, normalization.method = "LogNormalize", scale.factor = scale_factor)


## ==== Find HVGs with Seurat
cd4_seurat <- FindVariableFeatures(cd4_seurat, selection.method = "vst", nfeatures = n_hvgs)
seurat_hvg <- VariableFeatures(cd4_seurat)

## Remove clonality and unwanted genes
clonality_genes <- getClonalityGenes(cd4_seurat)
unwanted_genes  <- getUnwantedGenes(cd4_seurat)

seurat_hvg      <- seurat_hvg[!seurat_hvg %in% clonality_genes]
seurat_hvg      <- seurat_hvg[!seurat_hvg %in% unwanted_genes]
seurat_hvg      <- seurat_hvg[!seurat_hvg %in% grep("^MOUSE", seurat_hvg, value = T)]

write.table(seurat_hvg, "results/cd4/seurat_hvg.txt", sep = "\t", quote = F, row.names = F)

VariableFeatures(cd4_seurat) <- seurat_hvg

## Scale data based on the HVGs
cd4_seurat <- ScaleData(cd4_seurat, features = seurat_hvg)

## Perform  PCA
cd4_seurat <- RunPCA(cd4_seurat, features = seurat_hvg, npcs = 50)

## Choose nPC based on stdev
nPCs <- sum(cd4_seurat[["pca"]]@stdev > n_pc_stdev)
message(paste("nPCs:", nPCs, "\n"))
write.table(nPCs, "results/cd4/nPCs.txt", sep = "\t", quote = F, row.names = F)

## Run tSNE
cd4_seurat <- RunTSNE(cd4_seurat, dims = 1:nPCs, method = "FIt-SNE")

## RunUMAP-function does not work, meanwhile try something hacky-ish
umap_df <- cd4_seurat[["pca"]]@cell.embeddings[,1:nPCs] %>% umapr::umap() %>% select(UMAP1:UMAP2)
umap_df <- CreateDimReducObject(key = "umap", embeddings = as.matrix(x = umap_df))
cd4_seurat[['umap']] <- umap_df

saveRDS(cd4_seurat, file = "results/cd4/cd4_seurat.rds")

DimPlot(cd4_seurat, reduction = "pca", dims = c(1,2))
DimPlot(cd4_seurat, reduction = "umap")
DimPlot(cd4_seurat, reduction = "tsne")

FeaturePlot(cd4_seurat, reduction = "tsne", features = c("GZMB", "PRF1"), cols = c("lightgrey", "salmon"))
FeaturePlot(cd4_seurat, reduction = "umap", features = c("GZMB", "PRF1"), cols = c("lightgrey", "salmon"))







# # standard log-normalization
# cite_seurat <- cite_seurat %>% NormalizeData %>% FindVariableFeatures %>% ScaleData %>% RunPCA
# 
# # nPCs <- sum(cite_seurat@reductions$pca@stdev > 2)
# nPCs <- 25
# cite_seurat <- FindNeighbors(cite_seurat, dims = 1:nPCs)
# cite_seurat <- FindClusters(cite_seurat, resolution = 0.8)
# 
# umap_df <- cite_seurat[["pca"]]@cell.embeddings[,1:nPCs] %>% umapr::umap() %>% select(UMAP1:UMAP2)
# umap_df <- CreateDimReducObject(key = "umap", embeddings = as.matrix(x = umap_df))
# cite_seurat[['umap']] <- umap_df
# 
# 
# # Find the markers that define each cluster, and use these to annotate the clusters, we use max.cells.per.ident to speed up the process
# cbmc.rna.markers <- FindAllMarkers(cite_seurat, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)
# 
# 
# # Note, for simplicity we are merging two CD14+ Monocyte clusters (that differ in expression of HLA-DR genes) and NK clusters (that differ in cell cycle stage)
# new.cluster.ids <- c("Memory CD4 T", "CD14+ Mono", "Naive CD4 T", "NK", "CD14+ Mono", "Mouse", "B",  "CD8 T", "CD16+ Mono", "T/Mono doublets", "NK", "CD34+", "Multiplets", "Mouse", "Eryth", "Mk", "Mouse", "DC", "pDCs")
# names(new.cluster.ids) <- levels(cite_seurat)
# cite_seurat <- RenameIdents(cite_seurat, new.cluster.ids)



