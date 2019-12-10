
tenx_seurat.data <- Read10X(data.dir = "data/10x/filtered_feature_bc_matrix/")
rownames(x = tenx_seurat.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "", x = rownames(x = tenx_seurat.data[["Antibody Capture"]]))

tenx_seurat <- CreateSeuratObject(counts = tenx_seurat.data[["Gene Expression"]], min.cells = 3, min.features = 200)


# We will define an ADT assay, and store raw counts for it
tenx_seurat[["ADT"]] <- CreateAssayObject(tenx_seurat.data[["Antibody Capture"]][, colnames(x = tenx_seurat)])
tenx_seurat <- tenx_seurat %>% NormalizeData(assay = "ADT", normalization.method = "CLR") %>% ScaleData(assay = "ADT")
adt.data    <- GetAssayData(tenx_seurat, slot = "data", assay = "ADT")
adt.dist    <- dist(t((adt.data)))

# Now, we rerun tSNE using our distance matrix defined only on ADT (protein) levels.
tenx_seurat[["tsne_adt"]] <- RunTSNE(adt.dist, assay = "ADT", reduction.key = "adtTSNE_")
tenx_seurat[["adt_snn"]]  <- FindNeighbors(adt.dist)$snn
tenx_seurat               <- FindClusters(tenx_seurat, resolution = 0.2, graph.name = "adt_snn")

tenx_adt_markers <- FindAllMarkers(tenx_seurat, assay = "ADT", max.cells.per.ident = 100)
write.table(tenx_adt_markers, "results/10x/adt_markers.txt", sep = "\t", quote = F, row.names = F)

tenx_adt_markers %>% group_by(cluster) %>% top_n(3, wt = avg_logFC) %>% filter(avg_logFC > 0.5)
new.cluster.ids <- c("CD14+ Mono", "CD4 CD45RO+ T", "CD8 T", "B", "CD4 T", "NK", "pDC", "CD16+ Mono", "T/Mono doublets")
names(new.cluster.ids) <- levels(tenx_seurat)
tenx_seurat <- RenameIdents(tenx_seurat, new.cluster.ids)

DimPlot(tenx_seurat, reduction = "tsne_adt", label = T, cols = getPalette3(9)) + NoLegend()
ggsave("results/10x/adt_tsne.png", width = 8, height = 6)

saveRDS(tenx_seurat, "results/tenx_seurat.rds")




############################################################################################################################################


####################

scale_factor = 10e3     ## Scaling factor
n_hvgs       = 3000     ## How many HVGs to calculate (minus the clonality and unwanted variation genes)
n_pc_stdev   = 2        ## The minimun stdev value for PC to be accepted

####################

###############################################
## Recluster cd4 based on gene-expression
###############################################

cd4_tenx_seurat <- subset(tenx_seurat, idents = c("CD4 T", "CD4 CD45RO+ T"))

cd4_tenx_seurat <- NormalizeData(cd4_tenx_seurat, normalization.method = "LogNormalize", scale.factor = scale_factor)
cd4_tenx_seurat <- FindVariableFeatures(cd4_tenx_seurat, selection.method = "vst", nfeatures = n_hvgs)

seurat_hvg      <- VariableFeatures(cd4_tenx_seurat)
clonality_genes <- getClonalityGenes(cd4_tenx_seurat)
unwanted_genes  <- getUnwantedGenes(cd4_tenx_seurat)

seurat_hvg      <- seurat_hvg[!seurat_hvg %in% clonality_genes]
seurat_hvg      <- seurat_hvg[!seurat_hvg %in% unwanted_genes]
seurat_hvg      <- seurat_hvg[!seurat_hvg %in% grep("^MOUSE", seurat_hvg, value = T)]

VariableFeatures(cd4_tenx_seurat) <- seurat_hvg

write.table(seurat_hvg, "results/10x/cd4_seurat_hvg.txt", sep = "\t", quote = F, row.names = F)


cd4_tenx_seurat <- ScaleData(cd4_tenx_seurat, features = seurat_hvg)

cd4_tenx_seurat <- RunPCA(cd4_tenx_seurat, features = seurat_hvg, npcs = 50)
nPCs            <- sum(cd4_seurat[["pca"]]@stdev > n_pc_stdev)
write.table(nPCs, "results/10x/cd4_nPCs.txt", sep = "\t", quote = F, row.names = F)


## Run tSNE
cd4_tenx_seurat <- RunTSNE(cd4_tenx_seurat, dims = 1:nPCs, method = "FIt-SNE")

## RunUMAP-function does not work, meanwhile try something hacky-ish
umap_df <- cd4_tenx_seurat[["pca"]]@cell.embeddings[,1:nPCs] %>% umapr::umap() %>% select(UMAP1:UMAP2)
umap_df <- CreateDimReducObject(key = "umap", embeddings = as.matrix(x = umap_df))
cd4_tenx_seurat[['umap']] <- umap_df

DimPlot(cd4_tenx_seurat, reduction = "umap", label = T, cols = c("salmon", "dodgerblue")) + NoLegend()
ggsave("results/10x/umap_cd4.png", width = 8, height = 6)

p <- FeaturePlot(cd4_tenx_seurat, reduction = "umap", features = cytotoxic_markers, cols = c("lightgrey", "salmon"), ncol = 5)
ggsave(plot = p, "results/10x/umap_cd4_cytotxic.png", width = 24, height = 8)

saveRDS(cd4_tenx_seurat, file = "results/10x/cd4_tenx_seurat.rds")

