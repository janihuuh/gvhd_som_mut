
## For scVI, make csv
dir.create("results/scvi/", showWarnings = F)


## === FHRB1641

## By time points, without HVGs (but remove genes with 0 expression across all cells)
gvhd_2013   <- gvhd@assays$RNA@counts[, gvhd$timepoint == "2013"]
gvhd_2013   <- gvhd_2013[Matrix::rowSums(gvhd_2013) != 0, ]  %>% as.data.frame

gvhd_2017   <- gvhd@assays$RNA@counts[, gvhd$timepoint == "2017"]
gvhd_2017   <- gvhd_2017[Matrix::rowSums(gvhd_2017) != 0, ]  %>% as.data.frame

genes.to.use <- intersect(rownames(gvhd_2013), rownames(gvhd_2017))
length(genes.to.use)

gvhd_2013   <- gvhd_2013[genes.to.use, ]
gvhd_2017   <- gvhd_2017[genes.to.use, ]

rownames(gvhd_2013) == rownames(gvhd_2017)

fwrite(gvhd_2013, file = "results/scvi/gvhd_2013.csv", sep = ",", quote = F, row.names = T)
fwrite(gvhd_2017, file = "results/scvi/gvhd_2017.csv", sep = ",", quote = F, row.names = T)



####### Load results

## Get latent representation; UMAP it
gvhd_latent      <- fread("results/scvi/batch_latent_best.csv", header = F)
gvhd_ident       <- fread("results/scvi/batch_indices_best.csv", header = F)
gvhd_latent_umap <- uwot::umap(gvhd_latent)
gvhd_latent_umap <- gvhd_latent_umap %>% as.data.frame() %>% dplyr::rename(UMAP1 = V1, UMAP2 = V2)
# gvhd_latent_umap <- fread("results/scvi/batch_latent_umap_best.csv") %>% select(UMAP1:UMAP2)

data.frame(gvhd_latent_umap, batch = gvhd_ident$V1) %>%
  ggplot(aes(UMAP1, UMAP2, color = as.factor(batch))) + geom_point(size = 0.3) + add_guide




## Put embeddings into Seurat object
gvhd_latent      <- as.matrix(gvhd_latent)
gvhd_latent_umap <- as.matrix(gvhd_latent_umap)

rownames(gvhd_latent)      <- colnames(gvhd)
rownames(gvhd_latent_umap) <- colnames(gvhd)

rownames(gvhd_latent) == rownames(gvhd_latent_umap)

latent_dim_red            <- CreateDimReducObject(key = "latent", embeddings = as.matrix(x = gvhd_latent))
latent_umap_dim_red       <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = gvhd_latent_umap))

gvhd[['latent']]      <- latent_dim_red
gvhd[['latent_umap']] <- latent_umap_dim_red

data.frame(gvhd[["latent_umap"]]@cell.embeddings, timepoint = gvhd$timepoint) %>%
  ggplot(aes(latent_umap_1 , latent_umap_2, color = as.factor(timepoint))) + geom_point(size = 0.3) + add_guide
ggsave("results/scvi/optimal.pdf", width = 12, height = 10)





## Clustering
res       <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
gvhd      <- FindNeighbors(gvhd, reduction = "latent", dims = 1:30)
gvhd      <- FindClusters(object = gvhd, resolution = res)
clustering_columns <- grep("res", colnames(gvhd@meta.data), value = T)



## Plot clustering results
p <- NULL
i <- 1

for(clustering_column in clustering_columns){
  
  message(clustering_column)
  nColors <- gvhd@meta.data[,clustering_column] %>% levels %>% length
  
  p[[i]] <- DimPlot(gvhd, reduction = "latent_umap", group.by = clustering_column, cols = getPalette(nColors), label = T) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(title = clustering_column)
  i <- i + 1
  
}

png("results/final_clusters/latent_umap_per_clustering_column.png", width = 1024, height = 1024)
do.call(grid.arrange, c(p, ncol = 4))
dev.off()




q <- NULL
i <- 1

for(clustering_column in clustering_columns){
  
  message(clustering_column)
  q[[i]] <- gvhd@meta.data[,clustering_column] %>% levels %>% length
  i <- i + 1
  
}

data.frame(res,q) %>% ggplot(aes((res),q)) + geom_point(shape = 21) + theme_bw()


DimPlot(gvhd, reduction = "latent_umap", group.by = "RNA_snn_res.0.4", cols = getPalette(10), label = T)


#### Decide clustering on res.0.4 as it agrees best with the latent umap and plateu
Idents(gvhd) <- gvhd$RNA_snn_res.0.4

saveRDS(gvhd, "results/gvhd_scvi.rds")
