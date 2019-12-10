

## Init output options
output_dir       <- paste0("results/final_clusters/")
dir.create(output_dir, showWarnings = F)

## Init visualisation object
viz_df <- cbind(gvhd@meta.data,
                gvhd[["pca"]]@cell.embeddings[,1:6],
                gvhd[["umap"]]@cell.embeddings,
                cluster = Idents(gvhd)) 
write.table(viz_df, paste0(output_dir, "viz_df.txt"), sep = "\t", quote = F, row.names = F)

umap_mean <- data.frame(aggregate(UMAP_1 ~ cluster, viz_df, median), UMAP_2 = aggregate(UMAP_2 ~ cluster, viz_df, median)[,2])
nClusters <- levels(viz_df$cluster) %>% length

viz_df <- fread("results/final_clusters/viz_df.txt")


### ===== Visualize UMAP QC

DimPlot(gvhd, label = T, repel = T, cols = getPalette3(15)) + theme_void() + theme(legend.position = "none")
ggsave(paste0(output_dir, "UMAP_clusters.png"), width = 8, height = 6)

viz_df %>% ggplot(aes(UMAP_1, UMAP_2, color = as.factor(timepoint))) + geom_point(size = 0.3) + scale_color_manual(values = getPalette(5)[c(2,5)]) + add_guide + labs(color = "time point") + theme_void()
ggsave(paste0(output_dir, "UMAP_timepoint.pdf"), width = 8, height = 6)

viz_df %>% ggplot(aes(UMAP_1, UMAP_2, color = !is.na(new_clonotypes_id))) + geom_point(size = 0.3, alpha = 0.5) + scale_color_manual(values = c("lightgrey", "salmon")) + add_guide + labs(color = "TCRab")
ggsave(paste0(output_dir, "UMAP_tcarb.pdf"), width = 8, height = 6)

viz_df %>% ggplot(aes(UMAP_1, UMAP_2, color = log10(nFeature_RNA))) + geom_point(size = 0.3, alpha = 0.5) + scale_color_viridis_c()
ggsave(paste0(output_dir, "UMAP_nFeature.pdf"), width = 8, height = 6)

viz_df %>% ggplot(aes(UMAP_1, UMAP_2, color = log10(nCount_RNA))) + geom_point(size = 0.3, alpha = 0.5) + scale_color_viridis_c()
ggsave(paste0(output_dir, "UMAP_nCount.pdf"), width = 8, height = 6)

viz_df %>% ggplot(aes(UMAP_1, UMAP_2, color = ifelse(nCount_RNA < 1000, "suspect", "nonsuspect"))) + geom_point(size = 0.3, alpha = 0.5) + labs(color = "nCount_RNA < 1000") + add_guide
ggsave(paste0(output_dir, "UMAP_nCount1000.pdf"), width = 8, height = 6)

viz_df %>% ggplot(aes(UMAP_1, UMAP_2, color = percent.mt)) + geom_point(size = 0.3, alpha = 0.5) + scale_color_viridis_c()
ggsave(paste0(output_dir, "UMAP_percmito.pdf"), width = 8, height = 6)

viz_df %>% ggplot(aes(UMAP_1, UMAP_2, color = percent.ribo)) + geom_point(size = 0.3, alpha = 0.5) + scale_color_viridis_c()
ggsave(paste0(output_dir, "UMAP_percribo.pdf"), width = 8, height = 6)


## Plot individual clusters
p <- NULL
i <- 1
colors = getPalette(nClusters)

for(cluster_temp in levels(viz_df$cluster)){
  
  p[[i]] <- ggplot() +
    geom_point(data = subset(viz_df, cluster != cluster_temp), aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.8) +
    geom_point(data = subset(viz_df, cluster == cluster_temp), aes(x = UMAP_1, y = UMAP_2), color = colors[i], size = 0.8) +
    theme_void() + theme(legend.position = "none")
  
  i <- i + 1
}

pdf(paste0(output_dir, "UMAP_per_cluster.pdf"), width = 8, height = 6)
p
dev.off()




## How much cells with TCRab per cluster?
viz_df %>% group_by(cluster) %>% summarise(has_tcrab = sum(!is.na(new_clonotypes_id)), n = n()) %>% mutate(freq = has_tcrab / n) %>%
  ggplot(aes(cluster, freq, fill = cluster)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = getPalette(nClusters)) +
  ylim(values = c(0,1)) + theme_bw() + theme(legend.position = "none") + labs(y = "cells with TCRa | TCRb | TCRab") + geom_hline(yintercept = 0.5, linetype = "dotted")
ggsave(paste0(output_dir, "bar_cluster_with_tcrab.png"), width = 8, height = 6)


## Visualize time point changes
p <- calculateFoldchange_2v1(gvhd) %>% plotBarFoldchange 
ggsave(plot=p, paste0(output_dir, "bar_timepoint_change_2_v_1.png"), width = 12, height = 6)







## Pie chart the clusters
p <- melt(table(Idents(gvhd))) %>%
  ggplot(aes(x = "", y = value/sum(value), label = round(value/sum(value), 3), fill = as.factor(Var1))) +
  geom_bar(stat = "identity", color = "lightgrey") +
  coord_polar("y", start = 0) +
  # geom_text() +
  scale_fill_manual(values = getPalette(nClusters)) +
  theme_minimal() +
  labs(x = "", y = "", fill = "Clusters")
ggsave(plot = p, paste0(output_dir, "pie_cluster.png"), width = 8, height = 8)


p <- melt(table(Idents(gvhd))) %>%
  ggplot(aes(as.factor(reorder(Var1, value)), value, label = value, fill = as.factor(Var1), label = value)) +
  geom_bar(stat = "identity", color = "lightgrey") + geom_text() +
  scale_fill_manual(values = getPalette(nClusters)) +
  theme_minimal() + theme_bw() + theme(legend.position = "none") + coord_flip() + labs(x = "", y = "nCells")
ggsave(plot = p, paste0(output_dir, "bar_cluster.png"), width = 8, height = 8)



p <- viz_df %>%
  ggplot(aes(cluster, nCount_RNA, fill = cluster)) + geom_violin(alpha = 0.3) + scale_fill_manual(values = getPalette(nClusters)) + theme_bw() + theme(legend.position = "none") + scale_y_log10() + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "violin_nCount_RNA.png"), width = 12, height = 6)

p <- viz_df %>%
  ggplot(aes(cluster, nFeature_RNA, fill = cluster)) + geom_violin(alpha = 0.3) + scale_fill_manual(values = getPalette(nClusters)) + theme_bw() + theme(legend.position = "none") + scale_y_log10() + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "violin_nFeature_RNA.png"), width = 12, height = 6)

p <- viz_df %>%
  ggplot(aes(cluster, percent.mt, fill = cluster)) + geom_violin(alpha = 0.3) + scale_fill_manual(values = getPalette(nClusters)) + theme_bw() + theme(legend.position = "none") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "violin_percent_mt.png"), width = 12, height = 6)

p <- viz_df %>%
  ggplot(aes(cluster, percent.ribo, fill = cluster)) + geom_violin(alpha = 0.3) + scale_fill_manual(values = getPalette(nClusters)) + theme_bw() + theme(legend.position = "none") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "violin_percent.ribo.png"), width = 12, height = 6)

p <- viz_df %>%
  ggplot(aes(cluster, percent.cycle, fill = cluster)) + geom_violin(alpha = 0.3) + scale_fill_manual(values = getPalette(nClusters)) + theme_bw() + theme(legend.position = "none") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "violin_percent.cycle.png"), width = 12, height = 6)


## How patient specific the clusters are
p <- viz_df %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(timepoint), freq, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(nClusters)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "bar_per_timepoint.png"), width = 6, height = 4)

p <- viz_df %>% group_by(cluster, timepoint) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(cluster), freq, fill = as.factor(timepoint))) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(3)) + theme_bw() + labs(x = "", fill = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "bar_per_cluster.png"), width = 12, height = 4)

p <- viz_df %>% mutate(prepost = ifelse(timepoint == "1", "pre", "post")) %>% group_by(cluster, prepost) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(cluster), freq, fill = as.factor(prepost))) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(2)) + theme_bw() + labs(x = "", fill = "") + ggpubr::rotate_x_text(45) +
  geom_hline(yintercept = 0.5, linetype = "dotted")
ggsave(plot = p, paste0(output_dir, "bar_per_prepost.png"), width = 12, height = 4)





######### Markers


## Init output options
output_dir       <- paste0("results/cluster_markers/")
dir.create(output_dir, showWarnings = F)


## Find markers
all_markers     <- FindAllMarkers(gvhd, test.use = "t", verbose = T, random.seed = 123)
all_top_markers <- all_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
fwrite(all_markers, "results/cluster_markers/all_markers.txt", sep = "\t", quote = F, row.names = F)
fwrite(all_top_markers, "results/cluster_markers/all_top_markers.txt", sep = "\t", quote = F, row.names = F)
clusters_top_markers <- all_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

## Dotplots
zhang_cd4_markers_genes <- zhang_cd4_markers %>% unlist %>% unique %>% rev
zheng_cd4_markers_genes <- zheng_cd4_markers %>% unlist %>% unique %>% rev

zhang_cd8_markers_genes <- zhang_cd8_markers %>% unlist %>% unique %>% rev
guo_markers_genes       <- guo_markers %>% unlist %>% unique %>% rev

p <- DotPlot(gvhd, features = zhang_cd4_markers_genes, cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_zhang_cd4.pdf", width = 18, height = 10)

p <- DotPlot(gvhd, features = zheng_cd4_markers_genes, cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_zheng_cd4.pdf", width = 18, height = 10)

p <- DotPlot(gvhd, features = zhang_cd8_markers_genes, cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_zhang_cd8.pdf", width = 18, height = 10)

p <- DotPlot(gvhd, features = guo_markers_genes, cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_guo.pdf", width = 18, height = 10)

p <- DotPlot(gvhd, features = rev(unique(clusters_top_markers$gene)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_de_genes.pdf", width = 18, height = 10)




