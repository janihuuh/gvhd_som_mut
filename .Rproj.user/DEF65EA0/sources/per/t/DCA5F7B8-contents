
dir.create("results/mutated_clonotype/", showWarnings = F)

gvhd$pathological_clonotype <- "no"
gvhd$pathological_clonotype[gvhd$new_clonotypes_id %in% c("gvhd_clonotype1", "gvhd_clonotype2")] <- "clonotype"

DimPlot(gvhd, reduction = "umap", group.by = "pathological_clonotype", label = T, repel = T, split.by = "pathological_clonotype", cols = getPalette(4)) + theme(legend.position = "none")
ggsave("results/manuscript/umap_vb20clonotypes.pdf", width = 12, height = )

## Where is the mutated clonotype
table(gvhd$RNA_snn_res.0.4, gvhd$pathological_clonotype) %>% melt %>% group_by(Var2) %>% mutate(prop = value / sum(value)) %>% 
  ggplot(aes(Var1, prop, fill = Var2)) + geom_bar(stat = "identity", position = "dodge")

## Subset to cluster 0 and 1, and analyze clonotype against other

## Compare cluster 0 and 1
cluster_0vs1_markers <- FindMarkers(object = gvhd, ident.1 = 0, ident.2 = 1, only.pos = F, min.pct = 0.01, logfc.threshold = 1e-5, verbose = T, test.use = "t", return.thresh = 1)
cluster_0vs1_markers$direction <- ifelse(cluster_0vs1_markers$avg_logFC > 0, "up", "down")
cluster_0vs1_markers$gene <- rownames(cluster_0vs1_markers)
write.table(cluster_0vs1_markers, file = "results/mutated_clonotype/cluster_markers_clonotype.txt", quote = F, row.names = F, sep = "\t")

rnk             <- data.frame(Name = cluster_0vs1_markers$gene, metric = cluster_0vs1_markers$avg_logFC) %>% arrange(desc(metric))
rnk             <- rnk[!rnk$Name %in% clonality_genes, ]
write.table(rnk, file = "results/mutated_clonotype/cluster_0v1_markers_clonotype_for_gsea.rnk", quote = F, row.names = F, sep = "\t")

ggplot() +
  geom_point(data = subset(cluster_0vs1_markers, p_val_adj > 0.05), aes(avg_logFC, -log(p_val, base = 100)), shape = 21, color = "lightgrey", alpha = 0.8) +
  geom_point(data = subset(cluster_0vs1_markers, p_val_adj < 0.05), aes(avg_logFC, -log(p_val, base = 100), color = direction)) + 
  ggrepel::geom_label_repel(data  = subset(cluster_0vs1_markers, p_val_adj < 1e-100 & !gene %in% clonality_genes & !gene %in% unwanted_genes), aes(avg_logFC, -log(p_val, base = 100), label = gene, color = direction), vjust = 2) +

  scale_color_manual(values = c("darkgrey", "darkred")) +
  labs(x = "avg logFC", y = "-log100(pval)") +
  theme_bw() + theme(legend.position = "none") + xlim(values = c(-1.25, 1.25))
ggsave("results/mutated_clonotype/volcano/clonotype_markers_volcanoplot.pdf", width = 10, height = 8)



## Analyse the cluster 0 and 1; i.e. the cluster with the mutatated clone
gvhd_cytotoxic <- subset(gvhd, idents = c(0,1))

## Add meta.data and add as idents
Idents(gvhd_cytotoxic) <- as.factor(gvhd_cytotoxic$pathological_clonotype)

## Test the differences
cluster_markers           <- FindMarkers(object = gvhd_cytotoxic, ident.1 = "clonotype", ident.2 = "no", only.pos = F, min.pct = 0.001, logfc.threshold = 0.01, do.print = T, test.use = "t", return.thresh = 1)
cluster_markers$direction <- ifelse(cluster_markers$avg_logFC > 0, "up", "down")
cluster_markers$gene      <- rownames(cluster_markers)
write.table(cluster_markers, file = "results/mutated_clonotype/cluster_markers_clonotype.txt", quote = F, row.names = F, sep = "\t")

ggplot() +
  geom_point(data = subset(cluster_markers, p_val_adj > 0.05), aes(avg_logFC, -log(p_val, base = 100)), shape = 21, color = "lightgrey", alpha = 0.8) +
  geom_point(data = subset(cluster_markers, p_val_adj < 0.05), aes(avg_logFC, -log(p_val, base = 100), color = direction)) + 
  ggrepel::geom_label_repel(data  = subset(cluster_markers, p_val_adj < 1e-100 & !gene %in% clonality_genes & !gene %in% unwanted_genes), aes(avg_logFC, -log(p_val, base = 100), label = gene, color = direction), vjust = 2) +

  scale_color_manual(values = c("darkgrey", "darkred")) +
  labs(x = "avg logFC", y = "-log100(pval)") +
  theme_bw() + theme(legend.position = "none") + xlim(values = c(-2, 2))
ggsave("results/mutated_clonotype/volcano/clonotype_markers_volcanoplot.pdf", width = 10, height = 8)


## For GSEA; remove clonality genes
rnk             <- data.frame(Name = cluster_markers$gene, metric = cluster_markers$avg_logFC) %>% arrange(desc(metric))
rnk             <- rnk[!rnk$Name %in% clonality_genes, ]
write.table(rnk, file = "results/mutated_clonotype/cluster_markers_clonotype_for_gsea.rnk", quote = F, row.names = F, sep = "\t")
