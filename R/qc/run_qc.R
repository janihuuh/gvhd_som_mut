
require(scran)
require(scater)
dir.create("results/qc/", showWarnings = F)


## Remove non-tcr cells as these could be the clonotype
cells.with.tcr <- gvhd@meta.data %>% filter(!is.na(new_clonotypes_id)) %>% pull(barcode)
gvhd           <- subset(gvhd, cells = cells.with.tcr)

### Crate new sce-object
gvhd_counts   <- GetAssayData(gvhd)
gvhd_meta     <- gvhd@meta.data
gvhd_sce      <- SingleCellExperiment(assays = list(counts = gvhd_counts), colData = gvhd_meta)
gvhd_sce      <- calculateQCMetrics(gvhd_sce)

### Do basic filtering first, to not cofound the later analyses
viz_df  <- colData(gvhd_sce) %>% as.data.frame()

###################

min_mito     <- 0
max_mito     <- 15

min_ribo     <- 0
max_ribo     <- Inf

min_features <- 200
max_features <- 4000

min_counts   <- 0
# min_counts   <- 2e3
max_counts   <- Inf

min_pct50    <- 0
max_pct50    <- Inf

###################


#### Violin plots
plot_qc_violin(viz_df, var_to_plot = "total_features_by_counts", grouping = "timepoint", min = min_features, max = max_features)
ggsave("results/qc/total_features.png", width = 6, height = 4)

plot_qc_violin(viz_df, var_to_plot = "total_counts", grouping = "timepoint", min = min_counts, max = max_counts) + scale_y_log10()
ggsave("results/qc/total_counts.png", width = 6, height = 4)

plot_qc_violin(viz_df, var_to_plot = "percent.mt", grouping = "timepoint", min = min_mito, max = max_mito)
ggsave("results/qc/percent_mt.png", width = 6, height = 4)

plot_qc_violin(viz_df, var_to_plot = "percent.ribo", grouping = "timepoint", min = min_ribo, max = max_ribo)
ggsave("results/qc/percent_ribo.png", width = 6, height = 4)

plot_qc_violin(viz_df, var_to_plot = "pct_counts_in_top_50_features", grouping = "timepoint", min = min_pct50, max = max_pct50)
ggsave("results/qc/pct_counts_in_top_50_features.png", width = 6, height = 4)


## Scatter plots
viz_df %>%
  ggplot(aes(total_counts, total_features_by_counts, color = timepoint)) + geom_point(size = 0.3, alpha = 0.5) + scale_x_log10() + scale_y_log10() +
  
  geom_vline(xintercept = min_counts, linetype = "dotted") +
  geom_vline(xintercept = max_counts, linetype = "dotted") +
  
  geom_hline(yintercept = min_features, linetype = "dotted") +
  geom_hline(yintercept = max_features, linetype = "dotted") + add_guide
ggsave("results/qc/counts_vs_genes.png", width = 6, height = 4)

viz_df %>%
  ggplot(aes(pct_counts_in_top_50_features, total_features_by_counts,color = timepoint)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
  scale_y_log10() +
  geom_hline(yintercept = min_features, linetype = "dotted") +
  geom_hline(yintercept = max_features, linetype = "dotted") +
  geom_vline(xintercept = min_pct50, linetype = "dotted") +
  geom_vline(xintercept = max_pct50, linetype = "dotted")
ggsave("results/qc/nGene_vs_nUMIvs_top50features.png", width = 6, height = 4)

viz_df %>%
  ggplot(aes(percent.ribo, total_features_by_counts, color = timepoint)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
  scale_y_log10() +
  geom_hline(yintercept = min_features, linetype = "dotted") +
  geom_hline(yintercept = max_features, linetype = "dotted") +
  geom_vline(xintercept = min_ribo, linetype = "dotted") +
  geom_vline(xintercept = max_ribo, linetype = "dotted")
ggsave("results/qc/ribo_vs_genes.png", width = 6, height = 4)

viz_df %>%
  ggplot(aes(percent.mt, total_features_by_counts, color = timepoint)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
  scale_y_log10() +
  geom_hline(yintercept = min_features, linetype = "dotted") +
  geom_hline(yintercept = max_features, linetype = "dotted") +
  geom_vline(xintercept = max_mito, linetype = "dotted")
ggsave("results/qc/mito_vs_genes.png", width = 6, height = 4)

viz_df %>%
  ggplot(aes(percent.ribo, percent.mt, color = timepoint)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
  geom_hline(yintercept = min_mito, linetype = "dotted") +
  geom_hline(yintercept = max_mito, linetype = "dotted") +
  geom_vline(xintercept = min_ribo, linetype = "dotted") +
  geom_vline(xintercept = max_ribo, linetype = "dotted")
ggsave("results/qc/ribo_vs_mito.png", width = 6, height = 4)

viz_df %>%
  ggplot(aes(percent.ribo, pct_counts_in_top_50_features, color = timepoint)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
  geom_vline(xintercept = min_ribo, linetype = "dotted") +
  geom_vline(xintercept = max_ribo, linetype = "dotted") +
  geom_hline(yintercept = max_pct50, linetype = "dotted") +
  geom_hline(yintercept = min_pct50, linetype = "dotted")
ggsave("results/qc/ribo_vs_top50features.png", width = 6, height = 4)

viz_df %>%
  ggplot(aes(percent.mt, pct_counts_in_top_50_features, color = timepoint)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
  geom_vline(xintercept = max_mito, linetype = "dotted") +
  geom_hline(yintercept = max_pct50, linetype = "dotted") +
  geom_hline(yintercept = min_pct50, linetype = "dotted")
ggsave("results/qc/mito_vs_top50features.png", width = 6, height = 4)

viz_df %>%
  ggplot(aes(log10_total_features_by_counts, pct_counts_in_top_50_features, color = timepoint)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
  geom_hline(yintercept = max_pct50, linetype = "dotted") +
  geom_hline(yintercept = min_pct50, linetype = "dotted")
ggsave("results/qc/nGene_vs_nUMIvs_top50features.png", width = 10, height = 4)



##############################################################################

## In total, we remove with the following conditions:
percent_mito_outlier <- viz_df %>% dplyr::filter(percent.mt > max_mito | percent.mt < min_mito)                                   %>% pull(barcode) %>% as.character()
percent_ribo_outlier <- viz_df %>% dplyr::filter(percent.ribo > max_ribo | percent.ribo < min_ribo)                               %>% pull(barcode) %>% as.character()
features_outlier     <- viz_df %>% dplyr::filter(total_features_by_counts < min_features | total_features_by_counts > max_features)                   %>% pull(barcode) %>% as.character()
umis_outlier         <- viz_df %>% dplyr::filter(total_counts > max_counts | total_counts < min_counts)                           %>% pull(barcode) %>% as.character()
pct_50_outlier       <- viz_df %>% dplyr::filter(pct_counts_in_top_50_features < min_pct50 | pct_counts_in_top_50_features > max_pct50) %>% pull(barcode) %>% as.character()


##############################################################################

outlier_cells        <- c(percent_mito_outlier,
                          percent_ribo_outlier,
                          features_outlier,
                          umis_outlier,
                          pct_50_outlier)

reason               <- c(rep("percent_mito_outlier", length(percent_mito_outlier)),
                          rep("percent_ribo_outlier", length(percent_ribo_outlier)),
                          rep("features_outlier",     length(features_outlier)),
                          rep("umis_outlier",         length(umis_outlier)),
                          rep("pct_50_outlier",       length(pct_50_outlier)))

outlier_df <- data.frame(barcode = outlier_cells, reason = reason) %>% dplyr::mutate(from = substr(barcode, 1, 4))

outlier_df %>% group_by(from,reason) %>% summarise(n = n()) %>%
  ggplot(aes(reorder(reason,n),n,fill=from,label=n)) + geom_bar(stat = "identity") + labs(x = "") + scale_fill_manual(values = brewer.pal(5, "Set1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_bw() + facet_wrap(~from)
ggsave("results/qc/outliers.png", width = 4, height = 3)

write.table(outlier_df, "results/qc/outliers.txt", sep = "\t", quote = F, row.names = F)






## Remove the cells from Seurat-object and save a new seurat-object
cells.to.use <- colnames(gvhd)[!colnames(gvhd) %in% outlier_df$barcode]
gvhd_2   <- subset(gvhd, cells = cells.to.use)

gvhd_meta <- gvhd@meta.data[gvhd@meta.data$barcode %in% cells.to.use, ]
gvhd_meta <- gvhd_meta[match(colnames(gvhd_2), gvhd_meta$barcode), ]

## Make sure that rownames and colnames match and they're in right order
rownames(gvhd_meta) == colnames(gvhd_2)
colnames(gvhd_2)    == gvhd_meta$barcode

gvhd_2@meta.data <- gvhd_meta
gvhd <- gvhd_2
rm(gvhd_2)
saveRDS(gvhd, "results/gvhd_qc.rds")

## mutate() and other dplyr functions can stop working with scater, and thus unloading scater might be needed
detach("package:scran", unload=TRUE)
detach("package:scater", unload=TRUE)

