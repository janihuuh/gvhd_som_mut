

plotLatentUmap <- function(viz_df, cluster){
  
  viz_df_temp <- data.frame(viz_df, "seurat_cluster" = cluster)
  nClusters   <- unique(cluster) %>% length
  
  ## Visualise
  umap_mean <- data.frame(aggregate(latent_umap_1 ~ seurat_cluster, viz_df_temp, median), latent_umap_2 = aggregate(latent_umap_2 ~ seurat_cluster, viz_df_temp, median)[,2])
  
  
  ## Plot UMAPs with TCRGP predictions highlighted
  ggplot() +
    geom_point(data = viz_df_temp, aes(x = latent_umap_1, y = latent_umap_2, color = seurat_cluster), size = 0.8) +
    
    # stat_ellipse(data = viz_df_temp, geom = "polygon", aes(x = latent_umap_1, y = latent_umap_2, color = seurat_cluster, fill = seurat_cluster), alpha = 0.1, lty = "dotted") +
    ggrepel::geom_label_repel(data = umap_mean, aes(x = latent_umap_1, y = latent_umap_2, color = seurat_cluster, label = seurat_cluster), size = 5, color = "black") +
    
    theme_void() + theme(legend.position = "none") +
    scale_color_manual(values = getPalette(nClusters)) +
    scale_fill_manual(values = getPalette(nClusters)) + labs()
  
  
}



ViolinHeatmap <- function(seurat_object, features.plot){
  
  hm <- data.frame(cluster = seurat_object@ident, t(as.matrix(seurat_object@data[features.plot, ])))
  hm <- melt(hm, id.vars = "cluster")
  hm <- hm %>% mutate(clustermean = paste0(cluster,variable))
  
  hm_means <- aggregate(value ~ cluster + variable, hm, mean) %>% mutate(clustermean = paste0(cluster,variable))
  
  
  hm2 <- merge(hm, hm_means, by = "clustermean")
  
  ggplot(hm2, aes(cluster.x,value.x,fill=value.y)) + geom_violin() + scale_y_log10() + theme_minimal() + labs(x = "", y = "log10(expression)") +
    # coord_capped_cart(bottom='both', left='both') + theme(panel.border=element_blank(), axis.line=element_line()) +
    # scale_fill_manual(values = getPalette(nClusters)[cd8_clusters + 1]) + theme_bw() +
    theme(strip.text.x = element_text(size=8, angle=0, hjust = 0),
          strip.text.y = element_text(size=12, face="italic", angle = 45),
          legend.position = "none",
          strip.background = element_blank(),
          strip.placement = "outside") + scale_fill_viridis_c(option = "A") +
    facet_wrap(~variable.x, ncol = 1, scales = "free_y")
  
  
}

calculateFoldchange_2v1 <- function(seurat_object){
  
  df_temp <- data.frame(cluster = Idents(seurat_object), seurat_object@meta.data) %>%
    tidyr::complete(timepoint, cluster, fill = list(z = 0)) %>%
    
    group_by(timepoint, cluster) %>%
    summarise(n = n()) %>% mutate(prop = n / sum(n))
  
  df_temp1 <- df_temp %>% filter(timepoint == "2013")
  df_temp2 <- df_temp %>% filter(timepoint == "2017")

  
  df_temp <- data.frame(cluster  = df_temp1$cluster,
                        "2013"      = 0,
                        "2017"      = log2(df_temp2$prop / df_temp1$prop)) %>% melt(id = "cluster")
  
  return(df_temp)
  
}



plotBarFoldchange <- function(df){
  
  df <- df %>%
    mutate(direction = "stable") %>%
    mutate(direction = ifelse(value > 1, "increase", direction)) %>%
    mutate(direction = ifelse(value < -1, "decrease", direction)) %>%
    mutate(direction = factor(direction, levels = c("increase", "stable", "decrease")))
  
  limits <- max(abs(df$value))
  
  ggplot(df, aes(x = cluster, value, fill = direction, label = value)) +
    geom_bar(stat = "identity") + coord_flip() +
    scale_fill_manual(values = c("salmon", "lightgrey", "dodgerblue")) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_hline(yintercept = -1, linetype = "dotted") + theme_bw() + 
    labs(y = "log2(fc)", x = "") +
    ylim(-limits, limits)
  
}
