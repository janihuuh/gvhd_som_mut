
extractName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub("\\_.*", "", str1)
}

extractFileName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub(".*\\/", "", str1)
}

extractClusterNumber <- function(strs){
  
  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][1]
    i <- i + 1
  }
  
  return(p)
  
}

## Plot quality control violins
plot_qc_violin <- function(viz_df, var_to_plot, grouping, min, max){
  
  ## Plot univariate violin plots with filter thresholds
  
  # @ params:
  # viz_df = df that contains qc-analysis results and covariates of interest
  # var_to_plot = char, a column name that contains the variable to plot
  # grouping = char, a column name that contains the x-axis grouping
  # min = num, min value for variable
  # max = num, max value for variable
  
  viz_df_temp <- viz_df %>% dplyr::select(var_to_plot)
  
  label_df_min <- ifelse(viz_df_temp > min, "above", "below") %>% table
  label_df_max <- ifelse(viz_df_temp < max, "above", "below") %>% table
  
  ggplot(data = viz_df, aes_string(x = grouping, y = var_to_plot, fill = grouping)) +
    geom_violin(alpha = 0.5) +
    # geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
    
    geom_hline(yintercept = min, linetype = "dotted") +
    geom_hline(yintercept = max, linetype = "dotted") +
    
    annotate(geom = "text", x = 2.5, y = min, label = paste("Below the line:\n", label_df_min[2]), fontface = "italic") +
    annotate(geom = "text", x = 2.5, y = max, label = paste("Above the line:\n", label_df_max[2]), fontface = "italic") +
    
    labs(x = "", title = var_to_plot) + theme(legend.position = "none")
  
}

## Plot HVGs
plotSeuratHVG <- function(object, top_n = 20){
  
  ## @ params
  # object = seurat object
  # top_n  = numeric, how many HVG:s to name
  
  var.features <- VariableFeatures(object)
  top_genes    <- var.features %>% head(n = top_n)
  
  hvf.info <- HVFInfo(object = object)  %>% add_rownames(var = "gene") %>% as.data.frame()
  hvf.info <- hvf.info %>% dplyr::mutate(var.status = ifelse(test = hvf.info$gene %in% var.features,  yes = "yes", no = "no"))
  
  ggplot() +
    geom_point(data = subset(hvf.info, var.status == "no"), aes(mean, variance.standardized), alpha = 1, color = "grey", size = 0.8) +
    geom_point(data = subset(hvf.info, var.status == "yes"), aes(mean, variance.standardized), alpha = 0.8, color = "dodgerblue", size = 0.8) +
    ggrepel::geom_text_repel(data = subset(hvf.info, gene %in% top_genes), aes(mean, variance.standardized, label = gene), fontface = "italic") +
    
    scale_x_log10() + scale_y_log10() +
    theme_bw() + labs(x = "Average expression", y = "Standardized variance")
  
}

## Change the cluster names


getClusterPhenotypes <- function(clusters){

  # input : vector of clusters 

  clusters <- plyr::revalue(clusters, replace   = c("0"  = "Cytotoxic anti-apoptotic",
                                                    "1"  = "Cytotoxic effector",
                                                    "2"  = "Naive",
                                                    "3"  = "Regulatory",
                                                    "4"  = "Cytotoxic memory",
                                                    "5"  = "Cycling"))

  clusters <- clusters %>% as.character() %>% as.factor() 
  return(clusters)

}


