# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  UMAP illustration of the tumor-cell clusters for four tumor samples, colored by the cluster name
# Section:      Results - Transcriptome-based tumor-cell subclusters may represent genomically distinct subclones
#=======================================================================================

# Load libraries  and set working directory  -----------------------------------
## load libraries
packages = c(
  "data.table",
  "stringr",
  "plyr",
  "dplyr",
  "ggplot2",
  "RColorBrewer",
  "ggrastr"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}
## set working directory to current file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# input -------------------------------------------------------------------
## barcode-UMAP info for tumor clusters
barcode2umap_df <- fread(data.table = F, input = "../../data/MetaData_TumorCellOnlyReclustered.20210805.v1.tsv.gz")
## barcode to tumor subcluster assignment
barcode2tumorsubcluster_df <- fread(input = "../../data/Barcode2TumorSubclusterId.20210805.v1.tsv.gz", data.table = F)
## scrublet information
barcode2scrublet_df <- fread(input = "../../data/scrublet.united_outputs.20210729.v1.tsv.gz", data.table = F)

# merge data --------------------------------------------------------------
barcode2umap_df <- merge(x = barcode2umap_df, y = barcode2tumorsubcluster_df, 
                         by.x = c("easy_id", "barcode_tumorcellreclustered"),
                         by.y = c("easy_id", "barcode"),
                         all.x = T)
barcode2umap_df <- barcode2umap_df %>%
  filter(easy_id != "C3L-00359-T1")

# plot by each sample ----------------------------------------------------
## make different output files
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
dir_out_now <- paste0(dir_out, "F4b_UMAP_tumorclusters", "/")
dir.create(dir_out_now)
colors_all <- Polychrome::dark.colors(n = 24)

# for (easy_id_tmp in "C3L-00010-T1") {
for (easy_id_tmp in unique(barcode2umap_df$easy_id)) {
  scrublets_df <- barcode2scrublet_df %>%
    filter(Aliquot_WU == easy_id_tmp) %>%
    filter(predicted_doublet)
  barcodes_doublet <- scrublets_df$Barcode; length(barcodes_doublet)
  
  plot_data_df <- barcode2umap_df %>%
    filter(easy_id == easy_id_tmp) %>%
    filter(!(barcode_tumorcellreclustered %in% barcodes_doublet)) %>%
    mutate(Name_TumorCluster = paste0("C", id_manual_cluster_w0+1))
  cellnumber_percluster_df <- plot_data_df %>%
    select(Name_TumorCluster) %>%
    table() %>%
    as.data.frame() %>%
    rename(Name_TumorCluster = ".")
  plot_data_df <- plot_data_df %>%
    mutate(Name_TumorCluster = ifelse(Name_TumorCluster == "CNA" | Name_TumorCluster %in% cellnumber_percluster_df$Name_TumorCluster[cellnumber_percluster_df$Freq < 50], "Minor cluster (<50 cells)", Name_TumorCluster))
  
  ## make color for each cluster
  names_cluster_tmp <- sort(unique(plot_data_df$Name_TumorCluster))
  length_clusters <- length(names_cluster_tmp)
  if("Minor cluster (<50 cells)" %in% names_cluster_tmp) {
    uniq_cluster_colors <- c(colors_all[1:(length_clusters-1)], "grey40")
  } else {
    uniq_cluster_colors <-colors_all[1:length_clusters]
  }
  names(uniq_cluster_colors) <- names_cluster_tmp
  
  ## make plot
  p <- ggplot()
  p <- p + geom_point_rast(data = plot_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = Name_TumorCluster), shape = 16, alpha = 0.8, size = 1)
  p <- p + scale_color_manual(values = uniq_cluster_colors, na.translate = T)
  p <- p + ggtitle(label = easy_id_tmp, subtitle = "Manually grouped tumor clusters")
  # p <- p + ggtitle(label = easy_id_tmp, subtitle = "Original analysis")
  p <- p + guides(colour = guide_legend(override.aes = list(size=3), title = NULL, nrow = 2, label.theme = element_text(size = 20)))
  p <- p + theme_void()
  p <- p + theme(legend.position = "bottom")
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), title = element_text(size = 25))
  # save output -------------------------------------------------------------
  file2write <- paste0(dir_out_now, easy_id_tmp, ".pdf")
  pdf(file2write, width = 5, height = 5, useDingbats = F)
  print(p)
  dev.off()
}


