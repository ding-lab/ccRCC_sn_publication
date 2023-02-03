# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  UMAP showing tumor clusters of three tumor samples from the same patient, colored by the copy number status of VHL and SQSTM1.
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

# specify function and parameters----------------------------------------------
map_infercnv_state2category <- function(copy_state) {
  cnv_cat <- vector(mode = "character", length = length(copy_state))
  cnv_cat[is.na(copy_state)] <- "Not Available"
  cnv_cat[copy_state == 1] <- "2 Copies"
  cnv_cat[copy_state == 0] <- "0 Copies"
  cnv_cat[copy_state == 0.5] <- "1 Copy"
  cnv_cat[copy_state == 1.5] <- "3 Copies"
  cnv_cat[copy_state == 2] <- "4 Copies"
  cnv_cat[copy_state == 3] <- ">4 Copies"
  return(cnv_cat)
}
easy_id_tmp <- "C3N-01200-T1"

# Loop: for each aliquot, input seurat object and infercnv output, plot important genes on UMAP ---------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
dir_out_now <- paste0(dir_out, "F4c_UMAP_CNV_tumorclusters", "/")
dir.create(dir_out_now)
dir_out1 <- paste0(dir_out_now, easy_id_tmp, "/")
dir.create(dir_out1)

for (gene_tmp in c("VHL", "SQSTM1")) {
  file2write <- paste(dir_out1, easy_id_tmp, ".", gene_tmp, ".pdf", sep="")
  ## input plot data
  point_data_df <- fread(data.table = F, input = paste0("../../plot_data/F4c.", gene_tmp, ".", easy_id_tmp, ".tsv"))
  
  ## make text data
  cellnumber_percluster_df <- point_data_df %>%
    select(Text_TumorCluster) %>%
    table() %>%
    as.data.frame() %>%
    rename(Text_TumorCluster = ".")
  text_data_df <- point_data_df %>%
    filter(Text_TumorCluster %in% cellnumber_percluster_df$Text_TumorCluster[cellnumber_percluster_df$Freq >= 50]) %>%
    group_by(Text_TumorCluster) %>%
    summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
  
  p <- ggplot() +
    geom_point_rast(data = point_data_df, mapping = aes(UMAP_1, UMAP_2, color=cnv_cat_simple), alpha = 0.8, size = 0.2) +
    scale_color_manual(values = copy_number_colors)
  p <- p + geom_text_repel(data = text_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, label = Text_TumorCluster))
  p <- p + theme_bw()
  p <- p + theme(panel.border = element_blank(), 
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "none")
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
  pdf(file2write, width = 2, height = 2, useDingbats = F)
  print(p)
  dev.off()
}



