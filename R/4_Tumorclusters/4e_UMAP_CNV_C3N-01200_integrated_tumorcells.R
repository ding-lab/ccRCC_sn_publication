# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  UMAP showing merged data for tumor cells from the above three tumor samples, colored by the original tumor cluster name
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
# ## input the UMAP info
# barcode2umap_df <- fread(data.table = F, input = "../../data/C3N-01200.Tumorcells.Integrated.UMAP_data.20220608.v1.tsv.gz")
# ## input CNV genes
# knowncnvgenes_df <- fread(data.table = F, input = "../../data/Known_CNV_genes.20200528.v1.csv")
# ## set genes to plot
# genes2plot <- knowncnvgenes_df$Gene_Symbol
# map_infercnv_state2category <- function(copy_state) {
#   cnv_cat <- vector(mode = "character", length = length(copy_state))
#   cnv_cat[is.na(copy_state)] <- "Not Available"
#   cnv_cat[copy_state == 1] <- "2 Copies"
#   cnv_cat[copy_state == 0] <- "0 Copies"
#   cnv_cat[copy_state == 0.5] <- "1 Copy"
#   cnv_cat[copy_state == 1.5] <- "3 Copies"
#   cnv_cat[copy_state == 2] <- "4 Copies"
#   cnv_cat[copy_state == 3] <- ">4 Copies"
#   return(cnv_cat)
# }

# pre-process -------------------------------------------------------------
## set infercnv result direcotory
dir_infercnv_output <- "../../data/InferCNV/outputs/"
## set colors
PuBu_colors <- RColorBrewer::brewer.pal(n = 9, name = "PuBu")
PuRd_colors <- RColorBrewer::brewer.pal(n = 9, name = "PuRd")
copy_number_colors <- c("loss" = PuBu_colors[5],
                        "gain" = PuRd_colors[5],
                        "neutral" = "grey50")
## set output directory
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
dir_out_now <- paste0(dir_out, "F4e_UMAP_C3N-01200_integrated_CNV", "/")
dir.create(dir_out_now)

# input CNV data ----------------------------------------------------
# cnv_state_all_df <- NULL
# for (snRNA_aliquot_id_tmp in c("CPT0075120002", "CPT0075130004", "CPT0075140002")) {
#   ## input infercnv CNV state results
#   obs_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, snRNA_aliquot_id_tmp, ".infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt.gz"), data.table = F)
#   ref_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, snRNA_aliquot_id_tmp, ".infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt.gz"), data.table = F)
#   cnv_state_df <- rbind(melt(obs_cnv_state_mat, id.vars = c("V1")), melt(ref_cnv_state_mat, id.vars = c("V1")))
#   cnv_state_df <- cnv_state_df %>%
#     filter(V1 %in% genes2plot) %>%
#     mutate(aliquot = snRNA_aliquot_id_tmp)
#   cnv_state_all_df <- rbind(cnv_state_all_df, cnv_state_df)
# }
# rm(obs_cnv_state_mat)
# rm(ref_cnv_state_mat)

# plot by gene ------------------------------------------------------------
for (gene_tmp in c("VHL", "SQSTM1")) {
  # plot_data_df <- merge(x = barcode2umap_df %>%
  #                         mutate(barcode_individual = str_split_fixed(string = barcode_merged, pattern = "_", n = 2)[,1]), 
  #                       y = cnv_state_all_df %>%
  #                         filter(V1 == gene_tmp), by.x = c("orig.ident", "barcode_individual"), by.y = c("aliquot", "variable"), all.x = T)
  # ## map CNV state value to text
  # plot_data_df$cnv_cat <- map_infercnv_state2category(copy_state = plot_data_df$value)
  # plot_data_df <- plot_data_df %>%
  #   mutate(cnv_cat_simple = ifelse(cnv_cat %in% c("0 Copy", "1 Copy"), "loss",
  #                                  ifelse(cnv_cat %in% c("3 Copies", "4 Copies", ">4 Copies"), "gain", "neutral"))) %>%
  #   select(UMAP_1, UMAP_2, cnv_cat_simple)
  # ## write plot data
  # write.table(x = plot_data_df, file = paste0("../../plot_data/F4e.", gene_tmp, ".SourceData.tsv"), quote = F, sep = "\t", row.names = F)
  
  ## input data
  plot_data_df <- fread(data.table = F, input = paste0("../../plot_data/F4e.", gene_tmp, ".SourceData.tsv"))

  p <- ggplot() +
    geom_point_rast(data = plot_data_df, mapping = aes(UMAP_1, UMAP_2, color=cnv_cat_simple), shape = 16, alpha = 0.8, size = 0.2) +
    scale_color_manual(values = copy_number_colors)
  p <- p + theme_void()
  p <- p + theme(panel.border = element_blank(), 
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "none")
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
  file2write <- paste0(dir_out_now, "C3N-01200.", gene_tmp, ".pdf")
  pdf(file2write, width = 2, height = 2, useDingbats = F)
  print(p)
  dev.off()
}
