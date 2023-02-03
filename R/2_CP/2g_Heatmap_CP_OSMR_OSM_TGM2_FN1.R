# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Bar plot showing normalized bulk gene expression of COL4A1, OSMR, and TGM2 in sh-CP-C1, sh-CP-C2, sh-NT1, and sh-NT2
# Section:      Results - CP in mediating tumor extracellular matrix and tumor-stroma interaction in ccRCC
#=======================================================================================

# Load libraries  and set working directory  -----------------------------------
## load libraries
packages = c(
  "data.table",
  "stringr",
  "plyr",
  "dplyr",
  "ggplot2",
  "ComplexHeatmap",
  "circlize"
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

# specify parameters -------------------------------------------------
genes_filter <- c("CP", "OSMR", "OSM", "TGM2", "FN1")

# make plot data ------------------------------
# ## input data
# avgexp_df <- fread(input = "../../data/35_aliquot_merged.avgexp.SCT.data.Cell_group_w_epithelialcelltypes.20210907.v1.tsv.gz", data.table = F)
# ## filtr the rows
# plot_data_df <- avgexp_df %>%
#   rename(gene = V1) %>%
#   filter(gene %in% genes_filter)
# ## filter the columns and make data matrix
# plot_data_raw_mat <- as.matrix(plot_data_df[,-1])
# ## add row names
# rownames(plot_data_raw_mat) <- plot_data_df$gene
# ## filter rows based on the average expression
# genes_plot <- rownames(plot_data_raw_mat)
# genes_plot <- genes_filter[genes_filter %in% genes_plot]
# plot_data_mat <- t(apply(plot_data_raw_mat, 1, scale))
# rownames(plot_data_mat) <- rownames(plot_data_raw_mat); plot_data_mat <- plot_data_mat[genes_plot,]
# colnames(plot_data_mat) <- colnames(plot_data_raw_mat)
# ## filter column
# colnames_plot <- colnames(plot_data_mat)
# celltypes_plot <- gsub(x = colnames_plot, pattern = "SCT\\.", replacement = "")
# colnames_plot <- colnames_plot[!(celltypes_plot %in% c("Unknown", "others", "Immune.others", "Proximal.tubule", "Intercalated.cells", "Loop.of.Henle", "Principle.cells", "Distal.convoluted.tubule", "Podocytes"))]
# plot_data_mat <- plot_data_mat[genes_filter, colnames_plot]
# 
# ## save plot data
# write.table(x = t(plot_data_mat), file = "../../plot_data/F2g.SourceData.tsv", quote = F, sep = "\t", row.names = T)

# input plot data ---------------------------------------------------------
plot_data_df <- fread(data.table = F, input = "../../plot_data/F2g.SourceData.tsv")
plot_data_mat <- as.matrix(plot_data_df); 
rownames(plot_data_mat) <- genes_filter

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
colors_heatmapbody = colorRamp2(c(-1.5, 
                                  0, 
                                  seq(from = 0.5, to = 2, by = 0.5)), 
                                c(color_blue, "white", RColorBrewer::brewer.pal(n = 4, name = "YlOrRd")))

# process column name labels ----------------------------------------------
cellgroup_label_df <- data.frame(cell_type13 = c("B-cells", "CD4+ T-cells", "CD8+ T-cells", "DC", "Endothelial cells", "Fibroblasts", "Immune others", "Macrophages", "NK cells", 
                                                 "Normal epithelial cells", "Tumor cells", "Unknown",
                                                 "Proximal tubule", "Loop of Henle", "Distal convoluted tubule", 'Principle cells', "Intercalated cells", "Podocytes"))
cellgroup_label_df <- cellgroup_label_df %>%
  mutate(cell_type13.columnname = gsub(x = cell_type13, pattern = "\\-|\\+| ", replacement = "."))
celltypes_plot <- gsub(x = colnames(plot_data_mat), pattern = "SCT\\.", replacement = "")
celltypelabels_plot <- mapvalues(x = celltypes_plot, from = cellgroup_label_df$cell_type13.columnname, to = as.vector(cellgroup_label_df$cell_type13))

# plot heatmap body -----------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = t(plot_data_mat), 
                             col = colors_heatmapbody, na_col = color_na, border = "black",
                             show_column_names = T, column_names_side = "top", column_names_gp = gpar(fontface = "italic", fontsize = 13),
                             show_column_dend = F, cluster_columns = F,
                             show_row_names = T, row_names_side = "left",
                             row_names_gp = gpar(fontsize = 13), row_labels = celltypelabels_plot,
                             show_row_dend = F, 
                             column_title = NULL,
                             show_heatmap_legend = F)
## make legend
list_lgd = list(
  Legend(title = "Scaled snRNA\nexpression", title_gp = gpar(fontsize = 13), labels_gp = gpar(fontsize = 13),
         col_fun = colors_heatmapbody, 
         legend_width = unit(2, "cm"),
         direction = "vertical"))

# write output ------------------------------------------------------------
## save plot
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out,"F2g.Heatmap.CP_OSMR_TGM2_FN1.pdf")
pdf(file2write, width = 4, height = 3.2)
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = list_lgd)
dev.off()


