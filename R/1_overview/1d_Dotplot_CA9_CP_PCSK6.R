# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Dot plot showing the expression levels of CA9, CP, and PCSK6 in each cell type and each sample (non-log space)
# Section:      Results - Single-cell-based ccRCC tumor marker discovery and epigenetic regulation of tumor markers
#=======================================================================================

# Load libraries  and set working directory  -----------------------------------
## load libraries
packages = c(
  "data.table",
  "stringr",
  "plyr",
  "dplyr",
  "ggplot2"
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
## specify genes to plot
genes_plot <- c("CA9", "CP", "PCSK6")
## make colors
colors_cellgroup <- c("#E7298A", "#E69F00", "#56B4E9", "#F0E442", "#D55E00", "#0072B2", "#FB9A99", "#B2DF8A", "#000000",
                               "#1B9E77", "#B15928", "#7570B3", "#90AD1C", "#AA0DFE", "#85660D", "#BDCDFF", "grey80", "grey50")
                               names(colors_cellgroup) <- c("Tumor cells", "CD4+ T-cells", "CD8+ T-cells", "Macrophages", "NK cells", "DC", "Fibroblasts", "Myofibroblasts",  "B-cells",
                                                            "Proximal tubule", "Loop of Henle", "Distal convoluted tubule", 
                                                            'Principle cells', "Intercalated cells", "Podocytes", "Endothelial cells", 
                                                            "Unknown", "Immune others")
                               
# # input data -------------------------------------------------------------------
# ## input the average expression
# exp_wide_df <- fread(data.table = F, input = "../../data/avgexp_sct_data_bycelltype.tsv.gz")
# ## input meta data
# idmetadata_df <- fread(data.table = F, input = "../../data/SampleID.MetaData.tsv")
# ## input barcode-cell type map
# barcode2celltype_df <- fread(input = "../../data/snRNA.UMAPCoordinate.ClusterID.CellType.ByBarcode.tsv.gz", data.table = F)
# 
# # preprocess --------------------------------------------------------------
# ## process cell type labels
# cellgroup_label_df <- data.frame(cell_type13 = c("B-cells", "CD4+ T-cells", "CD8+ T-cells", "DC", "Endothelial cells", "Fibroblasts", "Immune others", "Macrophages", "NK cells", 
#                                                  "Normal epithelial cells", "Tumor cells", "Unknown",
#                                                  "Proximal tubule", "Loop of Henle", "Distal convoluted tubule", 'Principle cells', "Intercalated cells", "Podocytes"))
# cellgroup_label_df <- cellgroup_label_df %>%
#   mutate(cell_type13.columnname = gsub(x = cell_type13, pattern = "\\-|\\+| ", replacement = "."))
# ## count cells
# cellcount_df <- barcode2celltype_df %>%
#   select(orig.ident, Cell_group_w_epithelialcelltypes) %>%
#   table() %>%
#   as.data.frame()
# cellcount_df$cell_group.columnname <- mapvalues(x = cellcount_df$Cell_group_w_epithelialcelltypes, from = cellgroup_label_df$cell_type13, to = as.vector(cellgroup_label_df$cell_type13.columnname))
# cellcount_df <- cellcount_df %>%
#   mutate(id_sample_cell_group = paste0(orig.ident, "_", cell_group.columnname)) %>%
#   mutate(keep = (Freq >= 10))
# ids_sample_cell_group_keep <- cellcount_df$id_sample_cell_group[cellcount_df$keep]
# 
# # make plot data ----------------------------------------------------------
# plotdata_wide_df <- exp_wide_df %>%
#   filter(V1 %in% genes_plot)
# plotdata_df <- melt(data = plotdata_wide_df)
# summary(plotdata_df$value)
# plotdata_df <- plotdata_df %>%
#   mutate(id_sample_cell_group = gsub(x = variable, pattern = "SCT\\.", replacement = "")) %>%
#   mutate(aliquot = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 2)[,1]) %>%
#   mutate(cell_group.columnname = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 2)[,2]) %>%
#   # mutate(x_plot = value)
#   mutate(x_plot = ifelse(V1 == "CP" & value > 10, 10, value)) %>%
#   filter(id_sample_cell_group %in% ids_sample_cell_group_keep)
# 
# plotdata_df$cell_group <- mapvalues(x = plotdata_df$cell_group.columnname, 
#                                     from = cellgroup_label_df$cell_type13.columnname,
#                                     to = as.vector(cellgroup_label_df$cell_type13))
# plotdata_df$id_sample <- mapvalues(x = plotdata_df$aliquot, 
#                                    from = idmetadata_df$Aliquot.snRNA,
#                                    to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
# ## add cell count
# plotdata_df <- plotdata_df %>%
#   filter(!(cell_group %in% c("Unknown", "Immune others", "Normal epithelial cells")))
# ## save plot data
# write.table(x = plotdata_df, file = "../../plot_data/F1d.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# input plot data -------------------------------------------------------
plotdata_df <- fread(data.table = F, input = "../../plot_data/F1d.SourceData.tsv")

# sort plot data ----------------------------------------------------------
ids_samples_sorted <- plotdata_df %>%
  filter(cell_group == "Tumor cells") %>%
  # filter(V1 == "CP") %>%
  filter(V1 == "CA9") %>%
  arrange(desc(x_plot))
ids_samples_sorted <- unique(ids_samples_sorted$id_sample)
ids_samples_sorted <- c(ids_samples_sorted[!grepl(x = ids_samples_sorted, pattern = "\\-N")], "C3L-00079-N", "C3N-00242-N", "C3N-01200-N", "C3L-00088-N")
plotdata_df$y_plot <- factor(x = plotdata_df$id_sample, levels = rev(ids_samples_sorted))
plotdata_df$V1 <- factor(x = plotdata_df$V1, levels = genes_plot)

# plot --------------------------------------------------------------------
p <- ggplot(data = plotdata_df, mapping = aes(x = y_plot, y = x_plot, fill = cell_group, color = as.character(cell_group == "Tumor cells")))
p <- p + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), alpha = 0.7, dotsize = 1.5)
p <- p + scale_fill_manual(values = colors_cellgroup[unique(plotdata_df$cell_group)])
p <- p + guides(fill = guide_legend(override.aes = list(size = 5)))
p <- p + scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA))
p <- p + theme_classic(base_size = 12)
p <- p + coord_flip()
p <- p + facet_grid(cols = vars(V1), scales = "free_x")
p <- p + ylab("Normalized expression")
p <- p + theme(panel.grid.major.y = element_line(size=.1, color="grey50" ), strip.text = element_text(size = 15), strip.background = element_blank())
p <- p + theme(axis.text.y = element_text(size = 10), axis.title.y = element_blank())
p <- p + theme(axis.text.x = element_text(size = 10), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")))

dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out,"F1d.Dotplot.CA9_CP_PCSK6.pdf")
pdf(file2write, width = 11, height = 6.4, useDingbats = F)
print(p)
dev.off()
