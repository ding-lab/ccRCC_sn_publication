# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Dot plot showing the fold changes of expression of tumor-cell markers in ccRCC
# Section:      Results - Single-cell-based ccRCC tumor marker discovery and epigenetic regulation of tumor markers
#=======================================================================================

# Load libraries  and set working directory  -----------------------------------
## load libraries
packages = c(
  "data.table",
  "stringr",
  "plyr",
  "dplyr",
  "ggplot2",
  "RColorBrewer"
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

# specify parameters ------------------------------------------------------
dataname_snrna <- "Tumor cells vs. non-tumor cells (snRNA-seq)"
dataname_snatac <- "Tumor cells vs. non-tumor cells (snATAC-seq)"
dataname_bulk_rna <- "Tumors vs. NATs (bulk RNA-seq)"
dataname_bulk_protein <- "Tumors vs. NATs (bulk proteomics)"
genes_filter <- c("UBE2D2", "COL23A1", "ENPP3", "EGFR", "NDRG1", 
                "PCSK6", "CP", "SHISA9", "SLC6A3", "PLEKHA1", 
                "ABLIM3", "EPHA6", "SEMA6A", "FTO", "KCTD3", 
                "PHKA2", "ABCC3", "PLIN2", "TGFA", "SNAP25", "CA9" )

# # make source data -------------------------------------------------------------------
# genes_process_df <- fread(data.table = F, input = "../../data/ccRCC_markers.Surface.20210824.v1.tsv")
# ## input ATAC fold changes
# geneactivity_fc_df <- fread(data.table = F, input = "../../data/FoldChange_ATACGeneActivity_Tumor_vs_AllOtherCells.20210924.tsv")
# ## make source data 
# genes_process_df <- merge(x = genes_process_df, 
#                           y = geneactivity_fc_df %>%
#                             rename(log2FC.snATAC = avg_log2FC), by = c("Gene"), all.x = T)
# plotdata_wide_df <- genes_process_df %>%
#   filter(Gene %in% genes_filter) %>%
#   dplyr::select(Gene, avg_log2FC.mean.TumorcellsvsNontumor, log2FC.bulkRNA, log2FC.bulkpro, log2FC.snATAC) %>%
#   arrange(desc(avg_log2FC.mean.TumorcellsvsNontumor))
# plotdata_df <- melt(plotdata_wide_df)
# plotdata_df <- plotdata_df %>%
#   mutate(data_type = ifelse(variable == "avg_log2FC.mean.TumorcellsvsNontumor", dataname_snrna,
#                             ifelse(variable == "log2FC.snATAC", dataname_snatac, 
#                                    ifelse(variable == "log2FC.bulkRNA", dataname_bulk_rna, dataname_bulk_protein)))) %>%
#   mutate(foldchange = 2^value) %>%
#   mutate(y_plot = Gene)
# summary(plotdata_df$foldchange)
# plotdata_df <- plotdata_df %>%
#   mutate(x_plot = ifelse(foldchange >= 10, 10, foldchange)) %>%
#   dplyr::select(x_plot, y_plot, Gene, data_type)
# ## save source data
# dir_out <- "../../plot_data/"; dir.create(dir_out)
# write.table(x = plotdata_df, file = "../../plot_data/F1c.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# input source data ------------------------------------------------------
plotdata_df <- fread(data.table = F, input = "../../plot_data/F1c.SourceData.tsv")

# edit plot data ----------------------------------------------------------
plotdata_df$y_plot <- factor(x = plotdata_df$Gene, levels = genes_filter)
plotdata_df$data_type <- factor(x = plotdata_df$data_type, levels = c(dataname_snrna, dataname_snatac, dataname_bulk_rna, dataname_bulk_protein))
## make colors
colors_datatype <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1, 3, 4, 5)]
names(colors_datatype) <- c(dataname_snrna, dataname_bulk_rna, dataname_bulk_protein, dataname_snatac)

# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_dotplot(data = plotdata_df, mapping = aes(x = y_plot, y = x_plot, fill = data_type, color = (data_type == dataname_snrna)),
                      binaxis='y', stackdir='center', position=position_dodge(0.6), alpha = 0.7)
p <- p + scale_fill_manual(values = colors_datatype)
p <- p + scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA))
# p <- p + geom_hline(yintercept = 1, linetype = 2, alpha = 0.5)
p <- p + theme_classic(base_size = 12)
p <- p + coord_flip()
p <- p + scale_y_continuous(breaks = seq(0, 10, 2))
p <- p + ylab("Fold change")
p <- p + theme(panel.grid.major.y = element_line(size=.1, color="black" ))
p <- p + theme(axis.text.y = element_text(size = 12, color = "black"), axis.title.y = element_blank())
p <- p + theme(axis.text.x = element_text(size = 12, color = "black"), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")))
p <- p + theme(legend.position = "top")
p <- p + guides(fill = guide_legend(override.aes = list(size=4), nrow = 4, title = NULL, label.theme = element_text(size = 12)))
## save plot
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out,"F1c.Dotplot.FoldChanges_for_TumorMarkers.pdf")
pdf(file2write, width = 4.25, height = 7, useDingbats = F)
print(p)
dev.off()


