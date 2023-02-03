# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Volcano plot showing differentially expressed genes between the EMT tumor clusters and Epi-H tumor clusters highlighted in a. Labels on the right denote known mesenchymal markers, while those on the left denote known markers for PT cells
# Section:      Results - Tumor subgroups with distinct epithelial and mesenchymal features
#=======================================================================================

# Load libraries  and set working directory  -----------------------------------
## load libraries
packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "ggrepel",
  "ggrastr",
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

# set plotting parameters -------------------------------------------------
genes_mesenchymal <- c("SERPINE1", "TGFBI", "VIM", "FN1", "WNT5B", "ITGA5", "JUN", "TWIST1")
genes_epithelal <- c("LRP2", "ABI3BP", "PTGER3", "FRMD3", "SLC28A1", "SLC6A3", "EPB41LA4", "NFIB", "NFIA", "HNF4A", "HNF4G", "CIT")
genes_highlight <- c(genes_mesenchymal, genes_epithelal)
genes_highlight
## set y bottom threshold
y_bottom <- -log10(0.05)
## colors
colors_up_degs <- RColorBrewer::brewer.pal(n = 12, name = "Set1")[1]
colors_down_degs <- RColorBrewer::brewer.pal(n = 6, name = "Set1")[2]

# make data for plotting --------------------------------------------------
# ## input degs
# deg_df <- fread(data.table = F, input = "../../data/Selected_2EMTclusters_vs_5Epithelialclusters.logfc.threshold0.min.pct0.1.min.diff.pct0.AssaySCT.tsv.gz")
# plot_data_df <- deg_df %>%
#   dplyr::rename(genesymbol = genesymbol_deg) %>%
#   mutate(Log10p_val_adj = -log10(x = p_val_adj)) %>%
#   mutate(x_plot = ifelse(avg_log2FC < -3, -3,
#                          ifelse(avg_log2FC > 3, 3, avg_log2FC)))
# ## cap y axis
# y_cap <- max(plot_data_df$Log10p_val_adj[!is.infinite(plot_data_df$Log10p_val_adj)])
# ## set x limits to distinguish colors
# x_pos <- log2(1.2)
# x_neg <- -log2(1.2)
# x_pos <- quantile(x = plot_data_df$avg_log2FC, 0.925)
# x_neg <- quantile(x = plot_data_df$avg_log2FC, 0.025)
# plot_data_df <- plot_data_df %>%
#   mutate(y_plot = ifelse(Log10p_val_adj >= y_cap, y_cap, Log10p_val_adj)) %>%
#   
#   mutate(text_gene = ifelse((y_plot >= y_bottom) & ((genesymbol %in% genes_mesenchymal) & (x_plot >= x_pos)) | ((genesymbol %in% genes_epithelal) & (x_plot <= x_neg)), genesymbol, NA)) %>%
#   mutate(color_plot = ifelse(p_val_adj < 0.05, 
#                              ifelse(x_plot >= 0, "FDR<0.05 (up)", "FDR<0.05 (down)"), "FDR>=0.05")) %>%
#   arrange(factor(color_plot, levels = c("FDR>=0.05", "FDR<0.05 (up)", "FDR<0.05 (down)"))) %>%
#   select(x_plot, y_plot, text_gene, color_plot)
# ## write plot data
# write.table(x = plot_data_df, file = "../../plot_data/F5b.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# input plot data ---------------------------------------------------------
plot_data_df <- fread(data.table = F, input = "../../plot_data/F5b.SourceData.tsv")

# plot--------------------------------------------------------------------
p <- ggplot(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot, label = text_gene, color = color_plot))
p <- p + geom_point_rast(alpha = 0.5, size = 0.5, shape = 16)
p <- p + scale_color_manual(values = c("FDR<0.05 (up)" = colors_up_degs, "FDR<0.05 (down)" = colors_down_degs, "FDR>=0.05" = "grey70"))
p <- p + geom_text_repel(data = subset(plot_data_df, x_plot > 0),
                         color = "black", force = 1, fontface = "italic", segment.alpha = 0.5, segment.size = 0.2, 
                         box.padding = 0.5, min.segment.length = 0,
                         size = 5, max.overlaps = Inf, xlim = c(0.5, 4.5))
p <- p + geom_text_repel(data = subset(plot_data_df,  x_plot < 0),
                         color = "black", force = 4, fontface = "italic", segment.alpha = 0.5, segment.size = 0.2,
                         box.padding = 0.5, min.segment.length = 0,
                         size = 5, max.overlaps = Inf, xlim = c(-3.9, -0.75))
p <- p + theme_classic()
p <- p + ylim(c(0, 350)) + xlim(c(-3.5, 4.5))
p <- p + xlab("Log2(Fold-Change)")
p <- p + ylab("-Log10(P-value-adjusted)")
p <- p + theme(axis.text = element_text(size = 15, color = "black"),
               axis.title = element_text(size = 15), legend.position = "none")

# write output ------------------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out, "F5b.Volcanoplot.EMT_DEGs", ".pdf")
pdf(file2write, width = 6, height = 5, useDingbats = F)
print(p)
dev.off()
