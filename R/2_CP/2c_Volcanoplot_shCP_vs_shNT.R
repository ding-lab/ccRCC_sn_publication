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
  "ggrastr",
  "ggrepel"
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
fdr_sig_cutoff <- 0.05
genes_highlight <- c("CP", "COL4A1", "OSMR", "TGM2")

# # make plot data --------------------------------------------------
# ## input data
# deg_df <- fread(data.table = F, input = "../../data/Caki1_CP_vs_Caki1_NT.DEGs.20220517.v1.tsv.gz")
# ## make texts for DEG group
# text_up <- paste0("Up (", length(which(deg_df$FDR < 0.05 & deg_df$logFC > 0)), ")")
# text_down <- paste0("Down (", length(which(deg_df$FDR < 0.05 & deg_df$logFC < 0)), ")")
# ## format data
# plot_data_df <- deg_df %>%
#   dplyr::mutate(log10FDR = -log10(FDR)) %>%
#   dplyr::mutate(foldchange_type = ifelse(FDR < fdr_sig_cutoff, ifelse(logFC > 0, text_up, text_down), "insignificant"))
# plot_data_df <- plot_data_df %>%
#   mutate(x_plot = logFC) %>%
#   mutate(y_plot = ifelse(log10FDR > 350, 350, log10FDR)) %>%
#   arrange(desc(foldchange_type)) %>%
#   select(x_plot, y_plot, foldchange_type, external_gene_name, FDR, logFC)
# ## save plot data
# write.table(x = plot_data_df, file = "../../plot_data/F2c.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# input plot data ---------------------------------------------------------
plot_data_df <- fread(data.table = F, input = "../../plot_data/F2c.SourceData.tsv")

# plot-------------------------------------------------------------------
x_pos_cutoff <- min(plot_data_df$logFC[plot_data_df$FDR < 0.05 & plot_data_df$logFC > 0])
x_neg_cutoff <- max(plot_data_df$logFC[plot_data_df$FDR < 0.05 & plot_data_df$logFC < 0])
## make colors
text_up <- unique(plot_data_df$foldchange_type[plot_data_df$FDR < fdr_sig_cutoff & plot_data_df$logFC > 0])
text_down <- unique(plot_data_df$foldchange_type[plot_data_df$FDR < fdr_sig_cutoff & plot_data_df$logFC < 0])
color_red <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[2]
colors_deg <- c(color_red, color_blue, "grey50")
names(colors_deg) <- c(text_up, text_down, "insignificant")
## plot
p <- ggplot()
p <- p + geom_vline(xintercept = x_pos_cutoff, linetype = 2, color = "grey70")
p <- p + geom_vline(xintercept = x_neg_cutoff, linetype = 2, color = "grey70")
p <- p + geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey70")
p <- p + geom_point_rast(data = subset(plot_data_df, foldchange_type == "insignificant"), 
                         mapping = aes(x = x_plot, y = y_plot, color = foldchange_type), alpha = 0.5, shape = 16, size = 0.5)
p <- p + geom_point_rast(data = subset(plot_data_df, foldchange_type != "insignificant"), 
                         mapping = aes(x = x_plot, y = y_plot,  color = foldchange_type), alpha = 0.5, shape = 16, size = 0.5)
p <- p + scale_color_manual(values = colors_deg)
p <- p + geom_text_repel(data = subset(plot_data_df, x_plot > 0 & (external_gene_name %in% genes_highlight)),
                         mapping = aes(x = x_plot, y = y_plot, label = external_gene_name),
                         color = "black", alpha = 1, size = 4.5, #fontface = "bold",
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0, 
                         xlim = c(0, NA), max.overlaps = Inf)
p <- p + geom_text_repel(data = subset(plot_data_df, (x_plot < 0) & (external_gene_name %in% genes_highlight)),
                         mapping = aes(x = x_plot, y = y_plot, label = external_gene_name),
                         color = "black", alpha = 1, size = 6, #fontface = "bold",
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0, box.padding = 0.5,
                         xlim = c(NA, 0), max.overlaps = Inf, force = 5)
p <- p + theme_classic()
p <- p + xlab("Log2(fold change)")
p <- p + ylab("-Log10FDR")
p <- p + guides(color = guide_legend(title = paste0("DEG (", length(which(plot_data_df$FDR < fdr_sig_cutoff)), ")"), 
                                     title.position = "top", title.theme = element_text(size = 14), 
                                     nrow = 4, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 14, color = "black"),
               axis.title = element_text(size = 14),
               legend.position = c(0.80, 0.25),
               legend.box = "horizontal")

## save plot
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out,"F2c.Volcanoplot.shCP_vs_shNT.pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
print(p)
dev.off()


