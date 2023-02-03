# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Volcano plot showing differentially expressed genes between tumor clusters with top 10% quantile inflammatory scores vs. those with the bottom 10% quantile inflammatory scores (annotated in f)
# Section:      Results - Transcriptome-based tumor-cell subclusters may represent genomically distinct subclones
#=======================================================================================

# Load libraries  and set working directory  -----------------------------------
## load libraries
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "data.table",
  "stringr",
  "ggplot2",
  "ggrepel",
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
p_val_adj_sig_cutoff <- 0.05
genes_highlight <- c("B2M", "C1R", "ENTPD1", "C1S", "C3")

# make data for plotting --------------------------------------------------
## input data
deg_df <- fread(data.table = F, input = "../../data/Inflammatory_score_top_vs_bottom_tumorclusters..logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.tsv.gz")
## make text
text_up <- paste0("Up (", length(which(deg_df$p_val_adj < 0.05 & deg_df$avg_log2FC > 0)), ")")
text_down <- paste0("Down (", length(which(deg_df$p_val_adj < 0.05 & deg_df$avg_log2FC < 0)), ")")
## format
plot_data_df <- deg_df %>%
  mutate(log10p_val_adj = -log10(p_val_adj)) %>%
  mutate(foldchange_type = ifelse(p_val_adj < p_val_adj_sig_cutoff, ifelse(avg_log2FC > 0, text_up, text_down), "insignificant"))
plot_data_df <- plot_data_df %>%
  # mutate(x_plot = ifelse(avg_log2FC < -3, -3, ifelse(avg_log2FC > 3, 3,  avg_log2FC))) %>%
  mutate(x_plot = avg_log2FC) %>%
  mutate(y_plot = ifelse(log10p_val_adj > 350, 350, log10p_val_adj)) %>%
  mutate(label_plot = ifelse(genesymbol_deg %in% genes_highlight,  genesymbol_deg, NA)) %>%
  arrange(desc(foldchange_type)) %>%
  select(x_plot, y_plot, p_val_adj, avg_log2FC, foldchange_type, label_plot)
## write plot data
write.table(x = plot_data_df, file = "../../plot_data/F4h.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# input plot data ---------------------------------------------------------
plot_data_df <- fread(data.table = F, input = "../../plot_data/F4h.SourceData.tsv")

# plot all markers size not scaled--------------------------------------------------------------------
x_pos_cutoff <- min(plot_data_df$avg_log2FC[plot_data_df$p_val_adj < 0.05 & plot_data_df$avg_log2FC > 0])
x_neg_cutoff <- max(plot_data_df$avg_log2FC[plot_data_df$p_val_adj < 0.05 & plot_data_df$avg_log2FC < 0])
## make colors
color_red <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[2]
colors_deg <- c(color_red, color_blue, "grey50")
names(colors_deg) <- c(unique(plot_data_df$foldchange_type[plot_data_df$p_val_adj < 0.05 & plot_data_df$avg_log2FC >0]), 
                       unique(plot_data_df$foldchange_type[plot_data_df$p_val_adj < 0.05 & plot_data_df$avg_log2FC < 0]), 
                       "insignificant")
fontsize_plot <- 16
## plot
p <- ggplot(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot, color = foldchange_type, label = label_plot))
p <- p + geom_vline(xintercept = x_pos_cutoff, linetype = 2, color = "grey70")
p <- p + geom_vline(xintercept = x_neg_cutoff, linetype = 2, color = "grey70")
p <- p + geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey70")
p <- p + geom_point_rast(alpha = 0.5, shape = 16, size = 0.5)
p <- p + scale_color_manual(values = colors_deg)
p <- p + geom_text_repel(min.segment.length = 0, box.padding = 0.5, force = 3,
                         alpha = 1, color = "black", fontface = "italic",
                         segment.size = 0.2, segment.alpha = 0.8, size = 5.5, ylim = c(30,300),
                         max.overlaps = Inf)
p <- p + theme_classic(base_size = fontsize_plot)
# p <- p + xlab("Log2(fold change) of tumor-cell sn gene expression for\nBAP1-mutant tumor vs. other tumors(n=19)")
p <- p + xlab("Log2(fold change)")
p <- p + ylab("-Log10(adjusted P value)")
# p <- p + xlim(c(-3, 3))
p <- p + guides(color = guide_legend(title = paste0("DEG direction"), 
                                     title.position = "top", title.theme = element_text(size = fontsize_plot), 
                                     nrow = 1, override.aes = aes(size = 3), label.theme = element_text(size = fontsize_plot)))
p <- p + theme(axis.text = element_text(size = fontsize_plot, color = "black"),
               axis.title = element_text(size = fontsize_plot),
               legend.position = "bottom",
               # legend.position = c(0.80, 0.20), 
               legend.box.background = element_blank(), legend.background = element_blank(),
               legend.box = "horizontal")

# write output ------------------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out, "F4h.Volcanoplot.inflamm_vs_noninflamm_tumorclusters_degs", ".pdf")
pdf(file2write, width = 4.5, height = 4.25, useDingbats = F)
print(p)
dev.off()
