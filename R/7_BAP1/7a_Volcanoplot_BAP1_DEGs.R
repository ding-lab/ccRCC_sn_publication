# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Volcano plot displaying the differentially expressed genes (DEGs) between the tumor cells of BAP1-mutated tumors (26,806 cells) vs. tumor cells of non-BAP1/PBRM1-mutated tumors (31,002 cells) by snRNA-seq data
# Section:      Results - Impact of BAP1 mutations on the transcriptional network in ccRCC
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

# input -------------------------------------------------------------------
deg_df <- fread(data.table = F, input = "../../data/BAP1_snRNA_DEGs.CNVcorrected.20210913.v1.tsv.gz")

# set plotting parameters -------------------------------------------------
color_red <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[2]
genes_highlight <- c("DLC1", #paste0("ARHGAP", c(24, 28, 32, 42)),
                     "PTPRJ", "CDH16", "CPEB3", "NR6A1", "ZBTB16",
                     "CES3", "PDK4", "SERPINA1", "SLC5A1", "TGFBR3",
                     "RAPGEF5", "MAPK9", "EPHA6", "EFNA5")
x_cutoff <- 2
y_cutoff <- 350
# make data for plotting --------------------------------------------------
plot_data_df <- deg_df %>%
  filter(!is.na(FDR.snRNA.cnvcorrected)) %>%
  mutate(log10FDR = -log10(FDR.snRNA.cnvcorrected)) %>%
  mutate(foldchange_type = ifelse(FDR.snRNA.cnvcorrected < 0.05, ifelse(avg_log2FC.snRNA > 0  & Num_down.snRNA == 0 & Num_sig_up.snRNA >= 5, "consistently higher in BAP1-mutants",
                                                                        ifelse(avg_log2FC.snRNA < 0  & Num_up.snRNA == 0 & Num_sig_down.snRNA >= 5, "consistently lower in BAP1-mutants", "insignificant")), "insignificant")) %>%
  mutate(x_plot = ifelse(avg_log2FC.snRNA < (-x_cutoff), -x_cutoff, ifelse(avg_log2FC.snRNA > x_cutoff, x_cutoff,  avg_log2FC.snRNA))) %>%
  mutate(y_plot = ifelse(log10FDR > y_cutoff, y_cutoff, log10FDR)) %>%
  mutate(label_plot = ifelse(genesymbol_deg %in% genes_highlight, genesymbol_deg, NA)) %>%
  arrange(desc(foldchange_type))
table(plot_data_df$foldchange_type)
# plot all markers, highlight specified genes--------------------------------------------------------------------
## plot
p <- ggplot(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot,  color = foldchange_type, label = label_plot))
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey70")
p <- p + geom_point_rast(alpha = 0.8, shape = 16)
p <- p + scale_color_manual(values = c("consistently higher in BAP1-mutants" = color_red, 
                                       "consistently lower in BAP1-mutants" = color_blue, 
                                       "insignificant" = "grey80"))
p <- p + geom_text_repel(data = subset(plot_data_df, x_plot < 0), 
                         color = "black", alpha = 1, size = 6, fontface = "italic",
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0, box.padding = 0.75,
                         xlim = c(-x_cutoff, -0.3), force = 4,
                         max.overlaps = Inf)
p <- p + geom_text_repel(data = subset(plot_data_df, x_plot > 0), 
                         color = "black", alpha = 1, size = 6, fontface = "italic",
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0, box.padding = 0.75,
                         xlim = c(0.2, x_cutoff), force = 4,
                         max.overlaps = Inf)
p <- p + scale_size_area(max_size = 4)
p <- p + theme_classic()
p <- p + xlab("Log2(fold change) of tumor-cell gene expression")
p <- p + ylab("-Log10FDR")
p <- p + ylim(c(0, y_cutoff)) + xlim(c(-x_cutoff, x_cutoff))
p <- p + guides(color = guide_legend(title = "Fold change type", title.position = "top", title.theme = element_text(size = 14), 
                                     nrow = 4, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 15, color = "black"),
               axis.title = element_text(size = 15),
               legend.position = "bottom", legend.box = "horizontal")

# write output ------------------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out, "F7a_Volcanoplot_BAP1_DEgs", ".pdf")
pdf(file2write, width = 5.5, height = 6.5, useDingbats = F)
print(p)
dev.off()

