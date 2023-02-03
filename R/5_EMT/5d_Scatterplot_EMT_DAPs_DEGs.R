# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Scatter plot displaying the Log2 transformed fold change for gene promoter accessibility versus Log2 transformed fold change for gene expression in EMT tumor clusters vs. Epi-H tumor clusters highlighted in a
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

# make plot data ----------------------------------------------------------
# ## input degs
# peaks2degs_df <- fread(data.table = F, input = "../../data/2EMTclusters_vs_5Epithelialclusters.DAP_overlaps_DEGs.20210927.v1.tsv.gz")
# ## format data
# plotdata_df <- peaks2degs_df %>%
#   filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA)) %>%
#   dplyr::select(avg_log2FC.snATAC, avg_log2FC.snRNA, Gene, peak) %>%
#   unique()
# nrow(plotdata_df) ## 1096 gene-peak pairs
# length(unique(plotdata_df$Gene)) ## 1030 genes
# length(unique(plotdata_df$peak)) ## 1096 promoter peaks
# plotdata_df <- plotdata_df %>%
#   mutate(highlight = ((avg_log2FC.snATAC >= 1 & avg_log2FC.snRNA >= 1) | (avg_log2FC.snATAC <= -1 & avg_log2FC.snRNA <= -1) | (Gene %in% c("VIM", "FN1", "CDH2", "WNT5B")))) %>%
#   select(avg_log2FC.snATAC, avg_log2FC.snRNA, highlight, Gene)
# ## write plot data
# write.table(x = plotdata_df, file = "../../plot_data/F5d.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# input plot data ---------------------------------------------------------
plotdata_df <- fread(data.table = F, input = "../../plot_data/F5d.SourceData.tsv")

# highlight genes with fold change > 1 ----------------------------------------------------
p <- ggplot() + geom_point_rast(data = plotdata_df, 
                                mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA),
                                alpha = 0.8, shape = 16, size = 1.5)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + stat_cor(method = "spearman", label.x = 1, label.y = 1.5, size = 4)
p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T & avg_log2FC.snATAC > 0),
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene),
                         max.overlaps = Inf, size = 5,
                         segment.size = 0.4, segment.alpha = 0.6, min.segment.length = 0, box.padding = 0.5, color = RColorBrewer::brewer.pal(n = 9, name = "Set1")[1])
p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T & avg_log2FC.snATAC < 0),
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene),
                         max.overlaps = Inf, size = 5, force = 3,
                         segment.size = 0.4, segment.alpha = 0.6, min.segment.length = 0, box.padding = 0.5, color = RColorBrewer::brewer.pal(n = 6, name = "Set1")[2], xlim = c(-4, -1))
p <- p + guides(color = guide_legend(nrow = 1, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme_classic()
p <- p + theme(axis.text = element_text(size = 18, color = "black"),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
p <- p + ylim(c(-4, 4)) + xlim(c(-4, 4))
p <- p + xlab("Log2(fold change of gene promoter accessiblity)")
p <- p + ylab("Log2(fold change of gene expression)")
test_result <- cor.test(plotdata_df$avg_log2FC.snATAC, plotdata_df$avg_log2FC.snRNA)
test_result$p.value
# write output ------------------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out, "F5d.Scatterplot.EMT_DAPs_overlap_DEGs", ".pdf")
pdf(file2write, width = 5, height = 5, useDingbats = F)
print(p)
dev.off()
