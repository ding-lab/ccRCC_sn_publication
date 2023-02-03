# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Scatter plot showing the positive correlation of chromatin accessibility and transcriptional changes
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
  "ggpubr",
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

# make plot data ----------------------------------------------------------
# ## input data
# peaks2degs_df <- fread(data.table = F, input = "../../data/BAP1_vs_NonMutant_DAP2DEG.20211011.v1.tsv.gz")
# ## format data
# plotdata_df <- peaks2degs_df %>%
#   filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA)) %>%
#   dplyr::select(avg_log2FC.snATAC, avg_log2FC.snRNA, Gene, peak2gene_type) %>%
#   unique()
# plotdata_df <- plotdata_df %>%
#   mutate(highlight = (Gene %in% c("PTPRJ", "DLC1", "DDIT4", "PEBP1")  | (Gene %in% c("SLC38A1", "RAPGEF5", "EPHA6", "DUSP1", "FABP6") & avg_log2FC.snATAC > 0) | (Gene == "CES3" & avg_log2FC.snATAC < -1.5))) %>%
#   select(avg_log2FC.snATAC, avg_log2FC.snRNA, peak2gene_type, Gene, highlight)
# ## write plot data
# write.table(x = plotdata_df, file = "../../plot_data/F7b.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# input plot data ---------------------------------------------------------
plotdata_df <- fread(data.table = F, input = "../../plot_data/F7b.SourceData.tsv")

# plot highlight genes ----------------------------------------------------
## make colors
colors_peak2genetype <- brewer.pal(n = 7, name = "Dark2")[c(4, 6)]
names(colors_peak2genetype) <- c("Promoter", "Enhancer")
p <- ggscatter(data = plotdata_df, x = "avg_log2FC.snATAC", y = "avg_log2FC.snRNA", color = "peak2gene_type", alpha = 0.8, shape = 16, size = 2.5,
               add = "reg.line", add.params = list(color = "grey50", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
# p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0, size = 4)
p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T), 
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene), 
                         max.overlaps = Inf, size = 6,
                         segment.size = 0.4, segment.alpha = 0.6, min.segment.length = 0, box.padding = 0.5)
p <-  p + scale_color_manual(values = colors_peak2genetype)
p <- p + guides(color = guide_legend(nrow = 1, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 18, color = "black"),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
test_result <- cor.test(plotdata_df$avg_log2FC.snATAC, plotdata_df$avg_log2FC.snRNA)
test_result$p.value

# write output ------------------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out, "F7a.Volcanoplot.BAP1_DEgs", ".pdf")
pdf(file2write, width = 5.5, height = 6, useDingbats = F)
print(p)
dev.off()

