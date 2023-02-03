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
  "ggpubr",
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

# # meka plot data -------------------------------------------------------------------
# ## input data
# motifscore_bybarcode_df <- fread(data.table = F, input = "../../data/KLF9.motif_score_bycell.tsv.gz", fill=TRUE)
# plotdata_df <- motifscore_bybarcode_df %>%
#   mutate(cell_group_text = ifelse(cell_group_plot %in% c("Tumor", "EMT tumor cells"), "Tumor cells", "PT cells")) %>%
#   group_by(cell_group_text) %>%
#   mutate(outlier = (motif_score > median(motif_score) + IQR(motif_score) * 1.5) | (motif_score < median(motif_score) - IQR(motif_score) * 1.5)) %>%
#   select(cell_group_text, motif_score, outlier)
# ## save plot data
# write.table(x = plotdata_df, file = "../../plot_data/F2h.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# input plot data ---------------------------------------------------------
plotdata_df <- fread(data.table = F, input = "../../plot_data/F2h.SourceData.tsv")

# plot --------------------------------------------------------------------
color_tumorcell <- RColorBrewer::brewer.pal(n = 9, name = "Dark2")[4]
color_pt <- RColorBrewer::brewer.pal(n = 9, name = "Dark2")[1]
p <- ggviolin(data = plotdata_df, x = "cell_group_text", y = "motif_score", fill = "cell_group_text", color = NA)
p <- p + stat_boxplot(geom = "errorbar", width = 0.2)
p <- p + geom_boxplot(width = 0.2)
p <- p + scale_fill_manual(values = c("Tumor cells" = color_tumorcell, "PT cells" = color_pt))
p <- p + stat_compare_means(method = "t.test", label = "p.format", label.y = 6, label.x = 1.25)
p <- p + ylab(label = paste0("KLF9 ", " motif enrichment"))
p <- p + theme(legend.position = "none", axis.title.x = element_blank(), 
               axis.title.y = element_text(size = 12, color = "black"), axis.text = element_text(color = "black", size = 12))

pvalue_ttest <- t.test(plotdata_df$motif_score[plotdata_df$cell_group_text == "Tumor cells"], plotdata_df$motif_score[plotdata_df$cell_group_text == "PT cells"])
pvalue_ttest$p.value
pvalue_ttest$estimate
pvalue_ttest

# write output ------------------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
## save plot
file2write <- paste0(dir_out,"F2h_Violinplot_KLF9_motif_scores.pdf")
pdf(file2write, width = 2.5, height = 2.5, useDingbats = F)
print(p)
dev.off()
