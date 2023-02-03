# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Volcano plot showing differentially enriched TF motifs between ccRCC (tumor) cells and the combined proximal tubule (PT) cells from the NATs
# Section:      Results - Transcription factors mediating glycolytic genes in ccRCC cells

#=======================================================================================
# Load libraries  and set working directory  -----------------------------------
## load libraries
packages = c(
  "data.table",
  "stringr",
  "plyr",
  "dplyr",
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
## colors
color_red <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[2]

# # make data for plotting --------------------------------------------------
# ## input data
# dam_df <- fread(data.table = F, input = "../../data/MotifScoredifference.Tumor_Normal_comparison.20210715.tsv.gz")
# ## format data
# plot_data_df <- dam_df %>%
#   mutate(log10_pvalue = -log10(FDR)) %>%
#   mutate(log10_pvalue_capped = ifelse(is.infinite(log10_pvalue), 150, log10_pvalue)) %>%
#   group_by(TF_Name) %>%
#   summarise(Num_sig_up = length(which(diff > 0 & FDR < 0.05)), Num_sig_down = length(which(diff < 0 & FDR < 0.05)),
#             Num_up = length(which(diff > 0)), Num_down = length(which(diff < 0)),
#             avg_log10_pvalue = mean(log10_pvalue_capped), avg_diff = mean(diff)) %>%
#   mutate(foldchange_type = ifelse(Num_down == 0 & Num_sig_up >= 12, "consistently higher in ccRCC",
#                                   ifelse(Num_up == 0 & Num_sig_down >= 12, "consistently lower in ccRCC", "insignificant"))) %>%
#   mutate(x_plot = ifelse(avg_diff < -2, -2, ifelse(avg_diff > 2, 2,  avg_diff))) %>%
#   mutate(y_plot = avg_log10_pvalue) %>%
#   arrange(desc(foldchange_type), desc(x_plot)) %>%
#   select(x_plot, y_plot, Num_sig_up, TF_Name, foldchange_type)
# ## save plot data
# write.table(x = plot_data_df, file = "../../plot_data/F3a.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# input plot data ---------------------------------------------------------
plot_data_df <- fread(data.table = F, "../../plot_data/F3a.SourceData.tsv")

# plot --------------------------------------------------------------------
## decide TFs to show
tfnames_show <- plot_data_df$TF_Name[plot_data_df$Num_sig_up == 24]
tfnames_show <- c(tfnames_show, "HNF4G", "HNF4A", "HNF1A", "HNF1B", "RXRB", "RXRG", "PPARA::RXRA")
## label TFS to show
plot_data_df <- plot_data_df %>%
  # mutate(TF_modified = gsub(x = TF_Name, pattern = "\\(var.2\\)", replacement = "*")) %>%
  mutate(text_tf = ifelse(TF_Name %in% tfnames_show, TF_Name, NA))
p <- ggplot(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot, color = foldchange_type, label = text_tf))
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey70")
p <- p + geom_point(alpha = 0.8, shape = 16)
p <- p + scale_color_manual(values = c("consistently higher in ccRCC" = color_red, 
                                       "consistently lower in ccRCC" = color_blue, 
                                       "insignificant" = "grey50"))
p <- p + geom_text_repel(alpha = 1, size = 5, force = 3,
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0, box.padding = 0.5,
                         max.overlaps = Inf)
p <- p + scale_size_area(max_size = 4)
p <- p + theme_classic()
p <- p + xlab("Motif score difference for\nccRCC cells vs. PT cells")
p <- p + ylab("-Log10FDR")
p <- p + theme(axis.text = element_text(size = 14, color = "black"),
               axis.title = element_text(size = 14),
               legend.position = "none")

# write output ------------------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out,"F3a.Volcanoplot.diff_motifs_tumorcells_vs_pt", ".pdf")
pdf(file2write, width = 5.5, height = 5, useDingbats = F)
print(p)
dev.off()





