# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Volcano plot showing differentially expressed genes between ccRCC (tumor) cells and the combined proximal tubule (PT) cells from the NATs
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

# input -------------------------------------------------------------------
deg_df <- fread(data.table = F, input = "../../data/Tumor_vs_PT_DEGs.with.CNVcorrection.20210824.v1.tsv.gz")

# set plotting parameters -------------------------------------------------
color_red <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[2]
genes_highlight <- c("HK2", "HK1", #"HK3",  
                     "GPI", "PFKP", #"PFKM", 
                     "ALDOB", #"ALDOA", "ALDOC",
                     "TPI1", #"GAPDH", 
                     "PGK1", #"PGK2", 
                     #"PGAM1", "PGAM2", 
                     "ENO2", "ENO1", #"ENO3",
                     "PKM", #"PKLR",
                     "LDHA", #"LDHB", "LDHC", "LDHD",
                     "PDK1", #"PDK2", "PDK3", "PDK4",
                     #"PGM1", "PGM2", "UGP2", "GYS1", 
                     "PYGL", #"GBE1",
                     # "FASN", "ACACA", "ACLY",
                     "FBP1", "PC", #"PCK1", 
                     "MXI1", "RBPJ" # "HIF1A", "REL", "RELA", "ZNF75D", "HSF2", "NEUROD1", "SREVF2", "NEUROG2", "RREB1", "TBXT"
)

# make data for plotting --------------------------------------------------
plot_data_df <- deg_df %>%
  filter(!is.na(FDR.CNVcorrected)) %>%
  mutate(log10FDR = -log10(FDR.CNVcorrected)) %>%
  mutate(foldchange_type = ifelse(FDR.CNVcorrected < 0.05, ifelse(avg_log2FC.allTumorcellsvsPT > 0  & Num_down == 0 & Num_sig_up >= 15, "consistently higher in ccRCC cells",
                                                                  ifelse(avg_log2FC.allTumorcellsvsPT < 0  & Num_up == 0 & Num_sig_down >= 15, "consistently lower in ccRCC cells", "insignificant")), "insignificant")) %>%
  mutate(x_plot = ifelse(avg_log2FC.allTumorcellsvsPT < -3, -3, ifelse(avg_log2FC.allTumorcellsvsPT > 3, 3,  avg_log2FC.allTumorcellsvsPT))) %>%
  mutate(y_plot = ifelse(log10FDR > 350, 350, log10FDR)) %>%
  mutate(label_plot = ifelse(genesymbol_deg %in% genes_highlight, genesymbol_deg, NA)) %>%
  arrange(desc(foldchange_type))

# plot all markers size not scaled--------------------------------------------------------------------
x_high <- 2.5
x_low <- -2.5
## plot
p <- ggplot(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot,  color = foldchange_type, label = label_plot))
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey70")
p <- p + geom_point(alpha = 0.8, shape = 16)
p <- p + geom_text_repel(data = subset(plot_data_df, x_plot > 0), alpha = 1, size = 4.5, fontface = "italic",
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0, box.padding = 0.5, xlim = c(1, 3), force = 5,
                         max.overlaps = Inf)
p <- p + geom_text_repel(data = subset(plot_data_df, x_plot < 0), alpha = 1, size = 4.5, fontface = "italic",
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0, box.padding = 0.5, xlim = c(-3, -1),
                         max.overlaps = Inf)
p <- p + scale_color_manual(values = c("consistently higher in ccRCC cells" = color_red, 
                                       "consistently lower in ccRCC cells" = color_blue, 
                                       "insignificant" = "grey50"))
p <- p + theme_classic()
p <- p + xlab("Log2(fold change) of sn gene expression for\nccRCC cells vs. PT cells")
p <- p + ylab("-Log10FDR")
p <- p + ylim(c(0, 450)) + xlim(c(-3, 3))
p <- p + guides(color = guide_legend(title = "Fold change type", title.position = "top", title.theme = element_text(size = 14), 
                                     nrow = 4, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 14, color = "black"),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")

# write output ------------------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out,"F3b_Volcanoplot_DEGs_tumorcells_vs_pt", ".pdf")
pdf(file2write, width = 5.5, height = 6.5, useDingbats = F)
print(p)
dev.off()





