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
  "ggplot2"
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

# # make plot data ----------------------------------------------------------
# ## input data
# enricher_out_df <- fread(data.table = F, input = "../../data/ora.msigdb_h_cp.2caki1_cp_vs_2caki1_nt.20220525.v1.tsv.gz")
# plotdata_df <- enricher_out_df %>%
#   filter(Keep) %>%
#   filter(test == "Caki1_cp_vs_nt.down") %>%
#   mutate(size_plot = Count) %>%
#   mutate(x_plot = GeneRatio_num*100) %>%
#   mutate(log10FDR = -log10(p.adjust)) %>%
#   mutate(y_plot = str_split_fixed(string = ID, pattern = "_", n = 2)[,2]) %>%
#   mutate(y_plot = ifelse(y_plot == "EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT", y_plot)) %>%
#   arrange(x_plot) %>%
#   select(x_plot, y_plot, size_plot, log10FDR)
# ## save plot data
# write.table(x = plotdata_df, file = "../../plot_data/F2d.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# input plot data ---------------------------------------------------------
plotdata_df <- fread(data.table = F, input = "../../plot_data/F2d.SourceData.tsv")

# plot enrichment map -----------------------------------------------------
plotdata_df$y_plot <- factor(x = plotdata_df$y_plot, levels = plotdata_df$y_plot)
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, size = size_plot, color = log10FDR))
# p <- p + scale_color_gradientn(colors = c("blue", "purple", "red"))
p <- p + scale_color_gradientn(colors = c("blue", "purple", "red"), 
                               breaks = c(10, 20, 30, 40),
                               guide = guide_colourbar(title = "-log10(p.adjust)", direction = "horizontal", title.position = "top", title.theme = element_text(size = 15),
                                                       label.theme = element_text(size = 15)))
p <- p + scale_size_continuous(breaks = c(60, 80, 100, 120), 
                               guide = guide_legend(direction = "horizontal", title = "Gene count", nrow = 2, byrow = T, title.position = "top", title.theme = element_text(size = 15),
                                                    label.theme = element_text(size = 15)))
p <- p + theme_light(base_size = 12)
p <- p + xlab(label = "Gene ratio (%)")
p <- p + xlim(c(0.0345, 0.082)*100)
p <- p + theme(axis.text = element_text(color = "black", size = 15),
               axis.title.y = element_blank(), axis.title.x = element_text(size = 15),
               legend.position = "right", legend.box = "vertical", legend.background = element_rect(fill = NA, color = NA))

# write output ------------------------------------------------------------
## save plot
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out,"F2d.Dotplot.down_DEGs_shCP_vs_shNT_ORA.pdf")
pdf(file2write, width = 6.25, height = 2, useDingbats = F)
print(p)
dev.off()


