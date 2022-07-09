# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Bubble plot showing the pathways over-represented in genes up-regulated (top) and down-regulated (bottom) in ccRCC cells compared to the PT cells
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
enricher_out_df <- fread(data.table = F, input = "../../data/ORA_tumorcells_vs_pt.down_degs.tsv.gz")

# select non-overlapping pathways -----------------------------------------
## similarity value < 0.1
pathways_selected <- c("REACTOME_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT", "WP_EICOSANOID_METABOLISM_VIA_LIPO_OXYGENASES_LOX", "NABA_ECM_AFFILIATED", "WP_ZINC_HOMEOSTASIS", "HALLMARK_ESTROGEN_RESPONSE_LATE")

# make plot data ----------------------------------------------------------
plotdata_df <- enricher_out_df
plotdata_df <- plotdata_df %>%
  filter(Description %in% pathways_selected) %>%
  mutate(size_plot = Count) %>%
  mutate(x_plot = (size_plot/768)*100) %>%
  mutate(log10FDR = -log10(p.adjust))
pathway_label_df <- data.frame(Description = pathways_selected,
                               pathway_label = c("SLC transmembrane trasnport", "Eicosanoid metabolism", "Extracellular matrix proteins", "Zinc homeostasis", "Late response to estrogen"))

plotdata_df$y_plot <- mapvalues(x = plotdata_df$Description, from = pathway_label_df$Description, to = as.vector(pathway_label_df$pathway_label))
plotdata_df <- plotdata_df %>%
  arrange(x_plot)
plotdata_df$y_plot <- factor(x = plotdata_df$y_plot, levels = plotdata_df$y_plot)

# plot enrichment map -----------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, size = size_plot, color = log10FDR))
p <- p + scale_color_gradientn(colors = c("blue", "purple", "red"), breaks = c(2, 4, 6), guide = guide_colourbar(direction = "horizontal", title = NULL))
p <- p + scale_size_continuous(limits = c(0, 40), breaks = c(10, 20, 30, 40), guide = guide_legend(direction = "horizontal", title = NULL, nrow = 2, byrow = T))
p <- p + theme_light(base_size = 12)
p <- p + xlab(label = "Gene ratio (%)")
p <- p + theme(axis.text = element_text(color = "black"),
               axis.title.y = element_blank(), axis.title.x = element_text(size = 10),
               legend.position = "right", legend.box = "vertical")

# write output ------------------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out,"F3c_Dotplot_ORA_down_DEGs_tumorcells_vs_pt", ".pdf")
pdf(file2write, width = 4.5, height = 1.4, useDingbats = F)
print(p)
dev.off()





