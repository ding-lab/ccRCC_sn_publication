# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  UMAP showing merged data for tumor cells from the above three tumor samples, colored by the original tumor cluster name
# Section:      Results - Transcriptome-based tumor-cell subclusters may represent genomically distinct subclones
#=======================================================================================

# Load libraries  and set working directory  -----------------------------------
## load libraries
packages = c(
  "data.table",
  "stringr",
  "plyr",
  "dplyr",
  "ggplot2",
  "RColorBrewer",
  "ggrastr"
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

# make plot data --------------------------------------------------------------
# ## input the UMAP info
# barcode2umap_df <- fread(data.table = F, input = "../../data/C3N-01200.Tumorcells.Integrated.UMAP_data.20220608.v1.tsv.gz")
# ## input barcode-tumor subcluster info
# barcode2tumorsubcluster_df <- fread(input = "../../data/Barcode2TumorSubclusterId.20210805.v1.tsv.gz", data.table = F)
# ## merge
# barcode2umap_df <- merge(x = barcode2umap_df %>%
#                            mutate(barcode = str_split_fixed(string = barcode_merged, pattern = "_", n = 2)[,1]), 
#                          y = barcode2tumorsubcluster_df, 
#                          by.x = c("orig.ident", "barcode"),
#                          by.y = c("orig.ident", "barcode"),
#                          all.x = T)
# ## make plot data
# plot_data_df <- barcode2umap_df %>%
#   mutate(Text_Cluster = gsub(x = Cluster_Name, pattern = "C3N\\-01200\\-", replacement = "")) %>%
#   select(UMAP_1, UMAP_2, Text_Cluster)
# ## write plot data
# write.table(x = plot_data_df, file = "../../plot_data/F4d.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# input plot data ---------------------------------------------------------
plot_data_df <- fread(data.table = F, input = "../../plot_data/F4d.SourceData.tsv")

# plot----------------------------------------------------
## make colors
texts_cluster_uniq <- sort(unique(plot_data_df$Text_Cluster))
colors_cluster <- RColorBrewer::brewer.pal(n = 8, name = "Accent")[c(1, 2, 5, 8, 6, 7, 3)]
names(colors_cluster) <- texts_cluster_uniq

## make plot
p <- ggplot()
p <- p + geom_point_rast(data = plot_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = Text_Cluster), shape = 16, alpha = 0.8, size = 1)
p <- p + scale_color_manual(values = colors_cluster, na.translate = T)
# p <- p + ggtitle(label = paste0("Tumor-cell subclusters for sample ", easyid_tmp), subtitle = "original")
p <- p + guides(colour = guide_legend(override.aes = list(size=4), title = NULL, ncol = 1))
p <- p + theme_void()
# p <- p + theme(legend.position = "bottom")
p <- p + theme(legend.position = c(0.87, 0.75))
p <- p + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank())
# p
## save plot
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out, "F4d.UMAP.C3N-01200_integrated", ".pdf")
pdf(file2write, width = 4, height = 4, useDingbats = F)
print(p)
dev.off()