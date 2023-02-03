# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Bar plot showing the bulk RNA-seq expression of RCC4 cells with or without KLF9 knockdown
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

# specify parameters ------------------------------------------------------
samples_plot <- c("rcc4_scrambled", "rcc4_klf9_c2")
sampletexts_plot <- c("sh-NC", "sh-KLF9")
genes_plot <- c("KLF9", "HK2", "PFKP", "ENO2")
colors_byline <- RColorBrewer::brewer.pal(n = 6, name = "Set2")[c(1, 2)]
names(colors_byline) <- c("sh-NC", "sh-KLF9")

# # format expression data --------------------------------------------------
# ## input data
# exp_df <- fread(data.table = F, input = "../../data/Knockout_celllines.gene_level.CPM.TMMnormalized.DividedByControl.20220516.v1.tsv.gz")
# ## pre-process
# colnames_value <- colnames(exp_df)[grepl(pattern = "sample", x = colnames(exp_df))]
# exp_df <- as.data.table(exp_df)
# ## format
# plot_data_df <- exp_df %>%
#   filter(external_gene_name %in% genes_plot) %>%
#   melt.data.table(measure.vars = colnames_value) %>%
#   mutate(sample_text = gsub(x = variable, pattern = "sample\\.||_e1", replacement = "")) %>%
#   filter(sample_text %in% samples_plot) %>%
#   select(external_gene_name, value, sample_text)
# plot_data_df$sample_text2 <- mapvalues(x = plot_data_df$sample_text, from = samples_plot, to = sampletexts_plot)
# ## write plot data
# write.table(x = plot_data_df, file = "../../plot_data/F3g.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# input plot data ---------------------------------------------------------
plot_data_df <- fread(data.table = F, input = "../../plot_data/F3g.SourceData.tsv")

# plot --------------------------------------------------------------------
plot_data_df$sample_text2 <- factor(x = plot_data_df$sample_text2, levels = sampletexts_plot)
plot_data_df$external_gene_name <- factor(x = plot_data_df$external_gene_name, levels = genes_plot)
p <- ggplot()
p <- p + geom_col(data = plot_data_df, mapping = aes(x = external_gene_name, y = value, fill = sample_text2), position=position_dodge(), color = "black")
p <- p + scale_fill_manual(values = colors_byline[sampletexts_plot])
p <- p + guides(fill = guide_legend(title = "Cell line", title.theme = element_text(size = 15), label.theme = element_text(size = 15)))
p <- p + theme_classic()
p <- p + ylab(label = "% CPM to sh-NC")
p <- p + ggtitle(label = paste0("RNA-seq expression"))
p <- p + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, color = "black", size = 15), axis.text.y = element_text(color = "black", size = 15))
p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(color = "black", size = 15), title = element_text(size = 15))
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out,"F3g.Barplot.gene_expression_KLF9_HK2_PFKP_ENO2", ".pdf")
pdf(file2write, width = 5, height = 3, useDingbats = F)
print(p)
dev.off()
