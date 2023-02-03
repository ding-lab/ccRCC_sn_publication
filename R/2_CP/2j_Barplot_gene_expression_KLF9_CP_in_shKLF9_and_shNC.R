# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Barplot showing western blot densitometry of KLF9 and CP proteins in RCC4 cells expressing KLF9 shRNA (sh-KLF9) and RCC4 cells expressing scrambled control (sh-NC)
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
genes_plot <- c("KLF9", "CP")
colors_byline <- RColorBrewer::brewer.pal(n = 6, name = "Set2")[c(1, 2)]
names(colors_byline) <- c("sh-NC", "sh-KLF9")

# # make plot data --------------------------------------------------
# ## input data
# exp_df <- fread(data.table = F, input = "../../data/Knockout_celllines.gene_level.CPM.TMMnormalized.DividedByControl.20220516.v1.tsv.gz")
# ## preprocess
# colnames_value <- colnames(exp_df)[grepl(pattern = "sample", x = colnames(exp_df))]
# exp_df <- as.data.table(exp_df)
# ## make plot data
# plot_data_df <- exp_df %>%
#   filter(external_gene_name %in% genes_plot) %>%
#   melt.data.table(measure.vars = colnames_value) %>%
#   mutate(sample_text = gsub(x = variable, pattern = "sample\\.||_e1", replacement = "")) %>%
#   filter(sample_text %in% samples_plot) %>%
#   dplyr::select(external_gene_name, value, sample_text)
# plot_data_df$sample_text2 <- mapvalues(x = plot_data_df$sample_text, from = samples_plot, to = sampletexts_plot)
# ## save plot data
# write.table(x = plot_data_df, file = "../../plot_data/F2j.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# input plot data ---------------------------------------------------------
plot_data_df <- fread(data.table = F, input = "../../plot_data/F2j.SourceData.tsv")

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
file2write <- paste0(dir_out,"F2j.Barplot.gene_expression_KLF9_CP", ".pdf")
pdf(file2write, width = 4.5, height = 3, useDingbats = F)
print(p)
dev.off()





