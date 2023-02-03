# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Bar plot showing the number of tumor-cell clusters per sample.
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

# make plotting data -------------------------------------------------------------------
# ## input data
# clusters_selected_df <- fread(data.table = F, input = "../../data/intrapatient_tumorclusters.selected.20220706.v1.tsv")
# ## make plottting data
# plot_data_df <- clusters_selected_df %>%
#   group_by(sampleid) %>%
#   summarise(Freq = n())
# ## save plot data
# write.table(x = plot_data_df, file = "../../plot_data/F4a.SourceData.tsv", quote = F, sep = "\t", row.names = F)


# input plot data ---------------------------------------------------------
plot_data_df <- fread(data.table = F, input = "../../plot_data/F4a.SourceData.tsv")

# plot by sample --------------------------------------------------------------------
## make color
p <- ggplot(data = plot_data_df, mapping = aes(x = Freq))
p <- p + geom_histogram(color="black",  binwidth = 1, fill = RColorBrewer::brewer.pal(n = 3, name = "Accent")[1])
p <- p + scale_x_continuous(breaks = 1:10)
p <- p + xlab("Number of Tumor Subclusters")
p <- p + ylab("Number of Samples")
p <- p + theme_bw()
p <- p + theme(axis.text.y = element_text(size = 15, color = "black"),
               axis.text.x = element_text(vjust = 0.1, size = 15, color = "black"),
               axis.title = element_text(size = 15))

# save output -------------------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out,"F4a.Histogram.number_tumorclusters", ".pdf")
pdf(file = file2write, width = 4, height = 3, useDingbats = F)
print(p)
dev.off()
