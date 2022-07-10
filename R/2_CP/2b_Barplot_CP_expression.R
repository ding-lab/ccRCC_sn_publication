# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Bar plot showing normalized bulk gene expression of CP.
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

# input -------------------------------------------------------------------
exp_df <- fread(input = "../../data/CPM.TMM_normalized.Manuscript_Cell_Lines.20220710.v1.tsv.gz", data.table = F)

# specify parameters ------------------------------------------------------
colnames_value <- colnames(exp_df)[grepl(pattern = "sample", x = colnames(exp_df))]
exp_df <- as.data.table(exp_df)
genes_plot <- c("CP")
samples_plot <- c("caki_1_control_e1", "dr_caki_1_rna", "caki_1_cp_c2_e1", "caki_1_cp_c1_e1")
sampletext_plot <- c("caki1_nt1", "caki1_nt2", "caki1_cp_c2", "caki1_cp_c1")
sampletext_plot <- c("sh-NT1", "sh-NT2", "sh-CP-C2", "sh-CP-C1")

# make colors -------------------------------------------------------------
color_nt <- RColorBrewer::brewer.pal(n = 6, name = "Set3")[6]
color_cp <- RColorBrewer::brewer.pal(n = 6, name = "Set3")[5]
colors_bysample <- c(color_nt, color_nt, color_cp, color_cp)
names(colors_bysample) <- c("sh-NT1", "sh-NT2", "sh-CP-C2", "sh-CP-C1")

# format expression data --------------------------------------------------
plot_data_long_df <- exp_df %>%
  filter(external_gene_name %in% genes_plot) %>%
  melt.data.table(measure.vars = colnames_value) %>%
  mutate(sample = gsub(x = variable, pattern = "sample\\.", replacement = "")) %>%
  filter(sample %in% samples_plot)
plot_data_long_df$sample_text <- mapvalues(x = plot_data_long_df$sample, from = samples_plot, to = sampletext_plot)
plot_data_long_df$sample_text <- factor(x = plot_data_long_df$sample_text, levels = rev(sampletext_plot))


# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_col(data = plot_data_long_df, mapping = aes(x = value, y = sample_text, fill = sample_text), position=position_dodge(), color = "black")
p <- p + scale_fill_manual(values = colors_bysample)
p <- p + scale_x_continuous(breaks = seq(0, 12, 4), labels = seq(0, 12, 4))
p <- p + theme_classic()
p <- p + xlab(label = paste0(genes_plot, " gene CPM"))
p <- p + theme(axis.text = element_text(size = 15, color = "black"))
p <- p + theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 15),
               axis.ticks.y = element_blank(), legend.position = "none")
## save plot
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out,"F2b_Barplot_CP_Expression.pdf")
pdf(file2write, width = 2.5, height = 1.4, useDingbats = F)
print(p)
dev.off()


