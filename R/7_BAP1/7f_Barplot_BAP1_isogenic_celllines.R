# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Bar plot showing BAP1 and CES3 gene expression in the BAP1-reconstituted and control SKRC-42 cells
# Section:      Results - Impact of BAP1 mutations on the transcriptional network in ccRCC
#=======================================================================================

# Load libraries  and set working directory  -----------------------------------
## load libraries
packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
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

# preprocess --------------------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)

# plot ------------------------------------------------------
exp_df <- fread(input = "../../data/CPM.TMM_normalized.Manuscript_Cell_Lines.20220710.v1.tsv.gz", data.table = F)
for (gene_plot in c("BAP1", "CES3")) {
  # format expression data --------------------------------------------------
  colnames_id <- colnames(exp_df)[!grepl(x = colnames(exp_df), pattern = "sample")]
  plot_data_df <- exp_df %>%
    filter(external_gene_name %in% gene_plot) %>%
    melt(id.vars = colnames_id) %>%
    # filter(variable %in% c("sample.786_o.sglacz_e1", "sample.786_o.sgbap1_e1", "sample.skrc42.bap1_e1", "sample.skrc42.emptyvector_e1")) %>%
    filter(variable %in% c("sample.skrc42.bap1_e1", "sample.skrc42.emptyvector_e1")) %>%
    mutate(parental_line = str_split_fixed(string = variable, pattern = "\\.", n = 3)[,2]) %>%
    mutate(parental_line_text = ifelse(parental_line == "786_o", "786-O", "SKRC-42")) %>%
    mutate(BAP1_status = ifelse(variable %in% c("sample.786_o.sgbap1_e1", "sample.skrc42.emptyvector_e1"), "BAP1 null", "BAP1 wt"))
  plot_data_df$sample_text <- as.vector(plot_data_df$variable)
  plot_data_df$sample_text[plot_data_df$sample_text == "sample.786_o.sgbap1_e1"] <- "sgBAP1"
  plot_data_df$sample_text[plot_data_df$sample_text == "sample.786_o.sglacz_e1"] <- "sgLacZ"
  plot_data_df$sample_text[plot_data_df$sample_text == "sample.skrc42.bap1_e1"] <- "BAP1"
  plot_data_df$sample_text[plot_data_df$sample_text == "sample.skrc42.emptyvector_e1"] <- "control"
  plot_data_df$sample_text <- factor(x = plot_data_df$sample_text, levels = c("sgLacZ", "sgBAP1", "BAP1", "control"))
  ## write plot data
  write.table(x = plot_data_df, file = paste0("../../plot_data/F7f.", gene_plot, ".SourceData.tsv"), quote = F, sep = "\t", row.names = F)
  
  # make barplot ------------------------------------------------------------
  p <- ggplot()
  p <- p + geom_col(data = plot_data_df, mapping = aes(x = sample_text, y = value, fill = BAP1_status), color = "black")
  p <- p + theme_classic()
  p <- p + ylab(label = paste0(gene_plot, " expression (TPM)"))
  p <- p + theme(strip.background = element_rect(fill = NA, color = NA),
                 panel.spacing = unit(0, "lines"))
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, color = "black"), axis.text.y = element_text(color = "black"))
  p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
  
  # write output ------------------------------------------------------------
  file2write <- paste0(dir_out, "F7f.Barplot.", gene_plot, "_expression", ".pdf")
  pdf(file2write, width = 2.1, height = 2, useDingbats = F)
  print(p)
  dev.off()
}
