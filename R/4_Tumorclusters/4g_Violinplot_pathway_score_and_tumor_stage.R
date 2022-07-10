# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Violin plot showing maximum pathway score per tumor sample, grouped by tumor stage.
# Section:      Results - Transcriptome-based tumor-cell subclusters may represent genomically distinct subclones
#=======================================================================================

# Load libraries  and set working directory  -----------------------------------
## load libraries
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "data.table",
  "stringr",
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

# input -------------------------------------------------------------------
## input pathway scores
scores_df <- fread(data.table = F, input = "../../data/Pathwayscores_byTumorcluster.MSigDB.Hallmark.tsv.gz")
## input the pathways to plot
genesets_plot_df <- fread(data.table = F, input = "../../data/Genesets_enriched_in_intrapatienttumorclusterDEGs.20220606.v1.tsv.gz")
## input clinical data
clinical_df <- fread(data.table = F, input = "../../data/Patient_Clinical_Info.20201125.v1.tsv")

# preprocess -------------------------------------------------
pos <- position_jitter(width = 0.2, seed = 1)
clinical_filtered_df <- clinical_df %>% 
  filter(Case != "C3L-00359") %>%
  select(Case, Tumor_Stage_Pathological) 
scoregroups_test <- paste0(gsub(x = genesets_plot_df$Description, pattern = "HALLMARK_", replacement = ""), "_Score")
## set output directory
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
dir_out_now <- paste0(dir_out, "F4g_Violinplot_pathway_score_and_tumor_stage/"); dir.create(dir_out_now)

# plot top score per tumor, compare high grade vs. low grade ------------------------------------------------------
wilcox_pvalues_vec <- NULL
ttest_pvalue_vec <- NULL
# for (scoregroup_tmp in "MTORC1_SIGNALING_Score") {
for (scoregroup_tmp in scoregroups_test) {
  ## make plot data
  scores_tmp_df <- scores_df[, c("cluster_name", scoregroup_tmp)]
  colnames(scores_tmp_df) <- c("cluster_name", "score.bycluster")
  scores_bycase_df <- scores_tmp_df %>%
    mutate(sample_id = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
    mutate(Aliquot_WU = gsub(x = sample_id, pattern = "\\.", replacement = "-")) %>%
    mutate(Case = str_split_fixed(string = Aliquot_WU, pattern = "\\-T|\\-N", n = 2)[,1]) %>%
    group_by(Case) %>%
    summarise(score.bysample = max(score.bycluster))
  plot_data_df <- merge(x = clinical_filtered_df,
                        y = scores_bycase_df,
                        by.x = c("Case"), by.y = c("Case"), all.x = T)
  plot_data_df <- plot_data_df %>%
    mutate(x_plot = ifelse(Tumor_Stage_Pathological %in% c("Stage I", "Stage II"), "Stage I/II", "Stage III/IV")) %>%
    mutate(y_plot = score.bysample)
  # plot_data_df <- plot_data_df[!duplicated(plot_data_df$Case),]
  
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  p = p + geom_boxplot(width=.1)
  p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 2)
  p = p + stat_compare_means(data = plot_data_df, 
                             mapping = aes(x = x_plot, y = y_plot), 
                             hide.ns = F, method = "wilcox.test", label = "p.format")
  # p <- p + ggtitle(label = paste0("Max ", scoregroup_tmp, " by tumor"))
  p <- p + ylab(scoregroup_tmp)
  p <- p + theme_classic()
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.title.x = element_blank(), axis.text = element_text(color = "black", size = 12))
  file2write <- paste0(dir_out_now, scoregroup_tmp, ".pdf")
  pdf(file2write, width = 2.5, height = 2, useDingbats = F)
  print(p)
  dev.off()
  
  ## wilcox test
  test_result <- wilcox.test(plot_data_df$y_plot[plot_data_df$x_plot == "Stage III/IV"], plot_data_df$y_plot[plot_data_df$x_plot == "Stage I/II"])
  wilcox_pvalues_vec <- c(wilcox_pvalues_vec, test_result$p.value)
  ## t test
  test_result <- t.test(plot_data_df$y_plot[plot_data_df$x_plot == "Stage III/IV"], plot_data_df$y_plot[plot_data_df$x_plot == "Stage I/II"])
  ttest_pvalue_vec <- c(ttest_pvalue_vec, test_result$p.value)
}
result_df <- data.frame(name_score = scoregroups_test, 
                        fdr.wilcox = p.adjust(p = wilcox_pvalues_vec, method = "fdr"), pvalue.wilcox = wilcox_pvalues_vec,
                        fdr.ttest = p.adjust(p = ttest_pvalue_vec, method = "fdr"), pvalue.ttest = ttest_pvalue_vec)

# save output -------------------------------------------------------------
file2write <- paste0(dir_out_now, "wilcox_ttest_results.tsv")
write.table(x = result_df, file = file2write, quote = F, sep = "\t", row.names = F)
