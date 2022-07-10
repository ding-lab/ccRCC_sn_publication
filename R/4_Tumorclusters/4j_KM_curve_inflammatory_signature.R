# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Kaplan-Meier survival analysis showing overall survival after initial pathological diagnosis
# Section:      Results - Transcriptome-based tumor-cell subclusters may represent genomically distinct subclones
#=======================================================================================

# Load libraries  and set working directory  -----------------------------------
## load libraries
packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "survival",
  "survminer"
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
## input scores
exp_df <- fread(data.table = F, input = "../../data/tumorcell_intrinsic_inflammation_signature_scores.20220615.v1.tsv.gz")
## input survival ddata
survival_df <- fread(data.table = F, input = "../../data/CPTAC_Discovery_ccRCC_Survival_Time20220317.v1.tsv.gz")

# set paramters -----------------------------------------------------------
fontsize_plot = 14

# preprocess ------------------------------------------------------
testdata_df <- merge(x = exp_df, y = survival_df, by.x = "case", by.y = c("CASE_ID"), all.x = T)
testdata_df <- testdata_df %>%
  mutate(Expression = tumorcellintrinsic_inflamm_score)
cutoff_exp_low <- quantile(x = testdata_df$Expression, probs = 0.25, na.rm = T); cutoff_exp_low
cutoff_exp_high <- quantile(x = testdata_df$Expression, probs = 0.75, na.rm = T); cutoff_exp_high
testdata_df <- testdata_df %>%
  mutate(Expression_group = ifelse(Expression <= cutoff_exp_low, "Low", ifelse(Expression >= cutoff_exp_high, "High", "Medium")))
table(testdata_df$Expression_group)

# test overall survival ---------------------------------------------------
## EFS_censor == 0 with event; == 1 without event
## test
testdata_comp_df <- testdata_df %>%
  mutate(surv_time = (OS_time + 9)/365)  %>%
  mutate(surv_status = ifelse(OS_status == "censored", 1, 2)) %>%
  filter(!is.na(surv_status) & !is.na(surv_time) & !is.na(Expression_group)) %>%
  filter(Expression_group != "Medium")
fit_efs <- survfit(Surv(surv_time, surv_status) ~ Expression_group, data = testdata_comp_df)

res <- ggsurvplot(fit_efs,
                  data = testdata_comp_df,
                  conf.int = TRUE,
                  surv.median.line = "hv", pval = TRUE,
                  legend.title = paste0("tumor-cell-intrinsic\ninflammation score\n(mRNA)"),
                  legend.labs = c("High", "Low"), 
                  legend = "top",
                  xlab = "Time (year)",
                  ylab = "Overall Survival",
                  palette = c("#800026", "#FEB24C"),
                  ggtheme = theme_survminer(base_size = fontsize_plot,
                                            base_family = "",
                                            font.main = c(fontsize_plot, "plain", "black"),
                                            font.submain = c(fontsize_plot, "plain", "black"),
                                            font.x = c(fontsize_plot, "plain", "black"),
                                            font.y = c(fontsize_plot, "plain", "black"),
                                            font.caption = c(fontsize_plot, "plain", "black"),
                                            font.tickslab = c(fontsize_plot, "plain", "black"),
                                            legend = c("top", "bottom", "left", "right", "none"),
                                            font.legend = c(fontsize_plot, "plain", "black")),
                  conf.int.alpha = 0.1, tables.height = 0.3,
                  risk.table = TRUE, # Add risk table
                  risk.table.col = "strata", # Change risk table color by groups
                  linetype = "strata") # Change line type by groups
res$table <- res$table + theme(axis.line = element_blank())

# write output ------------------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out, "F4j_KM_curve_inflammatory_signature", ".pdf")
pdf(file2write, width = 3.25, height = 5, useDingbats = F)
print(res)
dev.off()
