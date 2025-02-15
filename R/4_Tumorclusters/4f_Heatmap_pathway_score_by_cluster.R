# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Heatmap showing the gene set scores for 90 tumor subclusters (columns). Tumor subclusters are grouped by patient separated by white lines.
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
  "ComplexHeatmap",
  "circlize",
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
## input meta data
id_metadata_df <- fread(data.table = F, input = "../../data/SampleID.MetaData.tsv")
## input pathway scores
scores_df <- fread(data.table = F, input = "../../data/Pathwayscores_byTumorcluster.MSigDB.Hallmark.tsv.gz")
## input the pathways to plot
genesets_plot_df <- fread(data.table = F, input = "../../data/Genesets_enriched_in_intrapatienttumorclusterDEGs.20220606.v1.tsv.gz")
## input the annotation for the hallmark gene sets
hallmark_anno_df <- fread(data.table = F, input = "../../data/Hallmark_gene_sets_summary.csv")
## input CNV data
cnv_df <- fread(data.table = F, input = "../../data/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20210806.v1.tsv.gz")
## input mutation mapped
driver_mutation_bytumorcluster_df <- fread(data.table = T, input = "../../data/Driver_mutation_mapped_per_intrapatienttumorcluster.20220610.v1.tsv")
## input te bulk mutation result
bulk_mut_bycase_df <- fread(input = "../../data/PBRM1_BAP1_Mutation_Status_By_Case.20210412.v1.tsv", data.table = F)
## input cell fraction by cluster
count_bycluster_df <- fread(data.table = F, input = "../../data/Fraction_cluster_by_case.20220616.v1.tsv")

# preprocess --------------------------------------------------------------
genesets_plot_df <- genesets_plot_df %>%
  mutate(scoregroup_name = paste0(gsub(x = Description, pattern = "HALLMARK_", replacement = ""), "_Score")) %>%
  filter(Description != "HALLMARK_SPERMATOGENESIS")
genesets_plot <- genesets_plot_df$scoregroup_name
hallmark_anno_df$`Hallmark Name`[hallmark_anno_df$`Hallmark Name` == "UV_RESPONSE_DOWN"] <- "UV_RESPONSE_DN"
## make cluster id
clustername_df <- data.frame(cluster_name = scores_df$cluster_name)
clustername_df <- clustername_df %>%
  mutate(sampleid = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
  mutate(sampleid = gsub(x = sampleid, pattern = "\\.", replacement = "-"))
clustername_df$case <- mapvalues(x = clustername_df$sampleid, from = id_metadata_df$Aliquot.snRNA.WU, to = as.vector(id_metadata_df$Case))
clustername_df <- clustername_df %>%
  filter(case != "C3L-00359") %>%
  filter(!(cluster_name %in% c("C3N.00733.T2_C5", "C3L.01313.T1_C7" , "C3L.01287.T1_C2")))
rownames(scores_df) <- scores_df$cluster_name
## prepare CNV data
cnv_df <- cnv_df %>%
  mutate(tumor_subcluster.dataname = gsub(x = tumor_subcluster, pattern = "\\-", replacement = "."))
# format expression data --------------------------------------------------
## get dim names
plot_data_t_mat <- as.matrix(scores_df[,genesets_plot])
plot_data_mat <- t(plot_data_t_mat)
colnames(plot_data_mat) <- scores_df$cluster_name
## make row label
rownames_plot <- rownames(plot_data_mat)
rowlabels_plot <- gsub(x = rownames_plot, pattern = "_Score", replacement = "")

# make column order -------------------------------------------------------
plot_data_mat <- plot_data_mat[, clustername_df$cluster_name]
colnames_plot <- colnames(plot_data_mat)
clusternames_column <- gsub(x = colnames_plot, pattern = "\\.", replacement = "-")
sampleids_column <- str_split_fixed(string = clusternames_column, pattern = "_", n = 2)[,1]
caseids_column <- clustername_df$case
## write source data
write.table(x = plot_data_mat, file = "../../plot_data/F4f.SourceData.tsv", quote = F, sep = "\t", row.names = T)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## colors for gene set scores
summary(as.vector(unlist(plot_data_mat)))
quantile(as.vector(unlist(plot_data_mat)), 0.95)
quantile(as.vector(unlist(plot_data_mat)), 0.05)
colors_heatmapbody = colorRamp2(c(quantile(as.vector(unlist(plot_data_mat)), 0.05), 
                                  0, 
                                  quantile(as.vector(unlist(plot_data_mat)), 0.95)), 
                                c("purple", "black", "yellow"))
colors_truefalse <- c("black", "white")
names(colors_truefalse) <- c("TRUE", "FALSE")
color_gridline = "grey90"
colors_topbottom <- c("red", "blue", "grey90")
names(colors_topbottom) <- c("top", "bottom", "middle")
colors_loss_frac <- colorRamp2(seq(0, 1, 0.2), c("white", brewer.pal(n = 6, name = "Blues")[-1]))
colors_gain_frac <- colorRamp2(seq(0, 1, 0.2), c("white", brewer.pal(n = 6, name = "Reds")[-1]))

# make column annotation --------------------------------------------------
## prepare data
inflam_score_vec <- scores_df[colnames_plot, "INFLAMMATORY_RESPONSE_Score"]
inflam_assign_vec <- ifelse(inflam_score_vec >= quantile(inflam_score_vec, probs = 0.9), "top",
                            ifelse(inflam_score_vec <= quantile(inflam_score_vec, probs = 0.1), "bottom", "middle"))
VHL_Loss_frac_vec <- mapvalues(x = colnames_plot, from = cnv_df$tumor_subcluster.dataname[cnv_df$gene_symbol == "VHL" & cnv_df$cna_3state == "Loss"], to = as.vector(cnv_df$Fraction[cnv_df$gene_symbol == "VHL" & cnv_df$cna_3state == "Loss"]))
VHL_Loss_frac_vec[VHL_Loss_frac_vec == colnames_plot] <- "0"; VHL_Loss_frac_vec <- as.numeric(VHL_Loss_frac_vec)
setd2_loss_frac_vec <- mapvalues(x = colnames_plot, from = cnv_df$tumor_subcluster.dataname[cnv_df$gene_symbol == "SETD2" & cnv_df$cna_3state == "Loss"], to = as.vector(cnv_df$Fraction[cnv_df$gene_symbol == "SETD2" & cnv_df$cna_3state == "Loss"]))
setd2_loss_frac_vec[setd2_loss_frac_vec == colnames_plot] <- "0"; setd2_loss_frac_vec <- as.numeric(setd2_loss_frac_vec)
SQSTM1_Gain_frac_vec <- mapvalues(x = colnames_plot, from = cnv_df$tumor_subcluster.dataname[cnv_df$gene_symbol == "SQSTM1" & cnv_df$cna_3state == "Gain"], to = as.vector(cnv_df$Fraction[cnv_df$gene_symbol == "SQSTM1" & cnv_df$cna_3state == "Gain"]))
SQSTM1_Gain_frac_vec[SQSTM1_Gain_frac_vec == colnames_plot] <- "0"; SQSTM1_Gain_frac_vec <- as.numeric(SQSTM1_Gain_frac_vec)
GOLPH3_Gain_frac_vec <- mapvalues(x = colnames_plot, from = cnv_df$tumor_subcluster.dataname[cnv_df$gene_symbol == "GOLPH3" & cnv_df$cna_3state == "Gain"], to = as.vector(cnv_df$Fraction[cnv_df$gene_symbol == "GOLPH3" & cnv_df$cna_3state == "Gain"]))
GOLPH3_Gain_frac_vec[GOLPH3_Gain_frac_vec == colnames_plot] <- "0"; GOLPH3_Gain_frac_vec <- as.numeric(GOLPH3_Gain_frac_vec)
## prepare cluster-level mutatoin mappping result
mut_map_vec <- mapvalues(x = clusternames_column, from = driver_mutation_bytumorcluster_df$Cluster_Name, to = as.vector(driver_mutation_bytumorcluster_df$number_cells_w_driver_mutation))
mut_map_vec <- as.character(mut_map_vec != "0")
## prepare bulk mutation status
VHL_bycase_vec <- mapvalues(x = caseids_column, from = bulk_mut_bycase_df$Case, to = as.vector(bulk_mut_bycase_df$VHL)); VHL_bycase_vec <- as.character(VHL_bycase_vec != "")
PBRM1_bycase_vec <- mapvalues(x = caseids_column, from = bulk_mut_bycase_df$Case, to = as.vector(bulk_mut_bycase_df$PBRM1)); PBRM1_bycase_vec <- as.character(PBRM1_bycase_vec != "")
BAP1_bycase_vec <- mapvalues(x = caseids_column, from = bulk_mut_bycase_df$Case, to = as.vector(bulk_mut_bycase_df$BAP1)); BAP1_bycase_vec <- as.character(BAP1_bycase_vec != "")
SETD2_bycase_vec <- mapvalues(x = caseids_column, from = bulk_mut_bycase_df$Case, to = as.vector(bulk_mut_bycase_df$SETD2)); SETD2_bycase_vec <- as.character(SETD2_bycase_vec != "")
## prepare data for cell fraction
fraction_bycase_vec <- mapvalues(x = clusternames_column, from = count_bycluster_df$Cluster_Name, to = as.vector(count_bycluster_df$frac_cluster_count_bycase)); fraction_bycase_vec <- as.numeric(fraction_bycase_vec)
colors_cellfrac = colorRamp2(c(quantile(fraction_bycase_vec, 0.1), 
                               quantile(fraction_bycase_vec, 0.3), 
                               quantile(fraction_bycase_vec, 0.5),
                               quantile(fraction_bycase_vec, 0.7),
                               quantile(fraction_bycase_vec, 0.9)), 
                             brewer.pal(n = 5, "YlGn"))
## make column annotation object
fontsize_plot <- 21
colanno_obj <- HeatmapAnnotation(fraction_in_patient = anno_simple(x = fraction_bycase_vec, col = colors_cellfrac, height = unit(0.7, "cm")),
                                 chr3p_VHL_loss_bycluster = anno_simple(x = VHL_Loss_frac_vec, col = colors_loss_frac, height = unit(0.7, "cm")),
                                 chr3p_SETD2_loss_bycluster = anno_simple(x = setd2_loss_frac_vec, col = colors_loss_frac, height = unit(0.7, "cm")),
                                 chr5q_SQSTM1_gain_bycluster = anno_simple(x = SQSTM1_Gain_frac_vec, col = colors_gain_frac, height = unit(0.7, "cm")),
                                 chr5q_GOLPH3_gain_bycluster = anno_simple(x = GOLPH3_Gain_frac_vec, col = colors_gain_frac, height = unit(0.7, "cm")),
                                 #chr3q_MECOM_gain_bycluster = anno_simple(x = MECOM_Gain_frac_vec, col = colors_gain_frac),
                                 driver_mutation_bycluster = anno_simple(x = mut_map_vec, col = colors_truefalse[mut_map_vec], height = unit(0.7, "cm")),
                                 VHL_mutated_bycase = anno_simple(x = VHL_bycase_vec, col = colors_truefalse[VHL_bycase_vec], height = unit(0.7, "cm")),
                                 PBRM1_mutated_bycase = anno_simple(x = PBRM1_bycase_vec, col = colors_truefalse[PBRM1_bycase_vec], height = unit(0.7, "cm")),
                                 BAP1_mutated_bycase = anno_simple(x = BAP1_bycase_vec, col = colors_truefalse[BAP1_bycase_vec], height = unit(0.7, "cm")),
                                 SETD2_mutated_bycase = anno_simple(x = SETD2_bycase_vec, col = colors_truefalse[SETD2_bycase_vec], height = unit(0.7, "cm")),
                                 Inflammatory_score_group = anno_simple(x = inflam_assign_vec, col = colors_topbottom[inflam_assign_vec], height = unit(0.7, "cm")),
                                 annotation_name_side = "left", annotation_name_gp = gpar(fontsize = fontsize_plot))
# ## merge data
# colanno_df <- enrich_plot_df %>%
#   mutate(EMT_Module_Enriched = (cluster_name %in% emt_group_df$cluster_name[emt_group_df$epithelial_group == "EMT"]))
# rownames(colanno_df) <- colanno_df$cluster_name
# colanno_df <- colanno_df[colnames_plot,-1]
# ## make colors
# colors_scores_list <- list()
# for (colname_tmp in colnames(enrich_plot_df)[-1]) {
#   colanno_df[, colname_tmp] <- as.character(colanno_df[, colname_tmp])
#   colors_tmp <- c(colors_isenriched["FALSE"], colors_enrich_type[colname_tmp])
#   names(colors_tmp) <- c("FALSE", "TRUE")
#   colors_scores_list[[colname_tmp]] <- colors_tmp
# }
# colanno_obj = HeatmapAnnotation(df = colanno_df, col = colors_scores_list,
#                                 annotation_name_gp = gpar(fontsize = 14), annotation_name_side = "left", show_legend = F)


# make row annotation -----------------------------------------------------
freq_de_vec <- mapvalues(x = rownames_plot, from = genesets_plot_df$scoregroup_name, to = as.vector(genesets_plot_df$Freq)); freq_de_vec <- as.numeric(freq_de_vec)
rowanno_obj <- rowAnnotation(Freq_of_DE = anno_barplot(freq_de_vec, width = unit(1.75, "cm"), axis_param = list(gp = gpar(fontsize = fontsize_plot))), annotation_name_side = "top", annotation_name_gp = gpar(fontsize = fontsize_plot))


# make gene set split -----------------------------------------------------
row_split_vec <- mapvalues(x = rowlabels_plot, from = hallmark_anno_df$`Hallmark Name`, to = as.vector(hallmark_anno_df$`Process Category`))


# make column split -------------------------------------------------------
column_split_vec <- mapvalues(x = colnames_plot, from = clustername_df$cluster_name, to = as.vector(clustername_df$case))
clustercount_df <- clustername_df %>%
  group_by(case) %>%
  summarise(number_clusters = n()) %>%
  arrange(desc(number_clusters))
column_split_factor <- factor(x = column_split_vec, levels = clustercount_df$case)

# plot  ------------------------------------------------------------
p <- Heatmap(matrix = plot_data_mat, 
             col = colors_heatmapbody,
             na_col = color_na, #border = "black",
             cell_fun = function(j, i, x, y, w, h, fill) {
               if (plot_data_mat[i,j] >= (quantile(plot_data_mat[i,], 0.75)+1.5*IQR(plot_data_mat[i,]))) {
                 grid.text("*", x, y)
               }
               # if (plot_data_mat[i,j] >= max(plot_data_mat[,j])) {
               #   # if (plot_data_mat[i,j] >= min(tail(sort(plot_data_mat[i,]), 2))) {
               #   grid.rect(x = x, y = y, width = w, height = h,
               #             gp = gpar(col = "red", fill = NA))
               # }
             },
             width = ncol(plot_data_mat)*unit(4, "mm"), 
             # height = nrow(plot_data_mat)*unit(5, "mm"),
             ## row
             show_row_names = T, row_names_gp = gpar(fontsize = fontsize_plot), row_names_side = "right",
             show_row_dend = T, row_dend_side = "left", cluster_row_slices = T, 
             row_split = row_split_vec, row_title_side = "left", row_title_rot = 0, row_title_gp = gpar(fontsize = 25),
             row_labels = tolower(rowlabels_plot), 
             right_annotation = rowanno_obj,
             ## column
             show_column_dend = F, cluster_columns = T, 
             column_split = column_split_factor, cluster_column_slices = F, column_title = NULL,
             #column_title_rot = 90,
             top_annotation = colanno_obj, 
             show_column_names = F, #column_names_side = "top", column_names_gp = gpar(fontsize = 5),
             show_heatmap_legend = F)
list_lgd = list(
  Legend(col_fun = colors_heatmapbody, 
         title = "Gene set score", 
         title_gp = gpar(fontsize = 14),
         labels_gp = gpar(fontsize = 14),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_cellfrac, 
         title = "% cluster cells in\npatient tumor cells", 
         title_gp = gpar(fontsize = 14),
         labels_gp = gpar(fontsize = 14),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_loss_frac, 
         title = "% cells in cluster\nwith CN loss", 
         title_gp = gpar(fontsize = 14),
         labels_gp = gpar(fontsize = 14),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_gain_frac, 
         title = "% cells in cluster\nwith CN gain", 
         title_gp = gpar(fontsize = 14),
         labels_gp = gpar(fontsize = 14),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"),
  Legend(labels = names(colors_truefalse), labels_gp = gpar(fontsize = 14),
         title = "Mutation status", title_gp = gpar(fontsize = 14),
         legend_gp = gpar(fill = colors_truefalse)),
  Legend(labels = names(colors_topbottom), labels_gp = gpar(fontsize = 14),
         title = "Inflammatory  score group", title_gp = gpar(fontsize = 14),
         legend_gp = gpar(fill = colors_topbottom)))

# write output ------------------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out, "F4f_Heatmap_pathway_score_by_cluster", ".pdf")
pdf(file2write, width = 25, height = 14, useDingbats = F)
draw(object = p,
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()
