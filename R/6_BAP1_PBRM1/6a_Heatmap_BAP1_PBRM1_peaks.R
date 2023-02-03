# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Heatmap showing the relative changes in ATAC-peak accessibility for peaks differentially accessible between the tumor cells of BAP1-mutated tumors (6 tumors, including 2 BAP1- and PBRM1-mutated tumors, 29,366 cells) vs. non-BAP1/PBRM1-mutated tumors (8 tumors; non-mutants) and peaks differentially accessible between tumor cells of PBRM1-mutated tumors (9 tumors, 32,255 cells) vs. non-BAP1/PBRM1-mutated tumors (non-mutants)
# Section:      Results - Chromatin accessibility changes in BAP1 and PBRM1 mutant tumors
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
## input accessibility data
acc_pbrm1_down_df=fread(data.table = F, input = '../../data/DOWN_PBRM1mutants_vs_nonMutans.Accessibility.20211011.tsv.gz')
acc_pbrm1_up_df=fread(data.table = F, input ='../../data/UP_PBRM1mutants_vs_nonMutans.Accessibility.20211011.tsv.gz')
acc_bap1_down_df=fread(data.table = F, input = '../../data/DOWN_BAP1mutants_vs_nonMutans.Accessibility.20211011.tsv.gz')
acc_bap1_up_df=fread(data.table = F, input ='../../data/UP_BAP1mutants_vs_nonMutans.Accessibility.20211011.tsv.gz')
## input sample category
mut_df <- fread(data.table = F, input = "../../data/PBRM1_BAP1_Mutation_Status_By_Case.20210412.v1.tsv")

# make data matrix --------------------------------------------------------
## add row names
row.names(acc_pbrm1_down_df)=paste(acc_pbrm1_down_df$Group.1,acc_pbrm1_down_df$Group.2,sep='_')
row.names(acc_pbrm1_up_df)=paste(acc_pbrm1_up_df$Group.1,acc_pbrm1_up_df$Group.2,sep='_')
row.names(acc_bap1_down_df)=paste(acc_bap1_down_df$Group.1,acc_bap1_down_df$Group.2,sep='_')
row.names(acc_bap1_up_df)=paste(acc_bap1_up_df$Group.1,acc_bap1_up_df$Group.2,sep='_')
## make row names consistent
rownames_plot <- rownames(acc_pbrm1_down_df)
## cbind
plotdata_df <- cbind(acc_pbrm1_down_df, 
                     acc_pbrm1_up_df[rownames_plot, grepl(pattern = "chr", x = colnames(acc_pbrm1_up_df))],
                     acc_bap1_down_df[rownames_plot, grepl(pattern = "chr", x = colnames(acc_bap1_down_df))],
                     acc_bap1_up_df[rownames_plot, grepl(pattern = "chr", x = colnames(acc_bap1_up_df))]
)
plotdata_df <- plotdata_df[,!duplicated(colnames(plotdata_df))]
## change group names
plotdata_df$Group.1=gsub('_','-',plotdata_df$Group.1)
plotdata_df$Group.2=gsub('_','.',plotdata_df$Group.2)
## make ID column
plotdata_df$ID=paste(plotdata_df$Group.1,plotdata_df$Group.2,sep='_')
## change row names
rownames(plotdata_df)=plotdata_df$ID
## filter by row names
plotdata_df <- plotdata_df[!(grepl(x = rownames(plotdata_df), pattern = "\\-N") & grepl(x = rownames(plotdata_df), pattern = "Tumor")),]
plotdata_df <- plotdata_df[!(grepl(x = rownames(plotdata_df), pattern = "\\-T") & grepl(x = rownames(plotdata_df), pattern = "PT")),]
colnames(plotdata_df)[1:2]=c('Sample','Cell_type')
plotdata_df <- plotdata_df[plotdata_df$Cell_type %in% c('Tumor','PT') & !(plotdata_df$Sample %in% c("C3L-01287-T1")),]
plotdata_df2 <- plotdata_df[, colnames(plotdata_df)[grepl(pattern = "chr", x = colnames(plotdata_df))]]
plotdata_mat <- scale(x = plotdata_df2)
rownames(plotdata_mat) <- plotdata_df$Sample
## re-order rows
rownames_ordered <- c("C3L-00416-T2","C3L-00908-T1", 
                      "C3N-01200-T1", "C3L-01313-T1", "C3N-00317-T1","C3N-00437-T1", #"C3L-01287-T1",
                      "C3N-01213-T1","C3L-00079-T1", "C3L-00610-T1", "C3N-00242-T1","C3L-00790-T1","C3L-00583-T1", "C3L-00004-T1", "C3N-00733-T1","C3L-01302-T1", 
                      "C3L-00448-T1",  "C3L-00917-T1", "C3L-00088-T1", "C3L-00088-T2","C3L-00010-T1", "C3L-00026-T1", "C3N-00495-T1","C3L-00096-T1", 
                      "C3N-01200-N", "C3L-00088-N", "C3L-00079-N", "C3N-00242-N")
plotdata_mat <- plotdata_mat[rownames_ordered,]
## write plot data
write.table(x = plotdata_mat, file = "../../plot_data/F6a.SourceData.tsv", quote = F, sep = "\t", row.names = T)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## colors for gene set scores
summary(as.numeric(unlist(plotdata_mat)))
color_red <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[2]
colors_heatmapbody = colorRamp2(seq(-2, 2, 0.4), 
                                rev(brewer.pal(n = 11, name = "BrBG")))
colors_bap1_vaf <- colorRamp2(c(0, 0.5), c("white smoke", '#984EA3'))
colors_pbrm1_vaf <- colorRamp2(c(0, 0.5), c("white smoke", '#FF7F00'))
colors_peaktype <- c("Up-regulated" = "#E41A1C", "Down-regulated" = "#377EB8", "other" = "white smoke")

# make row annotation -----------------------------------------------------
## get BAP1 mutated cases
row_anno_df=plotdata_df %>%
  select(Sample, Cell_type) %>%
  mutate(Case = gsub(x = Sample, pattern = "\\-[A-Z][0-9]", replacement = ""))
row_anno_df$BAP1_mutation <- mapvalues(x = row_anno_df$Case, from = mut_df$Case, to = as.vector(mut_df$BAP1))
row_anno_df$BAP1_mutation[row_anno_df$BAP1_mutation == row_anno_df$Case] <- ""
row_anno_df$PBRM1_mutation <- mapvalues(x = row_anno_df$Case, from = mut_df$Case, to = as.vector(mut_df$PBRM1))
row_anno_df$PBRM1_mutation[row_anno_df$PBRM1_mutation == row_anno_df$Case] <- ""
row_anno_df$mutation_type <- mapvalues(x = row_anno_df$Case, from = mut_df$Case, to = as.vector(mut_df$mutation_category_sim))
row_anno_df$mutation_type[row_anno_df$mutation_type == row_anno_df$Case] <- "PT"
row_anno_df <- row_anno_df %>%
  mutate(BAP1_Mut_VAF = as.numeric(str_split_fixed(string = BAP1_mutation, pattern = "\\(|\\)", n = 3)[,2])) %>%
  mutate(PBRM1_Mut_VAF = as.numeric(str_split_fixed(string = PBRM1_mutation, pattern = "\\(|\\)", n = 3)[,2])) %>%
  mutate(BAP1_Mut_VAF = ifelse(is.na(BAP1_Mut_VAF), 0, BAP1_Mut_VAF)) %>%
  mutate(PBRM1_Mut_VAF = ifelse(is.na(PBRM1_Mut_VAF), 0, PBRM1_Mut_VAF))
rownames(row_anno_df) <- row_anno_df$Sample
row_anno_df <- row_anno_df[rownames_ordered,]
row_ha= rowAnnotation(#Cell_type=row_anno_df$Cell_type, 
  BAP1_Mut_VAF=row_anno_df$BAP1_Mut_VAF,
  PBRM1_Mut_VAF=row_anno_df$PBRM1_Mut_VAF,
  col=list(BAP1_Mut_VAF=colors_bap1_vaf,
           PBRM1_Mut_VAF=colors_pbrm1_vaf,
           Cell_type=c('PT'='#1B9E77','Tumor'='#E7298A')), 
  annotation_width = unit(2, "mm"), show_legend = F)

# make column annotation --------------------------------------------------
is_bap1_peaks_vec <- ifelse(colnames(plotdata_mat) %in% colnames(acc_bap1_down_df), "Down-regulated",
                            ifelse(colnames(plotdata_mat) %in% colnames(acc_bap1_up_df), "Up-regulated", "other"))
is_pbrm1_peaks_vec <- ifelse(colnames(plotdata_mat) %in% colnames(acc_pbrm1_down_df), "Down-regulated",
                             ifelse(colnames(plotdata_mat) %in% colnames(acc_pbrm1_up_df), "Up-regulated", "other"))
column_ha <- HeatmapAnnotation(BAP1_associated_peak = anno_simple(x = is_bap1_peaks_vec, col = colors_peaktype[is_bap1_peaks_vec], height = unit(0.6, "cm")),
                               PBRM1_associated_peak = anno_simple(x = is_pbrm1_peaks_vec, col = colors_peaktype[is_pbrm1_peaks_vec],  height = unit(0.6, "cm")),
                               annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 17))

# make column split -------------------------------------------------------
column_split_vec <- ifelse(colnames(plotdata_mat) %in% colnames(acc_pbrm1_down_df), "PBRM1 down peak",
                           ifelse(colnames(plotdata_mat) %in% colnames(acc_pbrm1_up_df), "PBRM1 up peak", 
                                  ifelse(colnames(plotdata_mat) %in% colnames(acc_bap1_down_df), "BAP1 down peak", "BAP1 up peak")))

# make row split -------------------------------------------------------
row_split_vec <- row_anno_df$mutation_type
row_split_vec[row_split_vec == "Both mutated"] <- "BAP1&PBRM1\nmutated"
row_split_factor <- factor(row_split_vec, levels = c("BAP1&PBRM1\nmutated", "BAP1 mutated", 'PBRM1 mutated', "Non-mutants", 'PT'))

# plot --------------------------------------------------------------------
list_lgd = list(
  Legend(title = "Peak accessibility direction vs. non-mutants", 
         title_gp = gpar(fontsize = 12),
         labels = names(colors_peaktype), legend_gp = gpar(fill = colors_peaktype), grid_width = unit(0.5, "cm"), grid_height = unit(0.75, "cm"),
         labels_gp =  gpar(fontsize = 12),direction = "horizontal", nrow = 1),
  Legend(col_fun = colors_heatmapbody, 
         title = "Relative peak\naccessibility",
         title_gp = gpar(fontsize = 12),
         labels_gp = gpar(fontsize = 12),
         legend_width = unit(3, "cm"), 
         direction = "horizontal"),
  Legend(col_fun = colors_bap1_vaf, 
         title = "BAP1 mutation VAF", 
         title_gp = gpar(fontsize = 12),
         labels_gp = gpar(fontsize = 12),
         legend_width = unit(3, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_pbrm1_vaf, 
         title = "PBRM1 mutation VAF", 
         title_gp = gpar(fontsize = 12),
         labels_gp = gpar(fontsize = 12),
         legend_width = unit(3, "cm"),
         direction = "horizontal"))

p=ComplexHeatmap::Heatmap(matrix = plotdata_mat, col = colors_heatmapbody, name = "Peak\naccessibility", 
                          ## rows
                          show_row_names = F,  row_names_gp = gpar(fontsize = 13),
                          show_row_dend=FALSE,  left_annotation=row_ha, row_split = row_split_factor, row_title_rot = 0, row_title_gp = gpar(fontsize = 17),
                          cluster_row_slices=F, cluster_rows = F,
                          ## columns
                          show_column_names = FALSE, show_column_dend=FALSE, column_title = NULL, top_annotation = column_ha,
                          column_split = column_split_vec, cluster_column_slices=F,
                          ## other
                          show_heatmap_legend = F, use_raster = T)

# write output ------------------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out, "F6a.Heatmap.BAP1_PBRM1_peaks", ".pdf")
pdf(file2write, width = 10, height=6.5, useDingbats = F)
draw(object = p,
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()
