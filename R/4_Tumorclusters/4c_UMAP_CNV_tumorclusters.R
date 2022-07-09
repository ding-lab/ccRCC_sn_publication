# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  UMAP showing tumor clusters of three tumor samples from the same patient, colored by the copy number status of VHL and SQSTM1.
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

# input -------------------------------------------------------------------
## barcode-UMAP info for tumor clusters
barcode2umap_df <- fread(data.table = F, input = "../../data/MetaData_TumorCellOnlyReclustered.20210805.v1.tsv.gz")
## barcode to tumor subcluster assignment
barcode2tumorsubcluster_df <- fread(input = "../../data/Barcode2TumorSubclusterId.20210805.v1.tsv.gz", data.table = F)
## scrublet information
barcode2scrublet_df <- fread(input = "../../data/scrublet.united_outputs.20210729.v1.tsv.gz", data.table = F)
## input CNV genes
knowncnvgenes_df <- fread(data.table = F, input = "../../data/Known_CNV_genes.20200528.v1.csv")

# specify samples to process ----------------------------------------------
aliquots2process <- unique(barcode2tumorsubcluster_df$orig.ident[barcode2tumorsubcluster_df$easy_id %in% c("C3N-01200-T1", "C3N-01200-T2","C3N-01200-T3")])
aliquots2process

# make output directory ---------------------------------------------------
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
dir_out_now <- paste0(dir_out, "F4c_UMAP_CNV_tumorclusters", "/")
dir.create(dir_out_now)

# pre-process --------------------------------------------------------------
## set genes to plot
genes2plot <- knowncnvgenes_df$Gene_Symbol
knowncnvgenes_df <- knowncnvgenes_df %>%
  mutate(chr_arm = paste0(Chromosome, Arm))
## set infercnv result direcotory
dir_infercnv_output <- "../../data/InferCNV/outputs/"
## set colors
PuBu_colors <- RColorBrewer::brewer.pal(n = 9, name = "PuBu")
PuRd_colors <- RColorBrewer::brewer.pal(n = 9, name = "PuRd")
copy_number_colors <- c("loss" = PuBu_colors[5],
                        "gain" = PuRd_colors[5],
                        "neutral" = "grey50")
barcode2umap_df <- merge(x = barcode2umap_df, y = barcode2tumorsubcluster_df, 
                         by.x = c("orig.ident", "barcode_tumorcellreclustered", "easy_id"), by.y = c("orig.ident", "barcode", "easy_id"), all.x = T)

# Loop: for each aliquot, input seurat object and infercnv output, plot important genes on UMAP ---------
for (snRNA_aliquot_id_tmp in aliquots2process) {
  ## get the case id for this aliquot to show in the title
  easy_id_tmp <- unique(barcode2umap_df$easy_id[barcode2umap_df$orig.ident == snRNA_aliquot_id_tmp])
  
  ## create output directory by aliquot
  dir_out1 <- paste0(dir_out_now, easy_id_tmp, "/")
  dir.create(dir_out1)
  
  scrublets_df <- barcode2scrublet_df %>%
    filter(Aliquot_WU == easy_id_tmp) %>%
    filter(predicted_doublet)
  barcodes_doublet <- scrublets_df$Barcode; length(barcodes_doublet)
  
  ## get umap coordates
  umap_df <- barcode2umap_df %>%
    filter(orig.ident == snRNA_aliquot_id_tmp) %>%
    filter(!(barcode_tumorcellreclustered %in% barcodes_doublet)) %>%
    mutate(barcode = barcode_tumorcellreclustered) %>%
    mutate(Text_TumorCluster = paste0("C", id_manual_cluster_w0+1))
  
  ## input infercnv CNV state results
  obs_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, snRNA_aliquot_id_tmp, ".infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt.gz"), data.table = F)
  ref_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, snRNA_aliquot_id_tmp, ".infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt.gz"), data.table = F)
  dim(obs_cnv_state_mat)
  dim(ref_cnv_state_mat)
  
  ## transform infercnv result wide data frame to long data frame
  cnv_state_df <- rbind(melt(obs_cnv_state_mat, id.vars = c("V1")), melt(ref_cnv_state_mat, id.vars = c("V1")))
  rm(obs_cnv_state_mat)
  rm(ref_cnv_state_mat)
  
  # for (gene_tmp in c("VHL", "SQSTM1")) {
  for (gene_tmp in genes2plot) {
    chr_arm_tmp <- knowncnvgenes_df$chr_arm[knowncnvgenes_df$Gene_Symbol == gene_tmp]
    ## create output directory by chromosome region
    dir_out2 <- paste0(dir_out1, chr_arm_tmp, "/")
    dir.create(dir_out2, showWarnings = F)
    
    # file2write <- paste(dir_out2, easy_id_tmp, ".", gene_tmp, ".png", sep="")
    file2write <- paste(dir_out2, easy_id_tmp, ".", gene_tmp, ".pdf", sep="")
    if ((gene_tmp %in% cnv_state_df$V1) & !file.exists(file2write)) {
      ## extract current gene related results from the infercnv result data frame
      infercnv_observe_gene_tab <- cnv_state_df %>%
        rename(gene_symbol = V1) %>%
        filter(gene_symbol == gene_tmp) %>%
        mutate(barcode = str_split_fixed(string = variable, pattern = "_", n = 2)[,1]) %>%
        rename(copy_state = value) %>%
        select(gene_symbol, barcode, copy_state)
      
      ## add CNV state to the barcode - UMAP coordidate data frame
      point_data_df <- merge(umap_df, infercnv_observe_gene_tab, by = c("barcode"), all.x = T)
      
      ## map CNV state value to text
      point_data_df$cnv_cat <- map_infercnv_state2category(copy_state = point_data_df$copy_state)
      point_data_df$cnv_cat %>% table()
      
      ## make cells with CNV appear on top
      point_data_df <- point_data_df %>%
        mutate(cnv_cat_simple = ifelse(cnv_cat %in% c("0 Copy", "1 Copy"), "loss",
                                       ifelse(cnv_cat %in% c("3 Copies", "4 Copies", ">4 Copies"), "gain", "neutral"))) %>%
        arrange(factor(cnv_cat, levels = c("2 Copies", ">4 Copies", "0 Copy", "1 Copy", "3 Copies")))
      # arrange(desc(cnv_cat))
      
      ## make text data
      cellnumber_percluster_df <- point_data_df %>%
        select(Text_TumorCluster) %>%
        table() %>%
        as.data.frame() %>%
        rename(Text_TumorCluster = ".")
      text_data_df <- point_data_df %>%
        filter(Text_TumorCluster %in% cellnumber_percluster_df$Text_TumorCluster[cellnumber_percluster_df$Freq >= 50]) %>%
        group_by(Text_TumorCluster) %>%
        summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      
      p <- ggplot() +
        geom_point_rast(data = point_data_df, mapping = aes(UMAP_1, UMAP_2, color=cnv_cat), alpha = 0.8, size = 0.2) +
        scale_color_manual(values = copy_number_colors)
      p <- p + geom_text_repel(data = text_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, label = Text_TumorCluster))
      p <- p + theme_bw()
      p <- p + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
      p <- p + theme(legend.position = "none")
      p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                              axis.text.y=element_blank(),axis.ticks=element_blank(),
                              axis.title.x=element_blank(),
                              axis.title.y=element_blank())
      pdf(file2write, width = 2, height = 2, useDingbats = F)
      print(p)
      dev.off()
    }
  }
}



