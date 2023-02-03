# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  Circos plot showing the genome-wide chromatin accessibility changes associated with BAP1 mutation
# Section:      Results - Chromatin accessibility changes in BAP1 and PBRM1 mutant tumors
#=======================================================================================

# Load libraries  and set working directory  -----------------------------------
## make sure to install old version (0.4.10) of circlize, otherwise the text is weirdly placed
remove.packages("circlize")
install.packages("remotes")
library(remotes)
install_version("circlize", "0.4.10")
## load libraries
packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
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
## input differentially accessible peaks associated with BAP1
daps_df <- fread(data.table = F, input = "../../data/BAP1_vs_NonMutant_DAP2Gene.EnhancerPromoter.20211011.v1.tsv.gz")
## input differentially expressed genes and proteins associated with BAP1
degs_df <- fread(data.table = F, input = "../../data/BAP1_vs_NonMutants.DEGs.Chromosome_Regions_Annotated.20211011.v1.tsv.gz")
## input the DEGs overlapping with DAPs
dap2deg_df <- fread(data.table = F, input = "../../data/BAP1_vs_NonMutant_DAP2DEG.20211011.v1.tsv.gz")

# set parameters ----------------------------------------------------------
cutoff_log2FC_snATAC <- 0.3
cutoff_log2FC_snRNA <- 0.3
cutoff_fdr_bulkRNA <- 0.0001
cutoff_log2FC_bulkRNA <- 1

# prepare plot data for DAPs-------------------------------------------------------
daps_bed1 <- daps_df %>%
  filter(!is.na(avg_log2FC)) %>%
  filter(DAP_direction == "Down") %>%
  mutate(chr = str_split_fixed(string = peak, pattern = "\\-", n = 3)[,1]) %>%
  mutate(start = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,2])) %>%
  mutate(end = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,3])) %>%
  mutate(value = ifelse(avg_log2FC < -1.5, -1.5, avg_log2FC)) %>%
  select(chr, start, end, value)

daps_bed2 <- daps_df %>%
  filter(!is.na(avg_log2FC)) %>%
  filter(DAP_direction == "Up") %>%
  mutate(chr = str_split_fixed(string = peak, pattern = "\\-", n = 3)[,1]) %>%
  mutate(start = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,2])) %>%
  mutate(end = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,3])) %>%
  mutate(value = ifelse(avg_log2FC > 1.5, 1.5, avg_log2FC)) %>%
  select(chr, start, end, value)
daps_bed_list = list(daps_bed1, daps_bed2)
## write plot data
write.table(x = daps_bed1, file = "../../plot_data/F6b.snATAC.Down.SourceData.tsv", quote = F, sep = "\t", row.names = F)
write.table(x = daps_bed2, file = "../../plot_data/F6b.snATAC.Up.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# prepare plot data for snRNA DEGs-------------------------------------------------------
degs_snrna_bed1 <- degs_df %>%
  filter(!is.na(Num_sig_up.snRNA)) %>%
  filter(!is.na(FDR.snRNA.cnvcorrected) & FDR.snRNA.cnvcorrected < 0.05) %>%
  # filter(avg_log2FC.snRNA < -cutoff_log2FC_snRNA) %>%
  filter(Num_sig_down.snRNA >=5) %>%
  # filter(Num_up == 0) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  rename(value = avg_log2FC.snRNA) %>%
  select(chr, start, end, value)
degs_snrna_bed2 <- degs_df %>%
  filter(!is.na(Num_sig_up.snRNA)) %>%
  filter(!is.na(FDR.snRNA.cnvcorrected) & FDR.snRNA.cnvcorrected < 0.05) %>%
  # filter(avg_log2FC.snRNA > cutoff_log2FC_snRNA) %>%
  filter(Num_sig_up.snRNA >=5) %>%
  # filter(Num_up == 0) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  rename(value = avg_log2FC.snRNA) %>%
  select(chr, start, end, value)
degs_snrna_bed_list = list(degs_snrna_bed1, degs_snrna_bed2)
## write plot data
write.table(x = degs_snrna_bed1, file = "../../plot_data/F6b.snRNA.Down.SourceData.tsv", quote = F, sep = "\t", row.names = F)
write.table(x = degs_snrna_bed2, file = "../../plot_data/F6b.snRNA.Up.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# prepare plot data for bulk protein-------------------------------------------------------
deps_bed1 <- degs_df %>%
  filter(!is.na(FDR.bulkpro) & FDR.bulkpro < 0.05) %>%
  filter(meddiff_exp.bulkpro < 0) %>%
  filter(!is.na(start_position)) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  mutate(value = ifelse(meddiff_exp.bulkpro < -0.5, -0.5, meddiff_exp.bulkpro)) %>%
  select(chr, start, end, value, genesymbol_deg)
deps_bed2 <- degs_df %>%
  filter(!is.na(FDR.bulkpro) & FDR.bulkpro < 0.05) %>%
  filter(meddiff_exp.bulkpro > 0) %>%
  filter(!is.na(start_position)) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  mutate(value = ifelse(meddiff_exp.bulkpro > 0.5, 0.5, meddiff_exp.bulkpro)) %>%
  select(chr, start, end, value, genesymbol_deg)
deps_bed_list = list(deps_bed1, deps_bed2)
## write plot data
write.table(x = deps_bed1, file = "../../plot_data/F6b.Protein.Down.SourceData.tsv", quote = F, sep = "\t", row.names = F)
write.table(x = deps_bed2, file = "../../plot_data/F6b.Protein.Up.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# prepare plot data for bulk RNA-------------------------------------------------------
degs_bulkrna_bed1 <- degs_df %>%
  filter(!is.na(FDR.bulkRNA) & FDR.bulkRNA < cutoff_fdr_bulkRNA) %>%
  filter(logFC.bulkRNA < -(cutoff_log2FC_bulkRNA)) %>%
  filter(!is.na(start_position)) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  mutate(value = ifelse(logFC.bulkRNA < -5, -5, logFC.bulkRNA)) %>%
  select(chr, start, end, value, genesymbol_deg)
degs_bulkrna_bed2 <- degs_df %>%
  filter(!is.na(FDR.bulkRNA) & FDR.bulkRNA < cutoff_fdr_bulkRNA) %>%
  filter(logFC.bulkRNA > cutoff_log2FC_bulkRNA) %>%
  filter(!is.na(start_position)) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  mutate(value = ifelse(logFC.bulkRNA > 5, 5, logFC.bulkRNA)) %>%
  select(chr, start, end, value, genesymbol_deg)
degs_bulkrna_bed_list = list(degs_bulkrna_bed1, degs_bulkrna_bed2)
## write plot data
write.table(x = degs_bulkrna_bed1, file = "../../plot_data/F6b.bulkRNA.Down.SourceData.tsv", quote = F, sep = "\t", row.names = F)
write.table(x = degs_bulkrna_bed2, file = "../../plot_data/F6b.bulkRNA.Up.SourceData.tsv", quote = F, sep = "\t", row.names = F)

# prepare plot data for highlighting genes-------------------------------------------------------
genes_highlight <- dap2deg_df %>%
  dplyr::filter(abs(avg_log2FC.snATAC) >= cutoff_log2FC_snATAC & abs(avg_log2FC.snRNA) >= cutoff_log2FC_snRNA) %>%
  dplyr::filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA) & (avg_log2FC.snATAC >= 0 & avg_log2FC.snRNA >= 0 | avg_log2FC.snATAC <= 0 & avg_log2FC.snRNA <= 0))
# genes_highlight <- dap2deg_df %>%
#   filter(abs(avg_log2FC.snATAC) >= 0.5 & abs(avg_log2FC.snRNA) >= 0.5)
genes_validatedbybulkRNA_df <- degs_df %>%
  filter(genesymbol_deg %in% genes_highlight$Gene) %>%
  filter(!is.na(Num_sig_up.snRNA)) %>%
  filter(!is.na(FDR.snRNA.cnvcorrected) & FDR.snRNA.cnvcorrected < 0.05) %>%
  filter(!is.na(FDR.bulkRNA) & FDR.bulkRNA < cutoff_fdr_bulkRNA) %>%
  filter((Num_sig_up.snRNA >=5 & logFC.bulkRNA > 0) | (Num_sig_down.snRNA >=5 & logFC.bulkRNA < 0)) %>%
  filter(!is.na(start_position)) %>%
  mutate(chr = paste0("chr", chromosome_name))
genes_validatedbybulkpro_df <- degs_df %>%
  filter(genesymbol_deg %in% genes_highlight$Gene) %>%
  filter(!is.na(Num_sig_up.snRNA)) %>%
  filter(!is.na(FDR.snRNA.cnvcorrected) & FDR.snRNA.cnvcorrected < 0.05) %>%
  filter(!is.na(FDR.bulkpro) & FDR.bulkpro < 0.05) %>%
  filter((Num_sig_up.snRNA >=5 & meddiff_exp.bulkpro > 0) | (Num_sig_down.snRNA >=5 & meddiff_exp.bulkpro < 0)) %>%
  filter(!is.na(start_position))

highlight_bed <- degs_df %>%
  filter(genesymbol_deg %in% genes_highlight$Gene) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  mutate(value = genesymbol_deg) %>%
  mutate(validated_by_bulkRNA = (genesymbol_deg %in% genes_validatedbybulkRNA_df$genesymbol_deg)) %>%
  mutate(validated_by_bulkprotein = (genesymbol_deg %in% genes_validatedbybulkpro_df$genesymbol_deg)) %>%
  mutate(deg_group = ifelse(validated_by_bulkRNA, 
                            ifelse(validated_by_bulkprotein, "bulk RNA & bulk protein", "bulk RNA only "),
                            ifelse(validated_by_bulkprotein, "bulk protein only", "not validated by bulk"))) %>%
  mutate(deg_color = ifelse(validated_by_bulkRNA | validated_by_bulkprotein, 
                            ifelse(avg_log2FC.snRNA > 0, "#E41A1C" , "#377EB8"),
                            ifelse(avg_log2FC.snRNA > 0, "#FB9A99" , "#A6CEE3"))) %>%
  mutate(deg_font = ifelse(!validated_by_bulkRNA & !validated_by_bulkprotein, 3, 4))  %>%
  select(chr, start, end, value, deg_group, deg_color, deg_font)

# plot all data tracks --------------------------------------------------------------------
colors_datatypes <- c("#99D8C9", "#FFFFB3", "#FDD0A2",  "#DECBE4")
names(colors_datatypes) <- c("snATAC", "snRNA", "bulkRNA", "bulkprotein")
## initialize
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out, "F6b_Circos_BAP1_DAP_DEG_DEP", ".pdf")
pdf(file2write, width = 7, height = 7, useDingbats = F)
circos.par(gap.degree = c(rep(1, 23), 40), start.degree = 90)
circos.initializeWithIdeogram(plotType = NULL)
## plot labels
# circos.genomicLabels(highlight_bed[1:10,], labels.column = 4, side = "outside")
circos.genomicLabels(highlight_bed, labels.column = 4, side = "outside", connection_height = 0.08, labels_height = 0.08, cex = 0.8, col = highlight_bed$deg_color, line_lwd = 0.5, font = highlight_bed$deg_font)
## plot track
circos.genomicTrack(daps_bed_list, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      color_i <- c("#377EB8", "#E41A1C")[i]
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
                    }, track.height = 0.14, bg.col = colors_datatypes["snATAC"], bg.border = NA)
circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8, at = c(-1.5, 0, 1.5))
circos.genomicTrack(degs_snrna_bed_list, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      color_i <- c("#377EB8", "#E41A1C")[i]
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
                    }, track.height = 0.14, bg.col = colors_datatypes["snRNA"], bg.border = NA)
circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8, at = c(-1.5, 0, 1.5))
circos.genomicTrack(degs_bulkrna_bed_list, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      color_i <- c("#377EB8", "#E41A1C")[i]
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
                    }, track.height = 0.14, bg.col = colors_datatypes["bulkRNA"], bg.border = NA)
circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8, at = c(-5, 0, 5))

circos.genomicTrack(deps_bed_list,
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      color_i <- c("#377EB8", "#E41A1C")[i]
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
                    }, track.height = 0.14, bg.col = colors_datatypes["bulkprotein"], bg.border = NA)
circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8, at = c(-0.5, 0, 0.5))
circos.clear()
dev.off()

