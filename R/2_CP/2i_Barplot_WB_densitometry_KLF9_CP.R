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
  "ggpubr",
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

# set plot parameters -----------------------------------------------------
lines_plot <- c("RCC4_scrambled", "RCC4_KLF9_C2")
test_plot <- "t.test"
dir_out <- paste0("../../outputs/"); dir.create(dir_out)
## make colors
colors_byline <- RColorBrewer::brewer.pal(n = 6, name = "Set2")[c(1:2)]
names(colors_byline) <- c("RCC4_scrambled", "RCC4_KLF9_C2")


# plot -----------------------------------------------------------------
for (gene_plot in c("CP", "KLF9")) {
  # ## input data
  # densitometry_df <- fread(data.table = F, input = "../../data/wb_densitometry_alltargets_normalized.10112022.v2.csv")
  # plotdata_df <- densitometry_df %>%
  #   mutate(Date = as.character(Date)) %>%
  #   mutate(Line = paste0(Parent_Line, "_", Construct_transfected)) %>%
  #   filter(Gene_symbol %in% gene_plot) %>%
  #   filter(Line %in% lines_plot) %>%
  #   filter(!is.na(Value_bytub_byscrambled))
  # countsamples_bydate <- table(plotdata_df %>%
  #                                select(Date, Line)) %>% rowSums()
  # dates_filtered <- names(countsamples_bydate)[countsamples_bydate == length(lines_plot)]
  # plotdata_df <- plotdata_df %>%
  #   filter(Date %in% dates_filtered)
  # ## save plot data
  # write.table(x = plotdata_df, file = paste0("../../plot_data/F2i.", gene_plot, ".SourceData.tsv"), quote = F, sep = "\t", row.names = F)

  ## input plot data
  plotdata_df <- fread(data.table = F, input = paste0("../../plot_data/F2i.", gene_plot, ".SourceData.tsv"))
  ## calculate limit for y-axis
  plotdata_sum_df <- plotdata_df %>%
    group_by(Line) %>%
    summarise(y_plot = mean(Value_bytub_byscrambled),
              sd_plot = sd(Value_bytub_byscrambled))
  ymax <- max(plotdata_sum_df$y_plot+plotdata_sum_df$sd_plot)
  
  # do test and plot -----------------------------------------------------------
  stat.test <- compare_means(
    Value_bytub_byscrambled ~ Line, data = plotdata_df,
    method = test_plot, ref.group = "RCC4_scrambled"
  )
  
  p <- ggbarplot(data = plotdata_df, x = "Line", y = "Value_bytub_byscrambled", 
                 add = c("mean_sd", "dotplot"), fill = "Line", color = "black", 
                 error.plot = "upper_linerange")
  p <- p + stat_pvalue_manual(stat.test, 
                              y.position = seq(ymax*1.05, ymax*(1+0.1*(length(lines_plot)-1)), length.out = (length(lines_plot)-1)), 
                              label = "p = {signif(p, digits = 2)}",
                              label.size = 5)
  p <- p + scale_fill_manual(values = colors_byline)
  p <- p + scale_x_discrete(labels = c("RCC4_scrambled" = "sh-NC", "RCC4_KLF9_C2" = "sh-KLF9"))
  p <- p  + scale_y_continuous(expand = c(0,0), limits = c(0, ymax*(1+0.1*(length(lines_plot))))) ## CP expression for MXI1 lines - 3 lines
  p <- p + theme_classic()
  p <- p + ylab(label = paste0("Relative ", gene_plot, " level"))
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, color = "black", size = 12), 
                 axis.text.y = element_text(color = "black", size = 12))
  p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(color = "black", size = 12))
  
  ## save plot
  file2write <- paste0(dir_out,"F2i_Barplot_WB_densitometry_", gene_plot,  ".pdf")
  pdf(file2write, width = 1.75, height = 2.5, useDingbats = F)
  print(p)
  dev.off()
}


