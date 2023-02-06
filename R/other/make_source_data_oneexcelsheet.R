# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  make source data excel sheet

#=======================================================================================
# Load libraries  and set working directory  -----------------------------------
## load libraries
packages = c(
  "data.table",
  "stringr",
  "plyr",
  "dplyr",
  "xlsx",
  "openxlsx"
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

# get all file names ------------------------------------------------------
library(xlsx)
# xlsx_filename <- "../../../202301_Nat_Comm/Source_Data/Source_Data.nonWB.020223.xlsx"
xlsx_filename <- "~/Documents/Project/ccRCC_snRNA/Source_Data/Source_Data.nonWB.020423.xlsx"
dir_sourcedatafiles <- "~/Documents/Project/ccRCC_snRNA/Source_Data/Individual_Source_Data/"
# filenames_process <- list.files(path = "../../plot_data/")
filenames_process <- list.files(path = dir_sourcedatafiles)
# Create a blank workbook
OUT <- createWorkbook()
for (filename_tmp in filenames_process) {
  dataframe_tmp <- fread(data.table = F, input = paste0(dir_sourcedatafiles, filename_tmp))
  sheetname_tmp <- gsub(x = filename_tmp, pattern = "\\.SourceData\\.tsv", replacement = "")
  
  # Add some sheets to the workbook
  addWorksheet(OUT, sheetname_tmp)
  
  # Write the data to the sheets
  writeData(OUT, sheet = sheetname_tmp, x = dataframe_tmp)
  
  # write.xlsx(dataframe_tmp, file = xlsx_filename, sheetName = sheetname_tmp, append=TRUE, row.names=FALSE)
}
# Export the file
saveWorkbook(OUT, xlsx_filename)