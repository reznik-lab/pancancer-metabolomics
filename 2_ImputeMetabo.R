#### Initialize ----

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# clear workspace
rm(list = ls())

library(tidyverse)
library(openxlsx)
library(maplet)

#### Define Functions -----

# write SE to xls
writeSE2xls <- function(D, name) {
  
  wb = createWorkbook()
  sheet = addWorksheet(wb, "data")
  writeData(wb, sheet=sheet, as.data.frame(assay(D)), rowNames = T, colNames = T)
  sheet = addWorksheet(wb, "metanno")
  writeData(wb, sheet=sheet, as.data.frame(rowData(D)), colNames = T)
  sheet = addWorksheet(wb, "sampleanno")
  writeData(wb, sheet=sheet, as.data.frame(colData(D)), colNames = T)
  saveWorkbook(wb, paste0(name,".xlsx", sep=""), overwrite = TRUE)
}

#### Load normalized data from file ----

load("results/Workspace_1_PreprocessMetabo.Rdata")

#### Impute missing entries ----

# impute missing entries with minimum value per metabolite
lapply(names(D), function(x){
  print(x)
  D[[x]] <<- D[[x]] %>%
    mt_pre_impute_min()
}) %>% invisible()

#### Save Imputed data to file ----

# add one sheet to preprocessed data files
lapply(names(D), function(x){
  print(x)
  # load Excel file
  wb <- loadWorkbook(file=sprintf("results/preprocessed_data/PreprocessedData_%s.xlsx",x))
  # add sheet
  sheet = addWorksheet(wb, "data_imputed")
  writeData(wb, sheet=sheet, as.data.frame(assay(D[[x]])), rowNames = T, colNames = T)
  # save workbook
  saveWorkbook(wb, sprintf("results/preprocessed_data/PreprocessedData_%s.xlsx",x), overwrite = TRUE)
}) %>% invisible

#### Save Workspace to file ----

save(D, metanno, mapping_file, 
     file = "results/Workspace_2_ImputeMetabo.Rdata")
