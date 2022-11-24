#### Initialize ----

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# clear workspace
rm(list = ls())

library(tidyverse)
library(openxlsx)
library(maplet)

#### Load mapping and annotation files ----

# load master mapping file
mapping_file <- read.csv2(file = "data/MasterMapping_MetImmune_03_16_2022_release.csv", sep=",")

# change name of one of the ccRCC2 dataset batches for readability
mapping_file$Dataset[mapping_file$RNAFile == "MultiRegionalRCC.tpm.gene_symbol.csv"] <- "ccRCC3"

mm <- mapping_file %>% dplyr::select(MetabFile,RNAFile,Dataset)  %>% 
  distinct() %>%
  dplyr::arrange(Dataset) %>%
  dplyr::mutate(Cohort=Dataset)

#### Load Metabolomics Data ----

filedir_data <- "results/preprocessed_data"
filelist <- list.files(filedir_data) 
filelist <- filelist %>% 
  grep(pattern = "PreprocessedData_", x = .,value = T) %>%
  {.[!(. %in% grep(pattern = "RNA", x = .,value = T))]}

cohorts <- sub(".xlsx","",sub("^[^_]*_", "", filelist))

# load preprocessed and imputed metabolomics data
met <- lapply(cohorts %>% {names(.)=.;.}, function(x){
  read.xlsx(sprintf("results/preprocessed_data/PreprocessedData_%s.xlsx",x), sheet = "data_imputed", rowNames = T)
})
# add second entry for ccRCC2 (ie. ccRCC3) 
# for this dataset there are two batches of RNA that should be treated separately
met$ccRCC3 <- met$ccRCC2
# reorder in alphabetical order
met <- met[sort(names(met))]

# load sample annotations
anno <- lapply(cohorts %>% {names(.)=.;.}, function(x){
  read.xlsx(sprintf("results/preprocessed_data/PreprocessedData_%s.xlsx",x), sheet = "sampleanno", rowNames = T, )
})
# add second entry for ccRCC2 (ie. ccRCC3) 
# for this dataset there are two batches of RNA that should be treated separately
anno$ccRCC3 <- anno$ccRCC2
# reorder in alphabetical order
anno <- anno[sort(names(anno))]

# add cohort names to mapping file
mapping_file %<>%
  dplyr::left_join(mm %>% dplyr::select(Dataset, Cohort), by="Dataset") 

mapping_file$Cohort %>% table
sapply(lapply(met, colnames), function(x){which(x %in% mapping_file$MetabID) %>% length}) %>% sum()

#### Subset data to only samples with both metabolomics and rna data ----

met <- lapply(names(met) %>% {names(.)=.;.}, function(x) {
  met[[x]][,colnames(met[[x]]) %in% mapping_file$MetabID[mapping_file$Cohort==x]]
})

#### Reorder samples so that they are matched correctly between metabolomics and rna data ----

# reorder samples to match with mapping file
met <- lapply(names(met) %>% {names(.)=.;.}, function(x){
  metorder <- mapping_file$MetabID[mapping_file$Cohort %in% x]
  # print(all.equal(metorder, colnames(met[[x]])[match(metorder, colnames(met[[x]]))]))
  met[[x]][,match(metorder, colnames(met[[x]]))] 
})

# check order
sapply(names(met), function(x){
  all(colnames(met[[x]]) == mapping_file$MetabID[mapping_file$Cohort %in% x])
}) %>% all

#### Check initial dimensions ----

sapply(met, ncol)

#### Separate Tumor and Normal Samples ----

met_T <- lapply(met, function(x) {
  x[,colnames(x) %in% mapping_file$MetabID[mapping_file$TN %in% c("Tumor", "TUMOR")]]
})
met_N <- lapply(met, function(x) {
  x[,colnames(x) %in% mapping_file$MetabID[mapping_file$TN %in% c("Normal", "NORMAL")]]
})

# remove empty elements from list of normal samples
met_N <- met_N[sapply(met_N, function(x) dim(x)[2]) > 0]

# check dimensions
sapply(met_T, dim)
sapply(met_N, dim)

#### Filter ----

# first remove variables with more than imp.max missing values
imp.max <- 0.8

countmiss <- function(x, y, imp.max){
  # find min
  m <- min(y, na.rm=T)
  # count number of equal values
  repeats <- which(y==m) %>% length 
  # find imputation percentage
  (repeats-1)/ncol(x)
}

filter_func <- function(x, imp.max){
  # compute variance
  zerovar <- apply(x, 1, var)
  # remove metabolites with zero variance
  x <- x[which(zerovar!=0 & !is.na(zerovar)),]
  # compute missingness percentage
  miss <- apply(x, 1, function(y){
    countmiss(x, y, imp.max)
  })
  # only keep metabolites with less or equal imp.max missing values
  x[which(miss<=imp.max),]
}

# filter
met_T <- lapply(met_T, function(x){filter_func(x, imp.max)})
met_N <- lapply(met_N, function(x){filter_func(x, imp.max)})

# check dimensions
sapply(met_T, dim)
sapply(met_N, dim)

#### Save filtered Data to File ----

# concatenate Tumor and Normal lists
met_all <- c(met_T,met_N)
names(met_all)[1:length(met_T)] <- sprintf("%s_Tumor",names(met_T))
names(met_all)[(length(met_T)+1):length(met_all)] <- sprintf("%s_Normal",names(met_N))

# create file for ccRCC3
wb <- loadWorkbook(file="results/preprocessed_data/PreprocessedData_ccRCC2.xlsx")
saveWorkbook(wb, file="results/preprocessed_data/PreprocessedData_ccRCC3.xlsx", overwrite = TRUE)

# add one sheet to preprocessed data files
lapply(names(met_all), function(x){
  print(x)
  # load Excel file
  wb <- loadWorkbook(file=sprintf("results/preprocessed_data/PreprocessedData_%s.xlsx",gsub(pattern = "_Tumor|_Normal", replacement = "", x)))
  # add sheet
  sheet = addWorksheet(wb, sprintf("metabo_imputed_filtered_%s",strsplit(x, "_(?!.*_)", perl=TRUE)[[1]][2]))
  writeData(wb, sheet=sheet, met_all[[x]], rowNames = T, colNames = T)
  # save workbook
  saveWorkbook(wb, sprintf("results/preprocessed_data/PreprocessedData_%s.xlsx",gsub(pattern = "_Tumor|_Normal", replacement = "", x)), overwrite = TRUE)
}) %>% invisible
  
#### Save Results ----

save(met_all, anno, mapping_file, mm,
     file = "results/Workspace_3_FilterMetabo.Rdata")
 