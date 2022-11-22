#### Initialize ----

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(openxlsx)
library(maplet)

#### Load mapping and annotation files ----

# load preprocessed metabolomics data
load("results/Workspace_3_FilterMetabo.Rdata")

#### Load Transcriptomics Data ----

# load rna data
rna <- lapply(1:length(mm$Cohort), function(x) {
  read.table(file=paste0("data/transcriptomics_processed/", mm$RNAFile[x], sep=""), sep = ",", row.names = 1, header = TRUE)
})
names(rna) <- mm$Cohort

#### Subset data to only samples with both metabolomics and rna data ----

rna <- lapply(names(rna) %>% {names(.)=.;.}, function(x) {
  rna[[x]][,colnames(rna[[x]]) %in% (mapping_file$RNAID[mapping_file$Cohort==x] %>% make.names)]
})

#### Reorder samples so that they are matched correctly between metabolomics and rna data ----

# reorder samples to match with mapping file
rna <- lapply(names(rna) %>% {names(.)=.;.}, function(x){
  rnaorder <- mapping_file$RNAID[mapping_file$Cohort %in% x] %>% make.names()
  rna[[x]][,match(rnaorder, colnames(rna[[x]]))] 
})

# check order
sapply(names(rna), function(x){
  all(colnames(rna[[x]]) == mapping_file$RNAID[mapping_file$Cohort %in% x] %>% make.names()) 
}) %>% all

#### Check initial dimensions ----

sapply(rna,dim)

#### Separate Tumor and Normal Samples ----

rna_T <- lapply(rna, function(x) {
  x[,colnames(x) %in% (mapping_file$RNAID[mapping_file$TN %in% c("Tumor", "TUMOR")] %>% make.names)]
})
rna_N <- lapply(rna, function(x) {
  x[,colnames(x) %in% (mapping_file$RNAID[mapping_file$TN %in% c("Normal", "NORMAL")] %>% make.names)]
})

# remove empty elements from list of normal samples
rna_N <- rna_N[sapply(rna_N, function(x) dim(x)[2]) > 0]

# check dimensions
sapply(rna_T, dim)
sapply(rna_N, dim)

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
  x <- x[which(zerovar!=0),]
  # compute missingness percentage
  miss <- apply(x, 1, function(y){
    countmiss(x, y, imp.max)
  })
  # only keep metabolites with less or equal imp.max missing values
  x[which(miss<=imp.max),]
}

# filter
rna_T <- lapply(rna_T, function(x){filter_func(x, imp.max)})
rna_N <- lapply(rna_N, function(x){filter_func(x, imp.max)})

# check dimensions
sapply(rna_T, dim)
sapply(rna_N, dim)

#### Save filtered Data to File ----

rna_all <- c(rna_T,rna_N)
names(rna_all)[1:length(rna_T)] <- sprintf("%s_Tumor",names(rna_T))
names(rna_all)[(length(rna_T)+1):length(rna_all)] <- sprintf("%s_Normal",names(rna_N))

# add one sheet to preprocessed data files
lapply(names(rna_all), function(x){
  print(x)
  # load Excel file
  wb <- loadWorkbook(file=sprintf("data/preprocessed/PreprocessedData_%s.xlsx",gsub(pattern = "_Tumor|_Normal", replacement = "", x)))
  # add sheet
  sheet = addWorksheet(wb, sprintf("rna_imputed_filtered_%s",strsplit(x, "_(?!.*_)", perl=TRUE)[[1]][2]))
  writeData(wb, sheet=sheet, rna_all[[x]], rowNames = T, colNames = T)
  # save workbook
  saveWorkbook(wb, sprintf("data/preprocessed/PreprocessedData_%s.xlsx",gsub(pattern = "_Tumor|_Normal", replacement = "", x)), overwrite = TRUE)
}) %>% invisible

#### Save Results ----

save(rna_all,
     file = "results/Workspace_4_FilterRNA.Rdata")
