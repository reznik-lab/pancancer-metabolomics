#### Initialize ---

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# clear workspace
rm(list = ls())

# load libraries
library(readxl)
library(data.table)
library(dplyr)
library(stringr)
library(dplyr)
library(purrr)
library(survival)
library(parallel)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(patchwork)

#### Load and clean data ----

#Load master mapping spreadsheet
master_mapping <- read.csv("data/MasterMapping_MetImmune_03_16_2022_release.csv")

#Load metabolomics data
met_fs = list.files(path='data/metabolomics_processed', full.names = TRUE)
met_data_tumor = lapply(met_fs, function(x) as.data.frame(read_excel(x, sheet = 5)))
tumor_names <- gsub("PreprocessedData_","",list.files(path='data/metabolomics_processed')) #remove everything before the dataset name
tumor_names <- gsub(".xlsx","",tumor_names) #remove the .xlsx extension
tumor_names <- paste(tumor_names,"Tumor",sep="_") #add tumor to end so we know what type of sample it is
names(met_data_tumor) <- tumor_names

met_data_normal = lapply(met_fs[c(1,5,6,7,9,11,14,15)], function(x) as.data.frame(read_excel(x, sheet = 6))) #these are the datasets with normal samples
normal_names <- gsub("PreprocessedData_","",list.files(path='data/metabolomics_processed')) #remove everything before the dataset name
normal_names <- gsub(".xlsx","",normal_names) #remove the .xlsx extension
normal_names <- paste(normal_names,"Normal",sep="_") #add tumor to end so we know what type of sample it is
names(met_data_normal) <- normal_names[c(1,5,6,7,9,11,14,15)]

# concatenate tumor and normal
met_data <- c(met_data_tumor,met_data_normal)
# order alphabetically
met_data <- met_data[order(names(met_data))]

# Make the first column the rownames
for (i in 1:23) {
  rownames(met_data[[i]]) <- met_data[[i]][,1] #set the first column (metabolites) as the rownames
  met_data[[i]] <- met_data[[i]][,-1] #remove the first column
}

# Check there are 988 samples
sum(sapply(1:23,function(i){ncol(met_data[[i]])}))

# Load metabolite information
metabolite_mapping = as.data.frame(read_excel('data/MetinfoFromData.xlsx', sheet = 2))

# Load TME data
tme_fs = list.files(path='data/TME_deconvolution_processed', full.names = TRUE)
tme_data = lapply(tme_fs, function(x) read.csv(file=x,row.names = 1, check.names = FALSE))
names(tme_data) <- list.files(path='data/TME_deconvolution_processed/')

# Harmonize TME data with metabolomics data
tme_data$BRCA2_Tumor <- tme_data$MergedDeconvolution.PMC4187326_BRCA.csv[master_mapping$ITHID[match(colnames(met_data$BRCA2_Tumor),master_mapping$MetabID)],]
tme_data$BRCA1_Normal <- tme_data$GSE37751.ssGSEA.UnNorm.Deconvolution.csv[master_mapping$ITHID[match(colnames(met_data$BRCA1_Normal),master_mapping$MetabID)],]
tme_data$BRCA1_Tumor <- tme_data$GSE37751.ssGSEA.UnNorm.Deconvolution.csv[master_mapping$ITHID[match(colnames(met_data$BRCA1_Tumor),master_mapping$MetabID)],]
tme_data$COAD_Normal <- tme_data$MergedDeconvolution.GSE89076.csv[master_mapping$ITHID[match(colnames(met_data$COAD_Normal),master_mapping$MetabID)],]
tme_data$COAD_Tumor <- tme_data$MergedDeconvolution.GSE89076.csv[master_mapping$ITHID[match(colnames(met_data$COAD_Tumor),master_mapping$MetabID)],]
tme_data$DLBCL_Tumor <- tme_data$MergedDeconvolution.Cornell_DLBCL.csv[master_mapping$ITHID[match(colnames(met_data$DLBCL_Tumor),master_mapping$MetabID)],]  
tme_data$GBM_Normal <- tme_data$MergedDeconvolution.CPTAC_GBM.csv[master_mapping$ITHID[match(colnames(met_data$GBM_Normal),master_mapping$MetabID)],] 
tme_data$GBM_Tumor <- tme_data$MergedDeconvolution.CPTAC_GBM.csv[master_mapping$ITHID[match(colnames(met_data$GBM_Tumor),master_mapping$MetabID)],] 
tme_data$HurthleCC_Normal <- tme_data$HCC.ssGSEA.UnNorm.Deconvolution.csv[master_mapping$ITHID[match(colnames(met_data$HurthleCC_Normal),master_mapping$MetabID)],]
tme_data$HurthleCC_Tumor <- tme_data$HCC.ssGSEA.UnNorm.Deconvolution.csv[master_mapping$ITHID[match(colnames(met_data$HurthleCC_Tumor),master_mapping$MetabID)],]
tme_data$HCC_Tumor <- tme_data$MergedDeconvolution.GSE76297.csv[master_mapping$ITHID[match(colnames(met_data$HCC_Tumor),master_mapping$MetabID)],]
tme_data$ICC_Tumor <- tme_data$MergedDeconvolution.GSE76297.csv[master_mapping$ITHID[match(colnames(met_data$ICC_Tumor),master_mapping$MetabID)],]
tme_data$OV_Tumor <- tme_data$MergedDeconvolution.GSE26193.csv[master_mapping$ITHID[match(colnames(met_data$OV_Tumor),master_mapping$MetabID)],]
tme_data$PDAC_Normal <- tme_data$MergedDeconvolution.GSE62452.csv[master_mapping$ITHID[match(colnames(met_data$PDAC_Normal),master_mapping$MetabID)],]
tme_data$PDAC_Tumor <- tme_data$MergedDeconvolution.GSE62452.csv[master_mapping$ITHID[match(colnames(met_data$PDAC_Tumor),master_mapping$MetabID)],]
tme_data$PRAD_Normal <- tme_data$MergedDeconvolution.Cornell_PROSTATE.csv[master_mapping$ITHID[match(colnames(met_data$PRAD_Normal),master_mapping$MetabID)],]
tme_data$PRAD_Tumor <- tme_data$MergedDeconvolution.Cornell_PROSTATE.csv[master_mapping$ITHID[match(colnames(met_data$PRAD_Tumor),master_mapping$MetabID)],]
tme_data$ccRCC1_Tumor <- tme_data$MergedDeconvolution.Project_08402_M.csv[master_mapping$ITHID[match(colnames(met_data$ccRCC1_Tumor),master_mapping$MetabID)],]
tme_data$ccRCC2_Tumor <- tme_data$MergedDeconvolution.Flow4Batch.csv[master_mapping$ITHID[match(colnames(met_data$ccRCC2_Tumor),master_mapping$MetabID)],]
tme_data$ccRCC3_Normal <- tme_data$Multiregions.2Batchs.ImmuneDeconvolution.csv[master_mapping$ITHID[match(colnames(met_data$ccRCC3_Normal),master_mapping$MetabID)],]
tme_data$ccRCC3_Tumor <- tme_data$Multiregions.2Batchs.ImmuneDeconvolution.csv[master_mapping$ITHID[match(colnames(met_data$ccRCC3_Tumor),master_mapping$MetabID)],]
tme_data$ccRCC4_Normal <- tme_data$Multiregions.2Batchs.ImmuneDeconvolution.csv[master_mapping$ITHID[match(colnames(met_data$ccRCC4_Normal),master_mapping$MetabID)],]
tme_data$ccRCC4_Tumor <- tme_data$Multiregions.2Batchs.ImmuneDeconvolution.csv[master_mapping$ITHID[match(colnames(met_data$ccRCC4_Tumor),master_mapping$MetabID)],]

# remove unnecessary data
tme_data <- tme_data[-c(1:13)] 
# order alphabetically
tme_data <- tme_data[order(names(tme_data))]

# verify that samples in metabolite and TME data are in same order 
sapply(1:23, function(i){all(master_mapping$CommonID[match(rownames(tme_data[[i]]),master_mapping$ITHID)]==
                               master_mapping$CommonID[match(colnames(met_data[[i]]),master_mapping$MetabID)])}) %>% all

# Check there are 988 samples
sum(sapply(1:23,function(i){nrow(tme_data[[i]])}))

# set list of signatures required
bindea_sigs <- c("Th17 cells","NK cells","Tcm cells","T helper cells","Eosinophils","pDC","Tem cells","NK CD56dim cells",
                 "Neutrophils","Macrophages","DC","iDC","NK CD56bright cells","Tfh cells","Treg cells","aDC","Th1 cells",
                 "CD8 T cells","Cytotoxic cells","T cells","B cells","Mast cells","Th2 cells")
imp_signatures <- c("ImmuneScore",bindea_sigs)

# Clean TME data so we're only including signatures that we use in our analyses
for (i in 1:23) {
  tme_data[[i]] <- tme_data[[i]][,imp_signatures] 
}

# Load RNA data (only need datasets that measured histamine, GBM and OV don't have histamine measured)
rna_fs = list.files(path='data/transcriptomics_processed', full.names = TRUE)
rna_data = lapply(rna_fs, function(x) read.csv(file=x,row.names = 1, check.names = FALSE))
names(rna_data) <- list.files(path='data/transcriptomics_processed/')
for(i in 1:23){
  #include samples that have RNA, TME and metabolomics data
  rna_data[[i]] <- rna_data[[i]][,sapply(1:23, function(ii){which(colnames(rna_data[[ii]]) %in% master_mapping$RNAID)})[[i]]] 
}

rna_data$BRCA2_Tumor <- rna_data$PMC4187326_BRCA.tpm.gene_symbol.csv[,master_mapping$RNAID[match(colnames(met_data$BRCA2_Tumor),master_mapping$MetabID)]]
rna_data$BRCA1_Tumor <- rna_data$GSE37751.hugene10st.gene_symbol.csv[,master_mapping$RNAID[match(colnames(met_data$BRCA1_Tumor),master_mapping$MetabID)]]
rna_data$COAD_Tumor <- rna_data$GSE89076.Agilent8x60K.log2_transformed.gene_symbol.csv[,master_mapping$RNAID[match(colnames(met_data$COAD_Tumor),master_mapping$MetabID)]]
rna_data$DLBCL_Tumor <- rna_data$Cornell_DLBCL.tpm.gene_symbol.csv[,master_mapping$RNAID[match(colnames(met_data$DLBCL_Tumor),master_mapping$MetabID)]]
rna_data$HurthleCC_Tumor <- rna_data$HCC.tpm.gene_symbol.csv[,master_mapping$RNAID[match(colnames(met_data$HurthleCC_Tumor),master_mapping$MetabID)]]
rna_data$HCC_Tumor <- rna_data$GSE76297.hugene20st.gene_symbol.csv[,master_mapping$RNAID[match(colnames(met_data$HCC_Tumor),master_mapping$MetabID)]]
rna_data$ICC_Tumor <- rna_data$GSE76297.hugene20st.gene_symbol.csv[,master_mapping$RNAID[match(colnames(met_data$ICC_Tumor),master_mapping$MetabID)]]
rna_data$PDAC_Tumor <- rna_data$GSE62452.hugene10st.gene_symbol.csv[,master_mapping$RNAID[match(colnames(met_data$PDAC_Tumor),master_mapping$MetabID)]]
rna_data$PRAD_Tumor <- rna_data$Cornell_PROSTATE.BatchAdj.hugene20st.gene_symbol.csv[,master_mapping$RNAID[match(colnames(met_data$PRAD_Tumor),master_mapping$MetabID)]]
rna_data$ccRCC1_Tumor <- rna_data$Project_08402_M.tpm.gene_symbol.csv[,master_mapping$RNAID[match(colnames(met_data$ccRCC1_Tumor),master_mapping$MetabID)]]
rna_data$ccRCC2_Tumor <- rna_data$Flow4Batch.tpm.gene_symbol.csv[,master_mapping$RNAID[match(colnames(met_data$ccRCC2_Tumor),master_mapping$MetabID)]]
rna_data$ccRCC3_Tumor <- rna_data$MultiRegionalRCC.tpm.gene_symbol.csv[,master_mapping$RNAID[match(colnames(met_data$ccRCC3_Tumor),master_mapping$MetabID)]]
rna_data$ccRCC4_Tumor <- rna_data$Multiregions.2Batchs.Sample.hg19KnownGene.tpm.gene_symbol.csv[,master_mapping$RNAID[match(colnames(met_data$ccRCC4_Tumor),master_mapping$MetabID)]]

# remove unnecessary data
rna_data <- rna_data[-c(1:23)]
# order alphabetically
rna_data <- rna_data[order(names(rna_data))]

# Load aDC signature without IDO1
adc_fs = list.files(path='data/Metabolism_Immune.aDC_exIDO1', full.names = TRUE)
adc_data = lapply(adc_fs, function(x) read.csv(file=x,row.names = 1, check.names = FALSE))
names(adc_data) <- list.files(path='data/Metabolism_Immune.aDC_exIDO1')
rownames(adc_data$CPTAC_GBM.aDC_exIDO1.csv) <- gsub("\\.","-",rownames(adc_data$CPTAC_GBM.aDC_exIDO1.csv))
rownames(adc_data$Project_08402_M.aDC_exIDO1.csv) <- gsub("\\.","-",rownames(adc_data$Project_08402_M.aDC_exIDO1.csv))
adc_data <- adc_data[-6]

adc_data_sub <- adc_data

# select samples that have RNA, TME and metabolomics data
for(i in 1:13){
  adc_data_sub[[i]] <- adc_data_sub[[i]][sapply(1:13, function(ii){which(rownames(adc_data_sub[[ii]]) %in% master_mapping$ITHID)})[[i]],] 
}

for(i in 1:13){
  adc_data_sub[[i]] <- as.data.frame(adc_data_sub[[i]])
}

for(i in 1:13){
  rownames(adc_data_sub[[i]]) <-  sapply(1:13,function(ii){rownames(adc_data[[ii]])[which(rownames(adc_data[[ii]]) %in% master_mapping$ITHID)]})[[i]]
}

adc_data_sub$BRCA2_Tumor <- adc_data_sub$PMC4187326_BRCA.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$BRCA2_Tumor),master_mapping$MetabID)],]
adc_data_sub$BRCA1_Tumor <- adc_data_sub$GSE37751.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$BRCA1_Tumor),master_mapping$MetabID)],]
adc_data_sub$COAD_Tumor <- adc_data_sub$GSE89076.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$COAD_Tumor),master_mapping$MetabID)],]
adc_data_sub$DLBCL_Tumor <- adc_data_sub$Cornell_DLBCL.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$DLBCL_Tumor),master_mapping$MetabID)],]  
adc_data_sub$GBM_Tumor <- adc_data_sub$CPTAC_GBM.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$GBM_Tumor),master_mapping$MetabID)],] 
adc_data_sub$HurthleCC_Tumor <- adc_data_sub$HCC.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$HurthleCC_Tumor),master_mapping$MetabID)],]
adc_data_sub$HCC_Tumor <- adc_data_sub$GSE76297.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$HCC_Tumor),master_mapping$MetabID)],]
adc_data_sub$ICC_Tumor <- adc_data_sub$GSE76297.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$ICC_Tumor),master_mapping$MetabID)],]
adc_data_sub$OV_Tumor <- adc_data_sub$GSE26193.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$OV_Tumor),master_mapping$MetabID)],]
adc_data_sub$PDAC_Tumor <- adc_data_sub$GSE62452.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$PDAC_Tumor),master_mapping$MetabID)],]
adc_data_sub$PRAD_Tumor <- adc_data_sub$Cornell_PROSTATE.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$PRAD_Tumor),master_mapping$MetabID)],]
adc_data_sub$ccRCC1_Tumor <- adc_data_sub$Project_08402_M.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$ccRCC1_Tumor),master_mapping$MetabID)],]
adc_data_sub$ccRCC2_Tumor <- adc_data_sub$Flow4Batch.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$ccRCC2_Tumor),master_mapping$MetabID)],]
adc_data_sub$ccRCC3_Tumor <- adc_data_sub$MultiRegionalRCC.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$ccRCC3_Tumor),master_mapping$MetabID)],]
adc_data_sub$ccRCC4_Tumor <- adc_data_sub$MultiRegionalRCC.aDC_exIDO1.csv[master_mapping$ITHID[match(colnames(met_data$ccRCC4_Tumor),master_mapping$MetabID)],]

# remove unnecessary datasets
adc_data_sub <- adc_data_sub[-c(1:13)]
# order alphabetically
adc_data_sub <- adc_data_sub[order(names(adc_data_sub))]

for(i in 1:15){
  adc_data_sub[[i]] <- scale(adc_data_sub[[i]])
}

# Load flow-sorted ovarian data from Kilgour et al. 2021
kilgour_data_raw <- read.csv("data/flow_sorted_ovarian_metabolomics/Kilgour_2021_Raw_Data.csv")
rownames(kilgour_data_raw) <- kilgour_data_raw[,1]
kilgour_data_raw <- kilgour_data_raw[,-1]

# Clean flow-sorted ovarian data
# Probabilistic Quotient Normalization of the raw data
normalizationPQN<-function(dat,group="all",refGroupName){
  # row: metabolite
  # column: sample
  
  # step 1: 
  # get a reference sample aboundance by each metabolite median across cohort
  # for each row (metabolite), calculate median number
  if(group =="all"){
    refSampleVector<-apply(dat,1,function(y) {median(y, na.rm=TRUE)})
  }else{
    subsetDat<-dat[,colnames(dat) %in% refGroupName]
    refSampleVector<-apply(subsetDat,1,function(y) {median(y, na.rm=TRUE)})
  }
  
  # step2:
  # 
  quotientMatrix<-apply(dat,2, function(y){y/refSampleVector})
  
  
  # step3:
  #
  dilutionFactor<-apply(quotientMatrix,2,function(y){median(y, na.rm=TRUE)})
  
  
  # step4:
  # normalization
  
  newDat<-apply(dat,1,function(y){y/dilutionFactor})
  newDat<-t(newDat)
  
  return(newDat)
  
}


kilgour_data_pqn <- normalizationPQN(kilgour_data_raw)
kilgour_data <- log2(kilgour_data_pqn) #log2 normalize data 

#Colors for legends
mapping_color <- data.frame(Name = c("BRCA1","BRCA2","COAD","GBM","DLBCL",
                                     "HurthleCC","HCC","ICC","OV","PDAC",
                                     "PRAD","ccRCC1","ccRCC2","ccRCC3","ccRCC4")) %>%
  dplyr::arrange(Name) %>%
  dplyr::mutate(Color = c("#e5b6a5","#c3c1a4","#e5e3c0","#bbe1b9","#88beab",
                          "#97b1ab","#abc5bf","#abe7e1","#c6e6e7","#8fcde1",
                          "#88aee1","#a6b7d3","#d2d5ed","#c8b7de","#e6bbcc"))




#### Run concordance analysis ----
#Z-Score signatures
tme_zscored <- lapply(1:23,function(i){apply(tme_data[[i]],2,scale)}) #Z-score
tme_zscored_t <- lapply(1:23,function(i){t(tme_zscored[[i]])}) #transpose so it's in the same format as the metabolomics data
tme_zscored_t <- lapply(1:23,function(i){as.data.frame(tme_zscored_t[[i]])})
names(tme_zscored_t) <- names(tme_data) #add names to data
for(i in 1:23){
  colnames(tme_zscored_t[[i]]) <- rownames(tme_data[[i]]) #add column names back
}

#Extract tumor samples
tumors <- which(str_detect(names(tme_data),"Tumor")) #which datasets are tumor datasets
met_data_tumors <- met_data[tumors] #running analysis only on tumor samples
tme_zscored_t_tumors <- tme_zscored_t[tumors] #running analysis only on tumor samples

#Build matricies for concordance analysis
joint_metabolite_df_t <- lapply(met_data_tumors, function(y){
  y %>% tibble::rownames_to_column("metabolite")
}) %>%
  purrr::reduce(full_join, by="metabolite") %>%
  tibble::column_to_rownames("metabolite")

joint_imm_sig_df_t <- lapply(tme_zscored_t_tumors, function(y){
  y %>% tibble::rownames_to_column("imm_sig")
}) %>%
  purrr::reduce(full_join, by="imm_sig") %>%
  tibble::column_to_rownames("imm_sig")

# Sample size of all cohorts
ns_t <- sapply(1:15,function(i){ncol(met_data_tumors[[i]])})
nsum_t <- cumsum(ns_t)
#Vector with dataset id for each sample in the joint dataframe
dataset_id_t <- lapply(seq(ns_t), function(i) rep(i, ns_t[i])) %>% unlist
#Weight vector for each sample in the joint dataframe
weights_t <- lapply(seq(ns_t), function(i) rep(1/ns_t[i], ns_t[i])) %>% unlist

#Concordance meta-analysis (all datasets)
ns_t <- sapply(1:15,function(i){ncol(met_data_tumors[[i]])})
# vector with dataset id for each sample in the joint dataframe
dataset_id_t <- lapply(seq(ns_t), function(i) rep(i, ns_t[i])) %>% unlist
# weight vector for each sample in the joint dataframe
weights_t <- lapply(seq(ns_t), function(i) rep(1/ns_t[i], ns_t[i])) %>% unlist

tumor_concordance <- mclapply(rownames(joint_metabolite_df_t)[1:2359], function(m){
  mclapply(rownames(joint_imm_sig_df_t)[1:24], function(g){ 
    
    # immune signature value
    imm_signature <- joint_imm_sig_df_t[g,] %>% unlist
    # metabolite values
    met <- joint_metabolite_df_t[m,] %>% unlist
    
    # remove samples with NA metabolite values
    # (NAs occur because not all datasets have the same metabolites)
    imm_sig_full <- imm_signature[!is.na(met)]
    met_full <- met[!is.na(met)]
    dataset_id_full <- dataset_id_t[!is.na(met)]
    weights_full <- weights_t[!is.na(met)]
    
    # compute concordance
    conc <- survival::concordance(met_full ~ imm_sig_full + strata(dataset_id_full), keepstrata = T, weights=weights_full)
    
    # extract concordance of each dataset
    if (class(conc$count)!="numeric"){
      individual_conc <- conc$count %>% {.[,1]/rowSums(.[,1:3])} %>% data.frame %>% t %>% data.frame
      colnames(individual_conc) <- names(met_data_tumors)[unique(dataset_id_full)]
    } else {
      individual_conc <- conc$count %>% as.matrix() %>% t %>% data.frame %>%
        {.[,1]/rowSums(.[,1:3])} %>% data.frame %>% t %>% data.frame
      colnames(individual_conc) <- names(met_data_tumors)[unique(dataset_id_full)]
    }
    
    X <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                  names(met_data_tumors)) %>%
      dplyr::full_join(individual_conc, by=colnames(individual_conc))
    
    # summarize results
    data.frame(metabolite = m,
               immune_signature = g,
               # number of datasets the metabolite was measured in
               n_dataset = dataset_id_full %>% unique %>% length,
               concordance = conc$concordance,
               variance = conc$var,
               zscore = abs(conc$concordance - 0.5)/sqrt(conc$var),
               p.value = 2*pnorm(-abs(abs(conc$concordance - 0.5)/sqrt(conc$var)))) %>%
      cbind(X)
    
  }, mc.cores = 1) %>% 
    # combine all results
    {do.call(rbind,.)}
}, mc.cores = detectCores(logical = FALSE)-1) %>% 
  # combine all results
  {do.call(rbind,.)} 

tumor_concordance$pair <- paste(tumor_concordance$metabolite,tumor_concordance$immune_signature, sep = "-")
tumor_concordance[,c(4,8:22)] <- (2*tumor_concordance[,c(4,8:22)])-1 #scale concordance between -1 and 1

#change names of columns to remove the "_Tumor"
concordance_old_names <- colnames(tumor_concordance)[8:22]
concordance_new_names <- word(concordance_old_names,1,sep = "\\_Tumor")
colnames(tumor_concordance)[8:22] <- concordance_new_names


#Immune Score concordance analysis, with p-values for each dataset
####only for metabolites in 8+ datasets for immune score
metabolites_8 <- unique(tumor_concordance[tumor_concordance$n_dataset >7,1]) #only include metabolites measured in 8+ datasets

#This function returns p-values
concordance_by_type8_is <- function(columns){
  res_T_tumors <- mclapply(metabolites_8, function(m){
    mclapply(rownames(joint_imm_sig_df_t)[1], function(g){ #1 because ImmuneScore is the first signature in the dataframe
      
      # immune score value
      imm_score <- joint_imm_sig_df_t[g,columns] %>% unlist
      # metabolite values
      met <- joint_metabolite_df_t[m,columns] %>% unlist
      
      # remove samples with NA metabolite values
      # (NAs occur because not all datasets have the same metabolites)
      imm_score_full <- imm_score[!is.na(met)]
      met_full <- met[!is.na(met)]
      dataset_id_full <- dataset_id_t[columns][!is.na(met)]
      weights_full <- weights_t[columns][!is.na(met)]
      
      # compute concordance
      conc <- tryCatch(survival::concordance(met_full ~ imm_score_full + strata(dataset_id_full), keepstrata = T, weights=weights_full),
                       error = function(e) return(NA))
      p.value <- tryCatch(2*pnorm(-abs(abs(conc$concordance - 0.5)/sqrt(conc$var))),error = function(e) return(NA))
      
    })})
}

#Calculate p-values for each dataset
BRCA2_columns <- match(colnames(met_data$BRCA2_Tumor),colnames(joint_metabolite_df_t))
BRCA2_concordance8 <- unlist(concordance_by_type8_is(BRCA2_columns))

BRCA1_columns <- match(colnames(met_data$BRCA1_Tumor),colnames(joint_metabolite_df_t))
BRCA1_concordance8 <- unlist(concordance_by_type8_is(BRCA1_columns))

COAD_columns <- match(colnames(met_data$COAD_Tumor),colnames(joint_metabolite_df_t))
COAD_concordance8 <- unlist(concordance_by_type8_is(COAD_columns)) 

DLBCL_columns <- match(colnames(met_data$DLBCL_Tumor),colnames(joint_metabolite_df_t))
DLBCL_concordance8 <- unlist(concordance_by_type8_is(DLBCL_columns)) 

GBM_columns <- match(colnames(met_data$GBM_Tumor),colnames(joint_metabolite_df_t))
GBM_concordance8 <- unlist(concordance_by_type8_is(GBM_columns)) 

HurthleCC_columns <- match(colnames(met_data$HurthleCC_Tumor),colnames(joint_metabolite_df_t))
HurthleCC_concordance8 <- unlist(concordance_by_type8_is(HurthleCC_columns)) 

HCC_columns <- match(colnames(met_data$HCC_Tumor),colnames(joint_metabolite_df_t))
HCC_concordance8 <- unlist(concordance_by_type8_is(HCC_columns)) 

ICC_columns <- match(colnames(met_data$ICC_Tumor),colnames(joint_metabolite_df_t))
ICC_concordance8 <- unlist(concordance_by_type8_is(ICC_columns)) 

OV_columns <- match(colnames(met_data$OV_Tumor),colnames(joint_metabolite_df_t))
OV_concordance8 <- unlist(concordance_by_type8_is(OV_columns)) 

PDAC_columns <- match(colnames(met_data$PDAC_Tumor),colnames(joint_metabolite_df_t))
PDAC_concordance8 <- unlist(concordance_by_type8_is(PDAC_columns))

PRAD_columns <- match(colnames(met_data$PRAD_Tumor),colnames(joint_metabolite_df_t))
PRAD_concordance8 <- unlist(concordance_by_type8_is(PRAD_columns)) 

ccRCC1_columns <- match(colnames(met_data$ccRCC1_Tumor),colnames(joint_metabolite_df_t))
ccRCC1_concordance8 <- unlist(concordance_by_type8_is(ccRCC1_columns)) 

ccRCC2_columns <- match(colnames(met_data$ccRCC2_Tumor),colnames(joint_metabolite_df_t))
ccRCC2_concordance8 <- unlist(concordance_by_type8_is(ccRCC2_columns)) 

ccRCC3_columns <- match(colnames(met_data$ccRCC3_Tumor),colnames(joint_metabolite_df_t))
ccRCC3_concordance8 <- unlist(concordance_by_type8_is(ccRCC3_columns)) 

ccRCC4_columns <- match(colnames(met_data$ccRCC4_Tumor),colnames(joint_metabolite_df_t))
ccRCC4_concordance8 <- unlist(concordance_by_type8_is(ccRCC4_columns)) 

#Create a dataframe with adjusted p-values
immunescore_qvalues <- data.frame(BRCA1 = p.adjust(BRCA1_concordance8,method = "BH"),
                                  BRCA2 = p.adjust(BRCA2_concordance8,method = "BH"),
                                  COAD = p.adjust(COAD_concordance8,method = "BH"),
                                  DLBCL = p.adjust(DLBCL_concordance8,method = "BH"),
                                  GBM = p.adjust(GBM_concordance8,method = "BH"),
                                  HurthleCC = p.adjust(HurthleCC_concordance8,method = "BH"),
                                  HCC = p.adjust(HCC_concordance8,method = "BH"),
                                  ICC = p.adjust(ICC_concordance8,method = "BH"),
                                  OV = p.adjust(OV_concordance8,method = "BH"),
                                  PDAC = p.adjust(PDAC_concordance8,method = "BH"),
                                  PRAD = p.adjust(PRAD_concordance8,method = "BH"),
                                  ccRCC1 = p.adjust(ccRCC1_concordance8,method = "BH"),
                                  ccRCC2 = p.adjust(ccRCC2_concordance8,method = "BH"),
                                  ccRCC3 = p.adjust(ccRCC3_concordance8,method = "BH"),
                                  ccRCC4 = p.adjust(ccRCC4_concordance8,method = "BH"))
immunescore_qvalues <- immunescore_qvalues[, sort(colnames(immunescore_qvalues))]
immunescore_qvalues$metabolite = metabolites_8
immunescore_qvalues$immune_signature = rep("ImmuneScore", times = 276)
immunescore_qvalues$pair <- paste(immunescore_qvalues$metabolite, immunescore_qvalues$immune_signature, sep = "-")


#Calculate concordance between kynurenine and aDC signature (without IDO1)
# sample size of all cohorts
ns_t <- sapply(1:15,function(i){ncol(met_data_tumors[[i]])})
nsum_t <- cumsum(ns_t)
# vector with dataset id for each sample in the joint dataframe
dataset_id_t <- lapply(seq(ns_t), function(i) rep(i, ns_t[i])) %>% unlist
# weight vector for each sample in the joint dataframe
weights_t <- lapply(seq(ns_t), function(i) rep(1/ns_t[i], ns_t[i])) %>% unlist

adc <- as.numeric(unlist(adc_data_sub))
which(rownames(joint_metabolite_df_t)=="kynurenine") #169
kynurenine <- as.numeric(joint_metabolite_df_t[169,])

adc_full <- adc[!is.na(kynurenine)]
kynurenine_full <- kynurenine[!is.na(kynurenine)]
dataset_id_full <- dataset_id_t[!is.na(kynurenine)]
weights_full <- weights_t[!is.na(kynurenine)]

# compute concordance
adc_kyn_conc <- tryCatch(survival::concordance(kynurenine_full ~ adc_full + strata(dataset_id_full), keepstrata = T, weights=weights_full),
                         error = function(e) return(NA))
adc_kyn_concordance <- adc_kyn_conc$concordance #0.5865804
adc_kyn_concordance_scaled <- (2*adc_kyn_concordance)-1 #0.1731607
adc_kyn_pvalue <- tryCatch(2*pnorm(-abs(abs(adc_kyn_conc$concordance - 0.5)/sqrt(adc_kyn_conc$var))),error = function(e) return(NA)) #2.332743e-06



#### Run analysis for Immune Score analysis ----
#Panel a
#Calculate the percentage of metabolites significantly associated with Immune Score in each dataset
immunescore_percentage_sig <- sapply(1:15,function(i){length(which(subset(immunescore_qvalues,immune_signature =="ImmuneScore")[,i] < 0.05))/
    length(which(!is.na(subset(immunescore_qvalues,immune_signature =="ImmuneScore"))[,i]))})
immunescore_percentage_sig <- data.frame(cancer = colnames(immunescore_qvalues)[1:15],
                                         percent = immunescore_percentage_sig)
immunescore_percentage_sig$n_samples <- sapply(1:15,function(i){ncol((met_data)[tumors][[i]])})
immunescore_percentage_sig$cancer <- factor(immunescore_percentage_sig$cancer, 
                                            levels = immunescore_percentage_sig[order(immunescore_percentage_sig$percent),1])

#Immune Score expression by dataset
n_tumors <- sapply(1:15,function(i){nrow(tme_data[[tumors[i]]])}) #sample size for each tumor dataset
tumor_names <- sub("_Tumor", "", names(tme_data)[tumors])  
is_level <- unlist(sapply(1:15,function(i){tme_data[[tumors[i]]][,"ImmuneScore"]})) #ImmuneScore expression of every sample

is_expression <- data.frame(dataset = rep(tumor_names,times=n_tumors),
                            expression = is_level)
is_expression$dataset <- factor(is_expression$dataset, levels = levels(immunescore_percentage_sig$cancer)) #same dataset order as barplot in panel a

#Table of concordance of every metabolite in every dataset with p-values
#only metabolites in 8+ datasets
immunescore_concordance_subset <- subset(tumor_concordance,immune_signature == "ImmuneScore" & n_dataset > 7)[,c(1,8:22)]
immunescore_concordance_values <- melt(immunescore_concordance_subset) #melt so it can be used in ggplot

#Melt the corrected p-values so it can be added to the concordance dataframe
immunescore_qvalues_melt <- melt(immunescore_qvalues[,1:15]) 

colnames(immunescore_concordance_values) <- c("metabolite","dataset","concordance")
immunescore_concordance_values$dataset <- as.factor(immunescore_concordance_values$dataset)
#reorder dataset levels to be the same as in previous panel a plots
immunescore_concordance_values$dataset <- factor(immunescore_concordance_values$dataset, levels = levels(immunescore_percentage_sig$cancer)) 
immunescore_concordance_values$q_value <- immunescore_qvalues_melt$value #add adjusted p-values
#if the concordance is not significant, set to 0
immunescore_concordance_values$concordance_new <- ifelse(immunescore_concordance_values$q_value < 0.05, immunescore_concordance_values$concordance,0) 

#Compare median Immune Score expression to the percentage of significantly associated metabolites
is_median <- sapply(immunescore_percentage_sig$cancer,function(i){median(is_expression[is_expression$dataset==i,2])})
cor.test(is_median, immunescore_percentage_sig$percent, method = "spearman") #p-value = 0.6934; rho = -0.1111118


#Panel b
#Find metabolites that are siginificant in only HCC or ICC but are measured in both datasets
IS_top2_table <- na.omit(immunescore_qvalues[immunescore_qvalues$immune_signature=="ImmuneScore",c(10,12,16)]) 
#metabolites in HCC that are most significant that aren't also significant in ICC
HCC_metabolites <- c("1-linoleoyl-GPC (18:2)","1-palmitoleoyl-GPC (16:1)") 
#metabolites in ICC that are most significant that aren't also significant in HCC
ICC_metabolites <- c("thymine","nicotinamide adenine dinucleotide reduced (NADH)","UDP-glucuronate") 

#Matrix of the top metabolites in HCC and ICC to make expression heatmaps
HCC_matrix <- met_data$HCC_Tumor[c(HCC_metabolites,ICC_metabolites),]
HCC_matrix <- HCC_matrix[,order(tme_data$HCC_Tumor[,1])] #sort by increasing ImmuneScore expression
HCC_matrix <- t(scale(t(HCC_matrix)))

HCC_matrix_melt <- melt(HCC_matrix)
HCC_matrix_melt$Var2 <- factor(HCC_matrix_melt$Var2, levels = colnames(HCC_matrix))
colnames(HCC_matrix_melt) <- c("metabolite","variable","value")

ICC_matrix <- met_data$ICC_Tumor[c(HCC_metabolites,ICC_metabolites),]
ICC_matrix <- ICC_matrix[,order(tme_data$ICC_Tumor[,1])] #sort by increasing ImmuneScore expression
ICC_matrix <- t(scale(t(ICC_matrix)))

ICC_matrix_melt <- melt(ICC_matrix)
ICC_matrix_melt$Var2 <- factor(ICC_matrix_melt$Var2, levels = colnames(ICC_matrix))
colnames(ICC_matrix_melt) <- c("metabolite","variable","value")

#Immune Score expression of samples in HCC and ICC
HCC_IS <- data.frame(sample = rownames(tme_data$HCC_Tumor)[order(tme_data$HCC_Tumor[,1])],
                     IS = tme_data$HCC_Tumor[order(tme_data$HCC_Tumor[,1]),1],
                     sig = rep("ImmuneScore",54))
HCC_IS$sample <- factor(HCC_IS$sample, levels = HCC_IS$sample) #sort by increasing ImmuneScore expression

ICC_IS <- data.frame(sample = rownames(tme_data$ICC_Tumor)[order(tme_data$ICC_Tumor[,1])],
                     IS = tme_data$ICC_Tumor[order(tme_data$ICC_Tumor[,1]),1],
                     sig = rep("ImmuneScore",86))
ICC_IS$sample <- factor(ICC_IS$sample, levels = ICC_IS$sample) #sort by increasing ImmuneScore expression

#Panel c
#Determine the metabolites significantly associated with Immune Score
immunescore_concordance <- subset(tumor_concordance, immune_signature == "ImmuneScore" & n_dataset > 7) #metabolites in 8+ datasets 
immunescore_concordance$p.adj <- p.adjust(immunescore_concordance$p.value, method = "BH") #adjust p-values
#Mmap the pathways the metabolites belong to
immunescore_concordance$pathway <- metabolite_mapping$H_SUB_PATHWAY[match(immunescore_concordance$metabolite, metabolite_mapping$BIOCHEMICAL)] 
is_metabolite_rank <- immunescore_concordance[order(-immunescore_concordance$p.adj),1] #ranking metabolites by significance
immunescore_concordance$metabolite <- factor(immunescore_concordance$metabolite, levels = is_metabolite_rank) 
immunescore_concordance <- immunescore_concordance[order(immunescore_concordance$p.adj),] #order by p-value
length(which(immunescore_concordance$p.adj < 0.05)) #23


#Panel d
#Quinolinate vs z-scored ImmuneScore expression
with_quinolinate <- which(sapply(1:23,function(i){which(rownames(met_data[[i]])=="quinolinate")})!=0) #which datasets measured quinolinate; 13
with_quinolinate_tumors <- intersect(tumors,with_quinolinate) #only interested in tumor datasets

quinolinate_is_exp <- data.frame(Dataset = unlist(sapply(1:10,function(i){rep(names(met_data)[with_quinolinate_tumors[i]], dim(met_data[[with_quinolinate_tumors[i]]])[2])})),
                                 Quinolinate = unlist(sapply(1:10,function(i){scale(as.numeric(met_data[[with_quinolinate_tumors[i]]]["quinolinate",]))})),
                                 ImmuneScore = as.numeric(unlist(sapply(1:10,function(i){scale(tme_data[[with_quinolinate_tumors[i]]][,1])}))))
quinolinate_is_exp$CancerType <- word(quinolinate_is_exp$Dataset,1,sep = "\\_Tumor") #remove _Tumor from names to allow for matching
quinolinate_is_exp$CancerType <- as.factor(quinolinate_is_exp$CancerType)

#NMN vs z-scored ImmuneScore expression
with_nmn <- which(sapply(1:23,function(i){which(rownames(met_data[[i]])=="nicotinamide ribonucleotide (NMN)")})!=0) #which datasets measured quinolinate; 8
with_nmn_tumors <- intersect(tumors,with_nmn)#only interested in tumor datasets

nmn_is_exp <- data.frame(Dataset = unlist(sapply(1:8,function(i){rep(names(met_data)[with_nmn_tumors[i]], dim(met_data[[with_nmn_tumors[i]]])[2])})),
                         NMN = unlist(sapply(1:8,function(i){scale(as.numeric(met_data[[with_nmn_tumors[i]]]["nicotinamide ribonucleotide (NMN)",]))})),
                         ImmuneScore = as.numeric(unlist(sapply(1:8,function(i){scale(tme_data[[with_nmn_tumors[i]]][,1])}))))
nmn_is_exp$CancerType <- word(nmn_is_exp$Dataset,1,sep = "\\_Tumor") #remove _Tumor from names to allow for matching
nmn_is_exp$CancerType <- as.factor(nmn_is_exp$CancerType)


#Panel e
#Based on absolute concordance, which pathways are the most associated with Immune Score?
immune_score_concordance_8 <- tumor_concordance[,c(1:4)]
immune_score_concordance_8 <- subset(immune_score_concordance_8, immune_signature == "ImmuneScore" & n_dataset > 7) #only keep metabolites in 8+ datasets
#Map each metabolite to a pathway
immune_score_concordance_8$pathway <- metabolite_mapping$H_SUB_PATHWAY[match(immune_score_concordance_8$metabolite,metabolite_mapping$BIOCHEMICAL)]

#Generate list of pathways represented in the 276 metabolites
intersecting_pathways <- unique(immune_score_concordance_8$pathway) #73 actual pathways because one is NA (the first one)
immune_score_concordance_8$abs_c <- abs(immune_score_concordance_8$concordance) #absolute concordance
#One-sided wilcox test to determine which pathways have significantly higher absolute concordance to Immune Score than other pathways
#One-sided Wilcox test because we're only interested in pathways that have higher than expected concordance values
is_pathway_analysis <- data.frame(pathway = intersecting_pathways[2:74],
                                  pvalue = sapply(2:74, function(i){wilcox.test(immune_score_concordance_8[immune_score_concordance_8$pathway == intersecting_pathways[i],6],
                                                                                immune_score_concordance_8[immune_score_concordance_8$pathway != intersecting_pathways[i],6], 
                                                                                alternative = "greater")$p.value}))
is_pathway_analysis$p.adj <- p.adjust(is_pathway_analysis$pvalue, method = "BH")
is_pathway_analysis <- is_pathway_analysis[is_pathway_analysis$pathway %in% is_pathway_analysis[order(is_pathway_analysis$p.adj)[1:10],1],] #cutoff bottom portion of plot
is_pathway_order <- is_pathway_analysis[order(-is_pathway_analysis$p.adj),1] #order so pathway levels are sorted from least to most significant
is_pathway_analysis$pathway <- factor(is_pathway_analysis$pathway, levels = is_pathway_order) 


#Panel f
#Find samples which are tumor and CD45 negative/CD4/8 positive
kilgour_T_45 <- colnames(kilgour_data)[grep("T_45", colnames(kilgour_data))]
kilgour_T_imm <- setdiff(colnames(kilgour_data)[grep("T.", colnames(kilgour_data))],kilgour_T_45)

#Only interested in NAD+
kiglour_metabolite_names <- metabolite_mapping$BIOCHEMICAL[match(rownames(kilgour_data),metabolite_mapping$H_HMDB)] #match HMDB IDs to metabolite names 
which(kiglour_metabolite_names == "nicotinamide adenine dinucleotide (NAD+)") #73

#Wilcox test to compare NAD+ expression in CD45- cells vs CD4/CD8+ cells in tumor samples
wilcox.test(kilgour_data[73,kilgour_T_45],kilgour_data[73,kilgour_T_imm]) #p-value = 0.0001696
mean(kilgour_data[73,kilgour_T_45]) #mean of non-immune cells = 19.26
mean(kilgour_data[73,kilgour_T_imm]) #mean of immune cells = 18.0435

#Prepare expression matrix for plotting of the Kilgour data
kilgour_nad <- data.frame(cell_type = c(rep("CD45",length(kilgour_T_45)),rep("CD4/8",length(kilgour_T_imm))),
                          expression = c(kilgour_data[73,kilgour_T_45],kilgour_data[73,kilgour_T_imm]))
kilgour_nad$cell_type <- factor(kilgour_nad$cell_type, levels = c("CD45","CD4/8"))

#Prepare expression matrix for plotting of the cAMP ovarian data
which(rownames(met_data$OV_Tumor)=="nicotinamide adenine dinucleotide (NAD+)") #256 

ov_nad <- data.frame(NAD = as.numeric(met_data$OV_Tumor[256,]),
                     ImmuneScore = tme_data$OV_Tumor[,1])



#### Run analysis for Immune Cell Populations analysis ----
#Tumor concordance matrix with only immune cell signatures
bindea_concordance <- tumor_concordance[tumor_concordance$immune_signature %in% bindea_sigs,]
bindea_concordance <- bindea_concordance[bindea_concordance$n_dataset > 7,] #only keep metabolites measured in 8+ datasets
bindea_concordance$p.adj <- p.adjust(bindea_concordance$p.value, method = "BH") #re-correcting p-values based on this set of signatures
bindea_concordance$pair <- paste(bindea_concordance$metabolite,bindea_concordance$immune_signature,sep="-")
bindea_concordance$color <- ifelse(bindea_concordance$p.adj < 0.05, "Other","Not Significant")

#Panel a
#List of pairs to be labelled in the volcano plot
label_pairs_bindea <- c("kynurenine-aDC","quinolinate-T cells","histamine-Mast cells")

#Prepare matrix to be used for rug plot
#Plot only pairs which are associated with quinolinate
bindea_concordance$col2<- NA
bindea_concordance$col2[which(bindea_concordance$metabolite =="quinolinate")] <- 1
bindea_concordance$col3 <- c("#B3CDE3")[bindea_concordance$col2]

#Panel b
bindea_concordance_dots <- bindea_concordance
#Find max concordance for each signature
max_concordance <- sapply(bindea_sigs,function(i){max(bindea_concordance_dots[bindea_concordance_dots$immune_signature==i,4])}) 
max_concordance <- sort(max_concordance,decreasing = TRUE)
#Re-order signature factor levels so so will it plot signatures from highest max to lowest
bindea_concordance_dots$immune_signature <- factor(bindea_concordance_dots$immune_signature,levels = names(max_concordance)) 

#Panel c
#Histmaine vs z-scored Mast Cell expression
with_histamine <- which(sapply(1:23,function(i){which(rownames(met_data[[i]])=="histamine")})!=0) #which datasets measured histamine; 13
with_histamine_tumors <- intersect(tumors,with_histamine) #only interested in tumor datasets; 13

histamine_mc_exp <- data.frame(Dataset = unlist(sapply(1:13,function(i){rep(names(met_data)[with_histamine_tumors[i]], dim(met_data[[with_histamine_tumors[i]]])[2])})),
                               Histamine = unlist(sapply(1:13,function(i){scale(as.numeric(met_data[[with_histamine_tumors[i]]]["histamine",]))})),
                               Mast.Cells = as.numeric(unlist(sapply(1:13,function(i){scale(tme_data[[with_histamine_tumors[i]]][,23])}))))
histamine_mc_exp$CancerType <- word(histamine_mc_exp$Dataset,1,sep = "\\_Tumor") #remove _Tumor from names to allow for matching
histamine_mc_exp$CancerType <- as.factor(histamine_mc_exp$CancerType)


#Panel d
#Concordance of metabolites with mast cells
mastcell_concordance <- subset(tumor_concordance, immune_signature == "Mast cells" & n_dataset > 4) #metabolites in 4+ datasets
mastcell_concordance$p.adj <- p.adjust(mastcell_concordance$p.value, method = "BH") #re-adjust to just mast cell concordances
#Re-order metabolite factor levels from largest to smallest
mastcell_concordance$metabolite <- factor(mastcell_concordance$metabolite, levels = mastcell_concordance[order(-mastcell_concordance$p.adj),1]) 
mastcell_concordance <- mastcell_concordance[order(mastcell_concordance$p.adj),]
length(which(mastcell_concordance$p.adj < 0.05)) #how many significant metabolites; 21


#Panel e
#Histmaine vs z-scored HDC expression
all(names(met_data)[with_histamine_tumors] == names(rna_data)) #make sure datasets are in the same order
histamine_hdc_exp <- data.frame(Dataset = unlist(sapply(1:13,function(i){rep(names(rna_data)[[i]], ncol(met_data[[with_histamine_tumors[i]]]))})),
                                Histamine = unlist(sapply(1:13,function(i){scale(as.numeric(met_data[[with_histamine_tumors[i]]]["histamine",]))})),
                                HDC = as.numeric(unlist(sapply(1:13,function(i){scale(as.numeric(rna_data[[i]]["HDC",]))}))))
histamine_hdc_exp$CancerType <- word(histamine_hdc_exp$Dataset,1,sep = "\\_Tumor") #remove _Tumor from names to allow for matching
histamine_hdc_exp$CancerType <- as.factor(histamine_hdc_exp$CancerType)

#calculate the range of histamine expression
histamine_range <- as.data.frame(sapply(1:13,function(i){range(as.numeric(met_data[[with_histamine_tumors[i]]]["histamine",]))})) #range of histamine in all datasets
mean(as.numeric(2^(abs(histamine_range[1,]-histamine_range[2,])))) #737.023


#Panel f
#Kynurenine vs z-scored aDC (without IDO1) expression
with_kynurenine <- which(sapply(1:15,function(i){which(rownames(met_data_tumors[[i]])=="kynurenine")})!=0) #12 datasets

adc_kyn_exp <- data.frame(Dataset = unlist(sapply(1:12,function(i){rep(names(met_data_tumors)[with_kynurenine[i]], dim(met_data_tumors[[with_kynurenine[i]]])[2])})),
                          Kynurenine = unlist(sapply(1:12,function(i){scale(as.numeric(met_data_tumors[[with_kynurenine[i]]]["kynurenine",]))})),
                          aDC = as.numeric(unlist(sapply(1:12,function(i){scale(adc_data_sub[[with_kynurenine[i]]])}))))
adc_kyn_exp$CancerType <- word(adc_kyn_exp$Dataset,1,sep = "\\_Tumor") #remove MetabImmune_tumor.csv from names to allow for matching
adc_kyn_exp$CancerType <- as.factor(adc_kyn_exp$CancerType)


#Supplemental figure 2
#boxplot of mast cell-histamine concordance for each dataset
histamine_concordance <- tumor_concordance[tumor_concordance$immune_signature=="Mast cells" & tumor_concordance$metabolite=="histamine",8:22]
histamine_concordance <- melt(histamine_concordance)
histamine_concordance$pair <- c("Histamine-Mast Cells")


#### Save data to file ----

save.image("Workspace_8_ImmuneAnalysis.Rdata")




