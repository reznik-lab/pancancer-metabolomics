#### Initialize ----

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(openxlsx)
library(reshape2)
library(maplet)

#### Initialize parameters ----

# define parameter for normalization
met_maxmiss <- 0.2
# define parameter for filtering
metMax <- 1

#### Define Functions ----

harmonize_metnames <- function(D){
  # get original metabolites names
  g <- data.frame(BIOCHEMICAL=rowData(D)$BIOCHEMICAL)
  g$BIOCHEMICAL <- as.character(g$BIOCHEMICAL)
  # get substitutions from external file
  p <- g %>%
    dplyr::left_join(new_names, by=c("BIOCHEMICAL"="OLD_NAME"))
  
  # substitute names
  g$BIOCHEMICAL[!is.na(p$NEW_NAME)] <- p$NEW_NAME[!is.na(p$NEW_NAME)]
  g$BIOCHEMICAL
}

harmonize_metanno <- function(D){
  # get original metabolite annotations
  g <- as.data.frame(rowData(D)) 
  g <- g %>% 
    dplyr::left_join(metanno, by=c("H_name"="BIOCHEMICAL")) %>%
    # remove asterisks from harmonized names
    dplyr::mutate(H_name=gsub(pattern="\\*", replacement = "", x = H_name))
  g
}

# Preprocessing Pipeline
preprocessing_pipeline <- function(D){
  D %>%
    # missingness plot
    mt_plots_missingness() %>%
    # filter
    mt_pre_filter_missingness(feat_max = metMax, group_col = "GROUP") %>% 
    # sample boxplot
    mt_plots_sample_boxplot(color=GROUP, title='before normalization', plot_logged=T) %>% 
    # quotient normalization based on controls and on metabolites with less than 20% missingness
    mt_pre_norm_quot(ref_samples=GROUP=="NORMAL", feat_max = met_maxmiss) %>% 
    # plot dilution coefficients
    mt_plots_dilution_factor(in_col="GROUP") %>%
    # sample boxplot
    mt_plots_sample_boxplot(color=GROUP, title='after quotient normalization', plot_logged=T) %>% 
    # log
    mt_pre_trans_log()
}

# Preprocessing Pipeline
preprocessing_pipeline_NoNormal <- function(D){
  D %>%
    # missingness plot
    mt_plots_missingness() %>%
    # filter
    mt_pre_filter_missingness(feat_max = metMax, group_col = "GROUP") %>% 
    # sample boxplot
    mt_plots_sample_boxplot(color=GROUP, title='before normalization', plot_logged=T) %>% 
    # quotient normalization based on controls and on metabolites with less than 20% missingness
    mt_pre_norm_quot(feat_max = met_maxmiss) %>% 
    # plot dilution coefficients
    mt_plots_dilution_factor(in_col="GROUP") %>%
    # sample boxplot
    mt_plots_sample_boxplot(color=GROUP, title='after quotient normalization', plot_logged=T) %>% 
    # log
    mt_pre_trans_log()
}

# write SE to xls
writeSE2xls <- function(D, name) {
  
  wb = createWorkbook()
  sheet = addWorksheet(wb, "data")
  writeData(wb, sheet=sheet, assay(D), colNames = T, rowNames = T)
  sheet = addWorksheet(wb, "metanno")
  writeData(wb, sheet=sheet, as.data.frame(rowData(D)), colNames = T, rowNames = F)
  sheet = addWorksheet(wb, "sampleanno")
  writeData(wb, sheet=sheet, as.data.frame(colData(D)), colNames = T, rowNames = F)
  saveWorkbook(wb, paste0(name,".xlsx", sep=""), overwrite = TRUE)
}

#### Load Data ----

# load files with name harmonizations and annotations
new_names <- read.xlsx("data/MetinfoFromData.xlsx", sheet = "Sheet2")
metanno <- read.xlsx("data/MetinfoFromData.xlsx", sheet = "Sheet3")

# load master mapping file
mapping_file <- read.csv2(file = "data/MasterMapping_MetImmune_03_16_2022_release.csv", sep=",") 

# collect all datasets into one list
D <- list()

# set data path
path_data ="data/metabolomics_original/"

# Load Breast Cancer (Terunuma) Data 
D$BRCA1 <-
  mt_load_metabolon_v1(file=paste(path_data,"BRCA1.xlsx",sep=""), sheet="OrigData") %>%
  mt_anno_mutate(anno_type = "samples",col_name = "SAMPLE_NAME", term = SAMPLE_ID) %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term = case_when(TISSUE.TYPE=="TUMOR" ~ "TUMOR",
                                                                            TISSUE.TYPE=="NORMAL" ~ "NORMAL") %>% as.factor)
colnames(D$BRCA1) <- as.character(colData(D$BRCA1)$SAMPLE_NAME)

# Load Breast Cancer (Tang) Data 
D$BRCA2 <- 
  mt_load_metabolon_v1(file=paste(path_data,"BRCA2.xlsx",sep=""), sheet="OrigScale") %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term = case_when(TCGA.Designation!="Normal Breast" ~ "TUMOR",
                                                                             TCGA.Designation=="Normal Breast" ~ "NORMAL") %>% as.factor)

# Load COAD Data
D$COAD <- 
  mt_load_xls(file= paste(path_data,"COAD.xlsx",sep=""),sheet = "data", samples_in_rows=T, id_col="SAMPLE_NAME") %>%
  mt_anno_xls(file= paste(path_data,"COAD.xlsx",sep=""), sheet="sampleinfo", anno_type="samples", anno_id_col="SAMPLE_NAME") %>% 
  mt_anno_xls(file= paste(path_data,"COAD.xlsx",sep=""), sheet="metinfo", anno_type="features", anno_id_col="BIOCHEMICAL", data_id_col="name") %>%
  mt_pre_zero_to_na() %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term = case_when(Tissue=="Tumour" ~ "TUMOR",
                                                                             Tissue=="Normal" ~ "NORMAL") %>% as.factor)

# Load DLBCL Data
D$DLBCL <- mt_load_xls(file= paste(path_data,"DLBCL.xlsx",sep=""),sheet = "data", samples_in_rows=T, id_col="SAMPLE_NAME") %>%
  mt_anno_xls(file= paste(path_data,"DLBCL.xlsx",sep=""), sheet="sampleanno", anno_type="samples", anno_id_col="SAMPLE_NAME") %>% 
  mt_anno_xls(file= paste(path_data,"DLBCL.xlsx",sep=""), sheet="metanno", anno_type="features", anno_id_col="BIOCHEMICAL", data_id_col="name") %>%
  mt_anno_mutate(anno_type = "features", col_name = "HMDb", term=HMDb_ID) %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term = case_when(GROUP_DESC!="control" ~ "TUMOR",
                                                                             GROUP_DESC=="control" ~ "NORMAL") %>% as.factor) 

# Load Glioblastoma Data
# load metabolomics data
GBMmet <- 
  mt_load_xls(file= paste(path_data,"GBM.xlsx",sep=""),sheet = "data", samples_in_rows=F, id_col="Metabolite") %>%
  mt_anno_xls(file= paste(path_data,"GBM.xlsx",sep=""), sheet="sampleanno", anno_type="samples", anno_id_col="case_id",data_id_col = "sample") %>%
  mt_anno_mutate(anno_type = "features", col_name = "BIOCHEMICAL", term = name) %>%
  mt_anno_mutate(anno_type = "samples", col_name = "SAMPLE_NAME", term = "sample") %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term = case_when(!is.na(tumor_laterality) ~ "TUMOR",
                                                                             is.na(tumor_laterality) ~ "NORMAL") %>% as.factor)
# load lipidomics1 data
GBMlipid1 <- 
  mt_load_xls(file= paste(path_data,"GBM.xlsx",sep=""),sheet = "lipidome_positive_normalized", samples_in_rows=F, id_col="Lipid") %>%
  mt_anno_xls(file= paste(path_data,"GBM.xlsx",sep=""), sheet="sampleanno", anno_type="samples", anno_id_col="case_id",data_id_col = "sample") %>%
  mt_anno_mutate(anno_type = "features", col_name = "BIOCHEMICAL", term = name) %>%
  mt_anno_mutate(anno_type = "samples", col_name = "SAMPLE_NAME", term = "sample") %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term = case_when(!is.na(tumor_laterality) ~ "TUMOR",
                                                                             is.na(tumor_laterality) ~ "NORMAL") %>% as.factor)
# load lipidomics2 data
GBMlipid2 <- 
  mt_load_xls(file= paste(path_data,"GBM.xlsx",sep=""),sheet = "lipidome_negative_normalized", samples_in_rows=F, id_col="Lipid") %>%
  mt_anno_xls(file= paste(path_data,"GBM.xlsx",sep=""), sheet="sampleanno", anno_type="samples", anno_id_col="case_id",data_id_col = "sample") %>%
  mt_anno_mutate(anno_type = "features", col_name = "BIOCHEMICAL", term = name) %>%
  mt_anno_mutate(anno_type = "samples", col_name = "SAMPLE_NAME", term = "sample") %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term = case_when(!is.na(tumor_laterality) ~ "TUMOR",
                                                                             is.na(tumor_laterality) ~ "NORMAL") %>% as.factor)

# put all together
rd <- GBMmet %>% rowData %>% as.data.frame %>% 
  rbind(GBMlipid1 %>% rowData %>% as.data.frame) %>%
  rbind(GBMlipid2 %>% rowData %>% as.data.frame)

cd <- GBMmet %>% colData %>% as.data.frame %>% 
  dplyr::left_join(GBMlipid1 %>% colData %>% as.data.frame) %>%
  dplyr::left_join(GBMlipid2 %>% colData %>% as.data.frame) %>%
  dplyr::mutate(SAMPLE_NAME=sample)

dt <- GBMmet %>% assay %>% t %>% as.data.frame %>%
  tibble::rownames_to_column("SAMPLE_NAME") %>%
  dplyr::left_join(GBMlipid1 %>% assay %>% t %>% as.data.frame %>% tibble::rownames_to_column("SAMPLE_NAME"),
                   by="SAMPLE_NAME") %>%
  dplyr::left_join(GBMlipid2 %>% assay %>% t %>% as.data.frame %>% tibble::rownames_to_column("SAMPLE_NAME"),
                   by="SAMPLE_NAME") %>%
  tibble::column_to_rownames("SAMPLE_NAME") %>%
  t %>% as.data.frame

# store SE
D$GBM <- SummarizedExperiment(assays = dt, rowData = rd, colData = cd)

# Load Thyroid Cancer Data 
D$HurthleCC <-
  mt_load_metabolon_v1(file=paste(path_data,"HurthleCC.xlsx",sep=""), sheet="OrigScale (Tissue)") %>%
  mt_anno_mutate(anno_type = "features", col_name = "PUBCHEM", term = as.numeric(PUBCHEM)) %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term = case_when(TISSUE=="TUMOR" ~ "TUMOR",
                                                                             TISSUE=="NORMAL" ~ "NORMAL") %>% as.factor)

# Load Liver Cancer Data (HCC)
D$HCC <- 
  mt_load_xls(file=paste(path_data,"LiverCancer.xlsx",sep=""), sheet="data_HCC", samples_in_rows=T, id_col="SAMPLE_NAME") %>% 
  mt_anno_xls(file=paste(path_data,"LiverCancer.xlsx",sep=""), sheet="sampleanno", anno_type="samples", anno_id_col="SAMPLE_NAME",data_id_col = "SAMPLE_NAME") %>%
  mt_anno_mutate(anno_type = "features", col_name = "BIOCHEMICAL", term = name) %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term=as.factor("TUMOR"))

# Load Liver Cancer Data (ICC)
D$ICC <- 
  mt_load_xls(file=paste(path_data,"LiverCancer.xlsx",sep=""), sheet="data_ICC", samples_in_rows=T, id_col="SAMPLE_NAME") %>% 
  mt_anno_xls(file=paste(path_data,"LiverCancer.xlsx",sep=""), sheet="sampleanno", anno_type="samples", anno_id_col="SAMPLE_NAME",data_id_col = "SAMPLE_NAME") %>%
  mt_anno_mutate(anno_type = "features", col_name = "BIOCHEMICAL", term = name) %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term=as.factor("TUMOR"))

# Load Ovarian Cancer Data
D$OV <- 
  mt_load_xls(file=paste(path_data,"OV.xlsx",sep=""), sheet="data", samples_in_rows=TRUE, id_col="SAMPLE_ID") %>%
  mt_anno_xls(file=paste(path_data,"OV.xlsx",sep=""), sheet="metinfo", anno_type="features", anno_id_col="BIOCHEMICAL", data_id_col="name") %>%
  mt_anno_mutate(anno_type = "samples", col_name = "SAMPLE_NAME", term = SAMPLE_ID) %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term=as.factor("TUMOR"))

# Load Pancreatic Cancer Data 
D$PDAC <- 
  mt_load_metabolon_v1(file=paste(path_data,"PDAC.xlsx",sep=""), sheet="OrigScale") %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term = case_when(Group!="Nontumor (N)" ~ "TUMOR",
                                                                             Group=="Nontumor (N)" ~ "NORMAL") %>% as.factor)

# Load Prostate Cancer Data 
D$PRAD <- 
  mt_load_xls(file=paste(path_data,"PRAD.xlsx",sep=""), sheet="data", samples_in_rows=T, id_col="SAMPLE_NAME") %>% 
  mt_anno_xls(file=paste(path_data,"PRAD.xlsx",sep=""), sheet="sampleinfo", anno_type="samples", anno_id_col="SAMPLE_NAME") %>% 
  mt_anno_xls(file=paste(path_data,"PRAD.xlsx",sep=""), sheet="metinfo", anno_type="features", anno_id_col="BIOCHEMICAL", data_id_col="name") %>%
  mt_anno_mutate(anno_type = "features", col_name = "HMDb", term=HMDb_ID) %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term = case_when(T.or.N.=="T" ~ "TUMOR",
                                                                             T.or.N.=="N" ~ "NORMAL") %>% as.factor) %>%
  mt_modify_filter_samples(filter = !is.na(GROUP))

# Load Kidney Cancer 1 Data
D$ccRCC1 <- 
  mt_load_metabolon_v1(file=paste(path_data,"ccRCC1.xlsx",sep=""), sheet="OrigScale") %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term = case_when(TISSUE.TYPE=="T" ~ "TUMOR",
                                                                             TISSUE.TYPE=="N" ~ "NORMAL") %>% as.factor)

# Load Kidney Cancer 2/3 Data
D$ccRCC2 <- 
  mt_load_metabolon_v1(file=paste(path_data,"ccRCC2.xlsx",sep=""), sheet="OrigScale") %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term = case_when(TISSUE_STATUS=="TISSUE_TUMOR" ~ "TUMOR",
                                                                             TISSUE_STATUS=="TISSUE_NORMAL" ~ "NORMAL") %>% as.factor)

# Load Kidney Cancer 4
D$ccRCC4 <- 
  mt_load_metabolon_v1(file=paste(path_data,"ccRCC4.xlsx",sep=""), sheet="OrigScale") %>%
  mt_anno_mutate(anno_type = "samples", col_name = "GROUP", term = case_when(GROUP_NUMBER %in% c(1,3,5) ~ "TUMOR",
                                                                             GROUP_NUMBER %in% c(2,4,6) ~ "NORMAL") %>% as.factor)

#### Preprocessing ----

# rearrange datasets alphabetically
D <- D[sort(names(D))]

# preprocess
lapply(names(D) %>% {names(.)=.;.}, function(x){
  print(x)

  # Glioblastoma, and Liver cancer data are already normalized
  if(!(x %in% c("GBM","HCC","ICC","OV"))) {
    # for PRAD we need batch correction
    if(x=="PRAD"){
      D[[x]] <<- D[[x]] %>%
        # missingness plot
        mt_plots_missingness() %>%
        # filter
        mt_pre_filter_missingness(feat_max = metMax, group_col = "GROUP") %>% 
        # sample boxplot
        mt_plots_sample_boxplot(color=GROUP, title='before normalization', plot_logged=T) %>%
        # create batch variable
        mt_anno_mutate(anno_type="samples",col_name="Batch_correction",term=ifelse(Metabolon.Date=="2019-12-12",sprintf("%s_%s",Metabolon.Date,RUN_DAY),Metabolon.Date)) %>%
        # batch correction
        mt_pre_batch_median(batch_col = "Batch_correction") %>%
        # sample boxplot
        mt_plots_sample_boxplot(color=GROUP, title='before normalization and after batch correction', plot_logged=T) %>% 
        # quotient normalization based on controls and on metabolites with less than 20% missingness
        mt_pre_norm_quot(ref_samples=GROUP=="NORMAL", feat_max = met_maxmiss) %>% 
        # plot dilution coefficients
        mt_plots_dilution_factor(in_col="GROUP") %>%
        # sample boxplot
        mt_plots_sample_boxplot(color=GROUP, title='after quotient normalization', plot_logged=T) %>% 
        # log
        mt_pre_trans_log()
    }
    if("NORMAL" %in% D[[x]]$GROUP & x!="PRAD" & x!="DLBCL"){
      D[[x]] <<- D[[x]] %>%
        preprocessing_pipeline()
    }
    if(!("NORMAL" %in% D[[x]]$GROUP) | x=="DLBCL"){
      D[[x]] <<- D[[x]] %>%
        preprocessing_pipeline_NoNormal()
    }
  }
  # OV is already normalized but not logged
  if(x=="OV"){
    D[[x]] <<- D[[x]] %>%
      # log
      mt_pre_trans_log()
  }
  
  # filter rowdata columns
  rowData(D[[x]]) <<- rowData(D[[x]])[,colnames(rowData(D[[x]])) %in% c("BIOCHEMICAL", "name")]
  # filter coldata columns
  colData(D[[x]]) <<- colData(D[[x]])[,colnames(colData(D[[x]])) %in% c("SAMPLE_NAME","GROUP")]
}) %>% invisible()

#### Harmonize Metabolite Names ----

# collect all initial unique metabolite names
metlist_initial <- lapply(D,function(x){rowData(x)$BIOCHEMICAL}) %>% {do.call(c,.)} %>% unique

# add harmonized metabolite names
lapply(names(D), function(x){rowData(D[[x]])$H_name <<- harmonize_metnames(D[[x]])}) %>% invisible

# collect all final unique metabolite names
metlist_final <- lapply(names(D),function(x){rowData(D[[x]])$H_name}) %>% {do.call(c,.)} %>% unique
# check that all metabolites have harmonized annotations
metlist_final[!(metlist_final %in% metanno$BIOCHEMICAL)]

# add harmonized annotations
lapply(names(D), function(x){rowData(D[[x]]) <<- harmonize_metanno(D[[x]])}) %>% invisible

#### Check Duplicates in Harmonized Names ----

# manually remove duplicated harmonized entries (same metabolite)
lapply(names(D), function(x){
  l <- length(which(duplicated(rowData(D[[x]])$H_name)))
  if(l>=1){
    R.utils::printf("%s duplicated: %d\n",x, l)
    R.utils::printf("%s\n",rowData(D[[x]])$H_name[duplicated(rowData(D[[x]])$H_name)])
  }
}) %>% invisible
# BRCA1 X - 13516, histidylleucine (previous unknowns)
# CPTAC_GBM 2-aminoadipate
# Glioma adenosine 2'-monophosphate (2'-AMP) (rename AMP)
# PRAD heme (rename heme2)
# RC20 adenosine 2'-monophosphate (2'-AMP) (rename AMP)

rowData(D$BRCA1)$H_name[duplicated(rowData(D$BRCA1)$H_name)] <- sprintf("%s _dupl", rowData(D$BRCA1)$H_name[duplicated(rowData(D$BRCA1)$H_name)])
rowData(D$PRAD)$H_name[rowData(D$PRAD)$name == "heme"] <- "heme2"
rowData(D$ccRCC4)$H_name[rowData(D$ccRCC4)$name == "AMP"] <- "AMP"
rowData(D$GBM)$H_name[rowData(D$GBM)$name=="DL-2-Aminoadipic acid (spectral match)"] <- sprintf("%s_dupl", rowData(D$GBM)$H_name[rowData(D$GBM)$name=="DL-2-Aminoadipic acid (spectral match)"])

# change rownames to harmonized names
lapply(names(D), function(x){rownames(D[[x]]) <<- rowData(D[[x]])$H_name}) %>% invisible()

# check if renaming worked
sapply(names(D), function(x){
  all.equal(rownames(D[[x]]), rowData(D[[x]])$H_name)
})

#### Write normalized data to file ----

# check if directory to save results exists, otherwise create
if(!dir.exists("results")) {
  dir.create("results")
}
# check if directory to save results exists, otherwise create
if(!dir.exists("results/preprocessed_data")) {
  dir.create("results/preprocessed_data")
}

lapply(names(D), function(x){
  writeSE2xls(D=D[[x]], name=sprintf("results/preprocessed_data/PreprocessedData_%s", x))
}) %>% invisible

#### Save Workspace to file ----

# check if directory to save results exists, otherwise create
if(!dir.exists("results")) {
  dir.create("results")
}
save(D, metanno, mapping_file,
     file = "results/Workspace_1_PreprocessMetabo.Rdata")
