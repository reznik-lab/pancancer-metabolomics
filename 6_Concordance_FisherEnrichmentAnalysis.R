#### Initialize ----

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# clear workspace
rm(list = ls())

library(tidyverse)
library(magrittr)
library(pheatmap)

#### Load Data ----

# load pathway annotations
load("data_for_scripts/Concordance/Annotations_Pathways.Rdata")
# load metabolite annotations
load("data_for_scripts/Concordance/Annotations_Metabolon_Metabolites.Rdata")

# load concordance meta-analysis results
if(file.exists("results/Workspace_5_ConcordanceMetaAnalysis.Rdata")){
  load("results/Workspace_5_ConcordanceMetaAnalysis.Rdata")
} else{
  load("data_for_scripts/Concordance/ConcordanceMetaAnalysis.Rdata")
}

#### Set global parameters ----

# significance threshold on the FDR adjusted p-values
pcut <- 0.01
# cutoff on the number of datasets a metabolite must be present in
n_data_cut <- 7

#### Prepare Data ----

# extract genes with at least one significant association
gene_sign <- conc %>%
  dplyr::filter(n_dataset>n_data_cut) %>%
  dplyr::filter(padj<pcut) %>%
  dplyr::select(gene) %>%
  distinct() %>%
  dplyr::pull(gene)

# extract all unique genes
gene_all <- conc %>%
  dplyr::filter(n_dataset>n_data_cut) %>%
  dplyr::select(gene) %>%
  distinct() %>%
  dplyr::pull(gene)

# extract significant metabolites
met_sign <- conc %>%
  dplyr::filter(n_dataset>n_data_cut) %>%
  dplyr::filter(padj<pcut) %>%
  dplyr::select(metabolite) %>%
  distinct() %>%
  dplyr::pull(metabolite)

# filter pathways with at least one measured gene in them
pw_filtered <- lapply(paths_list$kegg, function(x){
  x[x %in% gene_all]
})

#### Run Pathway Enrichment Analysis ----

# Fisher's test-based pathway enrichement analysis for each metabolite 
res_met <- lapply(met_sign %>% {names(.)=.;.}, function(y){
  print(y)
  gene_set <- conc %>%
    dplyr::filter(n_dataset>7) %>%
    dplyr::filter(padj<0.01) %>%
    dplyr::filter(metabolite==y) %>%
    dplyr::select(gene) %>%
    distinct() %>%
    dplyr::pull(gene)

  lapply(pw_filtered %>% names %>% {names(.)=.;.}, function(x){
    ft <- data.frame(gene=gene_all) %>%
      dplyr::mutate(in_pw=ifelse(gene %in% pw_filtered[[x]],1,0) %>% factor(levels=c(1,0))) %>%
      dplyr::mutate(in_dt=ifelse(gene %in% gene_set,1,0) %>% factor(levels=c(1,0))) %>%
      dplyr::select(-gene) %>%
      table %>%
      fisher.test(alternative = "greater")

    data.frame(metabolite=y,
               pathway=x,
               odds.ratio=ft$estimate %>% as.numeric(),
               p.value=ft$p.value)

  }) %>% {do.call(rbind,.)} %>% as.data.frame
}) %>% {do.call(rbind,.)} %>% as.data.frame %>%
  # add KEGG annotations
  dplyr::left_join(kegg_hierarchy, by=c("pathway"="match_name")) %>%
  # remove pathways classified as "Human Disease"
  dplyr::filter(Class != "Human Diseases") %>%
  # multiple testing correction
  dplyr::mutate(padj=p.adjust(p.value, method = "BH")) %>%
  dplyr::mutate(score=-log10(padj))

# save results to file
save(res_met, file="results/Workspace_6_ConcordancePathwayEnrichmentAnalysis.Rdata")
