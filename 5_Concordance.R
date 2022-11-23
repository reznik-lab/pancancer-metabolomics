#### Initialize ----

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(readxl)
library(tidyverse)
library(magrittr)
library(purrr)
library(parallel)
library(survival)

#### Load Data ----

# load preprocessed and filtered data
load("results/Workspace_3_FilterMetabo.Rdata")
load("results/Workspace_4_FilterRNA.Rdata")

# reorder lists in alphabetical order
met_all <- met_all[sort(names(met_all))]
rna_all <- rna_all[sort(names(rna_all))]

sapply(met_all, ncol)
sapply(rna_all, ncol)

#### Scale Data ----

met_scaled <- lapply(met_all, function(y){
  mat <- y %>% t %>% scale %>% t %>% as.data.frame
  rownames(mat) <- rownames(y)
  mat
})

rna_scaled <- lapply(rna_all, function(y){
  mat <- y %>% t %>% scale %>% t %>% as.data.frame
  rownames(mat) <- rownames(y)
  mat
})

#### Create auxiliary variables ----

# create a mapping dataframe
df <- met_all %>% sapply(ncol) %>% as.data.frame %>% 
  dplyr::rename(n=".") %>% 
  tibble::rownames_to_column("dataset") %>%
  dplyr::arrange(dataset) %>%
  dplyr::mutate(tissue=ifelse(grepl(pattern = "_Tumor", x=dataset, fixed = T),"Tumor","Normal"))

# sample size vector for tumor datasets
ns <- sapply(rna_scaled[df$dataset[df$tissue=="Tumor"]], ncol) %>% as.vector()
# dataset ids for tumor datasets
z <- lapply(seq(ns), function(i) rep(i, ns[i])) %>% unlist
# weight vector for tumor datasets
w <- lapply(seq(ns), function(i) rep(1/ns[i], ns[i])) %>% unlist

#### Define Gene Subset ----

# genes with measured in all datasets
gg <- lapply(rna_scaled, rownames) %>% 
  purrr::reduce(intersect)

#### Join Data ----

met_joint <- met_scaled[df$dataset[df$tissue=="Tumor"]] %>%
  lapply(function(y){
    y %>% tibble::rownames_to_column("metabolite")
  }) %>%
  reduce(full_join, by="metabolite") %>%
  tibble::column_to_rownames("metabolite")

rna_joint <- rna_scaled[df$dataset[df$tissue=="Tumor"]] %>%
  lapply(function(y){
    y %>% tibble::rownames_to_column("gene")
  }) %>%
  reduce(full_join, by="gene") %>%
  tibble::column_to_rownames("gene")
  
#### Clean up Workspace ----

remove(rna_all, met_all, rna_scaled, met_scaled)

#### Compute Concordance ----

warning("The concordance calculation is computationally intense and will take >12 hours on a typical laptop.")

tic()
# compute pairwise concordance
conc <- mclapply(gg[1:10], function(g){
  mclapply(rownames(met_joint), function(m){
    
    # gene values
    x <- rna_joint[g,] %>% unlist
    # metabolite values
    y <- met_joint[m,] %>% unlist
    
    # remove NA
    x_full <- x[!is.na(y)]
    y_full <- y[!is.na(y)]
    z_full <- z[!is.na(y)]
    w_full <- w[!is.na(y)]
    
    # compute concordance
    conc <- survival::concordance(y_full ~ x_full + strata(z_full), keepstrata = T, weights=w_full)
    
    # extract concordance values for each dataset
    if (class(conc$count)!="numeric"){
      individual_conc <- conc$count %>% {.[,1]/rowSums(.[,1:3])} %>% data.frame %>% t %>% data.frame
      colnames(individual_conc) <- df$dataset[df$tissue=="Tumor"][z_full %>% unique]
    } else {
      individual_conc <- conc$count %>% as.matrix() %>% t %>% data.frame %>%
        {.[,1]/rowSums(.[,1:3])} %>% data.frame %>% t %>% data.frame
      colnames(individual_conc) <- df$dataset[df$tissue=="Tumor"][z_full %>% unique]
    }
    
    # arrange dataframe
    X <- setNames(data.frame(matrix(ncol = df$dataset[df$tissue=="Tumor"] %>% length, nrow = 0)), 
                  df$dataset[df$tissue=="Tumor"]) %>%
      dplyr::full_join(individual_conc, by=colnames(individual_conc))
    
    # create result entry
    data.frame(metabolite = m,
               gene = g,
               n_met = y_full %>% length,
               n_dataset = z_full %>% unique %>% length,
               concordance = conc$concordance,
               variance = conc$var,
               zscore = (conc$concordance - 0.5)/sqrt(conc$var)) %>%
      cbind(X)
    
  }, mc.cores = 1) %>% {do.call(rbind,.)}
}, mc.cores = 4) %>% {do.call(rbind,.)}
toc()

#### Compute Pvalues ----

conc %<>% 
  dplyr::filter(n_dataset>=2) %>%
  dplyr::mutate(p.value=2*pnorm(-abs(zscore))) %>% 
  dplyr::mutate(padj=p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(padj,-n_dataset,-concordance) 

#### Add distance ----

# load precalculated distance 
load("data_for_scripts/Concordance/distance.Rdata")

# merge pathway distance with concordance results
conc %<>%
  dplyr::left_join(d_GEM_long, by=c("gene","metabolite"="name")) %>%
  dplyr::rename(distance=dist_GEM_min) %>%
  # reorder columns for better readability
  dplyr::select(metabolite, gene, distance,n_met, n_dataset, 
                concordance, variance, zscore, p.value, padj,
                everything())

#### Save Results ----

save(conc, file="results/Workspace_5_ConcordanceMetaAnalysis.Rdata")
