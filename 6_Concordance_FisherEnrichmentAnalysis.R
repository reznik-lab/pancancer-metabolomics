#### Initialize ----

zap()

library(tidyverse)
library(magrittr)
library(pheatmap)

setwd("/Users/elb4003/Work/PostDocWCM/Projects/MetaboPanCan/MetaboRNA/")

# check if directory to save results exists, otherwise create
if(!dir.exists(paste0(Sys.Date(),sep=""))) {
  dir.create(paste0(Sys.Date(),sep=""))
}

# #### Load Data ----
# 
# load("Files2Share/Annotations_Pathways.Rdata")
# load("Files2Share/Results_Concordance_MetaAnalysis.Rdata")
# load("Files2Share/Annotations_Metabolon_Metabolites.Rdata")
# 
# #### Prepare Data ----
# 
# gene_sign <- conc %>% 
#   dplyr::filter(n_dataset>7) %>%
#   dplyr::filter(padj<0.01) %>%
#   dplyr::select(gene) %>%
#   distinct() %>%
#   dplyr::pull(gene)
# 
# gene_all <-  conc %>% 
#   dplyr::filter(n_dataset>7) %>%
#   dplyr::select(gene) %>%
#   distinct() %>%
#   dplyr::pull(gene)
# 
# geneset$kegg %>% length
# 
# met_sign <- conc %>% 
#   dplyr::filter(n_dataset>7) %>%
#   dplyr::filter(padj<0.01) %>%
#   dplyr::select(metabolite) %>%
#   distinct() %>%
#   dplyr::pull(metabolite)
# 
# pw_filtered <- lapply(paths_list$kegg, function(x){
#   x[x %in% gene_all]
# })
# 
# sapply(pw_filtered, length) %>% as.data.frame %>% 
#   dplyr::rename(n_genes_restricted=".") %>%
#   tibble::rownames_to_column("pathways") %>%
#   dplyr::left_join(sapply(paths_list$kegg, length) %>% as.data.frame %>% 
#                      dplyr::rename(n_genes_full=".") %>%
#                      tibble::rownames_to_column("pathways"),
#                    by="pathways") %>%
#   dplyr::mutate(diff=n_genes_full-n_genes_restricted) %>% View
# 
# #### Run Pathway Enrichment Analysis ----
# 
# res <- lapply(pw_filtered %>% names %>% {names(.)=.;.}, function(x){
#   ft <- data.frame(gene=gene_all) %>%
#     dplyr::mutate(in_pw=ifelse(gene %in% pw_filtered[[x]],1,0) %>% factor(levels=c(1,0))) %>%
#     dplyr::mutate(in_dt=ifelse(gene %in% gene_sign,1,0) %>% factor(levels=c(1,0))) %>% 
#     dplyr::select(-gene) %>% 
#     table %>% 
#     fisher.test(alternative = "greater")
#   
#   data.frame(pathway=x,
#              odds.ratio=ft$estimate %>% as.numeric(),
#              p.value=ft$p.value)
# }) %>% {do.call(rbind,.)} %>% as.data.frame %>%
#   dplyr::left_join(kegg_hierarchy, by=c("pathway"="match_name")) %>%
#   dplyr::filter(Class != "Human Diseases") %>%
#   dplyr::mutate(padj=p.adjust(p.value, method = "BH")) %>%
#   dplyr::mutate(score=-log10(padj))
# 
# res_met <- lapply(met_sign %>% {names(.)=.;.}, function(y){
#   print(y)
#   gene_set <- conc %>% 
#     dplyr::filter(n_dataset>7) %>%
#     dplyr::filter(padj<0.01) %>%
#     dplyr::filter(metabolite==y) %>%
#     dplyr::select(gene) %>%
#     distinct() %>%
#     dplyr::pull(gene)
#   
#   lapply(pw_filtered %>% names %>% {names(.)=.;.}, function(x){
#     ft <- data.frame(gene=gene_all) %>%
#       dplyr::mutate(in_pw=ifelse(gene %in% pw_filtered[[x]],1,0) %>% factor(levels=c(1,0))) %>%
#       dplyr::mutate(in_dt=ifelse(gene %in% gene_set,1,0) %>% factor(levels=c(1,0))) %>% 
#       dplyr::select(-gene) %>% 
#       table %>% 
#       fisher.test(alternative = "greater")
#     
#     data.frame(metabolite=y,
#                pathway=x,
#                odds.ratio=ft$estimate %>% as.numeric(),
#                p.value=ft$p.value)
#     
#   }) %>% {do.call(rbind,.)} %>% as.data.frame 
# }) %>% {do.call(rbind,.)} %>% as.data.frame %>%
#   dplyr::left_join(kegg_hierarchy, by=c("pathway"="match_name")) %>%
#   dplyr::filter(Class != "Human Diseases") %>%
#   dplyr::mutate(padj=p.adjust(p.value, method = "BH")) %>%
#   dplyr::mutate(score=-log10(padj))
# 
# # save results to file
# save(res, res_met, file="Files2Share/Results_FisherPathwayEnrichment.Rdata")

#### Load Results ----

load(file="Files2Share/Results_FisherPathwayEnrichment.Rdata")

#### Plot Results ----

# remove metabolites with no significant results
metlist <- res_met %>%
  dplyr::filter(padj<0.01) %>%
  dplyr::select(metabolite) %>%
  table(dnn="metabolite") %>% as.data.frame %>%
  dplyr::filter(Freq>=1) 

# remove pathways with no hit
pwlist <- res_met %>%
  dplyr::filter(padj<0.01) %>%
  dplyr::select(pathway) %>%
  table(dnn="pathway") %>% as.data.frame %>%
  dplyr::filter(Freq>=1) 

annot_row <- metlist %>% 
  dplyr::select(-Freq) %>%
  dplyr::left_join(anno_data %>% dplyr::select(metabolite, H_SUPER_PATHWAY), 
                   by="metabolite") %>%
  dplyr::rename(SUPER_PATHWAY=H_SUPER_PATHWAY) %>%
  # add number of connections
  dplyr::left_join(conc %>% 
                     dplyr::filter(n_dataset>7) %>%
                     dplyr::filter(padj<0.01) %>%
                     dplyr::select(metabolite) %>% 
                     table(dnn="metabolite") %>% as.data.frame, by="metabolite") %>%
  dplyr::rename(degree=Freq) %>%
  # dplyr::mutate(degreeLog=log10(degree)) %>%
  # dplyr::select(-degree) %>%
  tibble::column_to_rownames("metabolite")

# create df for heatmap
phmat <- res_met %>%
  dplyr::filter(metabolite %in% metlist$metabolite) %>%
  dplyr::filter(pathway %in% pwlist$pathway) %>%
  dplyr::select(pathway, metabolite, score) %>%
  tidyr::spread(key = metabolite, value=score) %>% 
  dplyr::arrange(pathway) %>%
  tibble::column_to_rownames("pathway")

annot_col <- res_met %>%
  dplyr::mutate(pw_class=ifelse(Class=="Metabolism", "Metabolism",
                                ifelse(Group=="Immune system","Immune system", "Other"))) %>%
  dplyr::select(pathway,pw_class) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("pathway")

ph <- phmat %>% t %>% as.data.frame %>% 
  pheatmap(cluster_rows = T, cluster_cols = T, clustering_method = "ward.D2",
           cutree_rows = 4, cutree_cols = 7, show_colnames = F,
           annotation_row = annot_row,
           annotation_col = annot_col,
           annotation_colors = list(pw_class=c(Metabolism="thistle",
                                               "Immune system"="olivedrab3",
                                               Other="wheat"),
                                    degree=colorRampPalette(c("#E5E5E5","black"))(10)),
           breaks = c(0,2,5,10,15,20,25),
           color = c("#E5E5E5","#D4A4A4","#D4A4A4","#CB8383","#C36363","#BA4242","#B22222"))

# save to file
pdf(sprintf("%s/FisherEnrichment_Metabolite_Heatmap.pdf", Sys.Date()), width=12, height = 8)
print(ph)
dev.off()

# save pathway annotations to file
# save results to file
library(xlsx)
res_met %>%
  dplyr::select(Class,Group,ID,Name) %>% 
  distinct %>%
  write.xlsx2(file=sprintf("%s/TableS18.xlsx", Sys.Date()), sheetName = "KEGGpathways", 
            col.names = TRUE, row.names = F, append = FALSE)

#### Overview Barplot ----

# extract clusters and get annotations for Fisher's test
z <- data.frame(Cluster=sort(cutree(ph$tree_row, k=2))) %>% 
  tibble::rownames_to_column(var="metabolite")

overview_pea <- rbind(phmat %>% as.data.frame %>%
  dplyr::select(any_of(z %>% dplyr::filter(Cluster==2) %>% dplyr::pull(metabolite))) %>% 
  apply(1,function(x) sum(x > 2)) %>% as.data.frame %>% 
  dplyr::rename(Freq=".") %>%
  dplyr::arrange(-Freq) %>%
  dplyr::mutate(fraction_significant=Freq/nrow(z %>% dplyr::filter(Cluster==2))) %>%
  dplyr::mutate(Cluster=2) %>%
  tibble::rownames_to_column("pathway"),
  phmat %>% as.data.frame %>%
    dplyr::select(any_of(z %>% dplyr::filter(Cluster==1) %>% dplyr::pull(metabolite))) %>% 
    apply(1,function(x) sum(x > 2)) %>% as.data.frame %>% 
    dplyr::rename(Freq=".") %>%
    dplyr::arrange(-Freq) %>%
    dplyr::mutate(fraction_significant=Freq/nrow(z %>% dplyr::filter(Cluster==1))) %>%
    dplyr::mutate(Cluster=1) %>%
    tibble::rownames_to_column("pathway")) %>% 
  dplyr::left_join(annot_col %>% tibble::rownames_to_column("pathway"), by="pathway") %>% 
  dplyr::mutate(Cluster=factor(Cluster, levels=c(2,1))) %>% 
  dplyr::mutate(pw_class=factor(pw_class, levels=c("Immune system","Metabolism","Other"))) %>%
  ggplot(aes(x=reorder(pathway, fraction_significant), y=fraction_significant, group=Cluster, fill=pw_class, color=Cluster)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c(Metabolism="thistle",
                      "Immune system"="olivedrab3",
                      Other="wheat")) +
  scale_color_manual(values=c("1"="gray40","2"="firebrick")) +
  theme_bw() +
  ylab("Fraction of metabolites that display significant enrichment") +
  xlab("pathway") +
  coord_flip()

# save to file
pdf(sprintf("%s/FisherEnrichment_Overview.pdf", Sys.Date()), width=12, height = 8)
print(overview_pea)
dev.off()

#### Degree dataframe for ed -----

degree <- list()
degree$genes <- conc %>%
  dplyr::filter(n_dataset>7) %>%
  dplyr::mutate(sign_gmi=ifelse(padj<0.01, 1,0)) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(n_GMIs=sum(sign_gmi)) %>%
  dplyr::arrange(-n_GMIs)
degree$metabolites <- conc %>%
  dplyr::filter(n_dataset>7) %>%
  dplyr::mutate(sign_gmi=ifelse(padj<0.01, 1,0)) %>%
  dplyr::group_by(metabolite) %>%
  dplyr::summarise(n_GMIs=sum(sign_gmi)) %>%
  dplyr::arrange(-n_GMIs)
save(degree, file = sprintf("%s/DegreeOfGMIs.Rdata", Sys.Date()))
  
  
  
           