#### Initialize ----

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(dplyr)
library(magrittr)
library(xlsx)
library(ggplot2)
library(ggrepel)
library(writexl)

#### Global parameters ----

mapping_color <- data.frame(OldName= c("BRCA_Terunuma","BRCA_Tang","COAD","CPTAC_GBM","DLBCL",
                                       "HCC","LiCa1","LiCa2","OV","PDAC",
                                       "PRAD","RC12","RC18","RC18_2","RC20"),
                            NewName = c("BRCA1","BRCA2","COAD","GBM","DLBCL",
                                        "HurthleCC","HCC","ICC","OV","PDAC",
                                        "PRAD","ccRCC1","ccRCC2","ccRCC3","ccRCC4")) %>%
  dplyr::arrange(NewName) %>%
  dplyr::mutate(Color = c("#e5b6a5","#c3c1a4","#e5e3c0","#bbe1b9","#88beab",
                          "#97b1ab","#abc5bf","#abe7e1","#c6e6e7","#8fcde1",
                          "#88aee1","#a6b7d3","#d2d5ed","#c8b7de","#e6bbcc"))

#### Load Data ----

# load RDS files
ss_gene <- readRDS("results/Workspace_7_TumorVsNormalAnalysis/tumor_vs_normal_aggregated_summary/tumor_vs_normal_gene_aggregated_summary.rds")
ss_met <- readRDS("results/Workspace_7_TumorVsNormalAnalysis/tumor_vs_normal_aggregated_summary/tumor_vs_normal_metabolite_aggregated_summary.rds")

load("data_for_scripts/TumorVsNormalAnalysis/Annotations_Pathways.Rdata")

#### Arrange Data ----

# assign names to list items
nm <- sapply(ss_gene, function(x){x$pathway %>% unique})
names(ss_gene) <- nm
nm <- sapply(ss_met, function(x){x$pathway %>% unique})
names(ss_met) <- nm

# merge and clean
dt_new <- lapply(ss_gene, function(x){
  x %>% 
    dplyr::select(-name, -value) %>%
    distinct()
}) %>% {do.call(rbind,.)} %>%
  rbind(lapply(ss_met, function(x){
    x %>% 
      dplyr::select(-name, -value) %>%
      distinct()
  }) %>% {do.call(rbind,.)} 
  ) 

# find pathways with no metabolites or gene measured in 2 dataset or more
pw_remove <- dt_new %>%
  # dplyr::filter(omics=="metabolite") %>% 
  dplyr::group_by(pathway, omics) %>% 
  dplyr::summarise(var=sum(numOfItemsMeasured!=0)) %>% 
  dplyr::filter(var<=5) %>%
  dplyr::pull(pathway)

# remove those pathways
dt_new %<>%
  dplyr::filter(!(pathway %in% pw_remove))

# create score per dataset
dt_new %<>%
  dplyr::mutate(score_DF=ifelse(numOfItemsMeasured!=0,(numOfItemsUp+numOfItemsDown)/numOfItemsMeasured,0))

# compute correlation between gene and metabolite scores per pathways
dt_new_cor <- dt_new %>%
  dplyr::select(-numOfTotalItemsInPathway,-numOfItemsMeasured,
                -numOfItemsUp,-numOfItemsDown) %>%
  tidyr::pivot_wider(names_from = "omics", values_from = "score_DF") %>% 
  dplyr::group_by(pathway) %>%
  dplyr::mutate(cor=cor.test(gene,metabolite, method="spearman")$estimate , 
                p.value=cor.test(gene,metabolite, method="spearman")$p.value,
                ci_lower=DescTools::SpearmanRho(gene,metabolite, conf.level=0.95)[2] %>% as.vector,
                ci_upper=DescTools::SpearmanRho(gene,metabolite, conf.level=0.95)[3] %>% as.vector) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(p.value)

# add p-value adjustment
dt_new_cor %<>%
  dplyr::left_join(dt_new_cor %>%
                     dplyr::select(pathway, cor, p.value, ci_lower, ci_upper) %>%
                     distinct %>%
                     dplyr::mutate(padj=p.adjust(p.value, method = "BH"))) %>%
  mutate_at(vars(pathway), funs(factor(.,levels=unique(.)))) 

#### Making heatmap ----

genesum = data.frame()
for (ii in 1:length(ss_gene)){
  temp = ss_gene[[ii]][,c(3,5:9)]
  if (any(duplicated(temp))){
    temp = temp[-which(duplicated(temp)),]
  }
  genesum = rbind(genesum,temp)
}

genesum$DF = (genesum$numOfItemsUp-genesum$numOfItemsDown)/genesum$numOfItemsMeasured
genesum$omics = 'Transcript'

# score each pathway by its mean
pwayscore = c()
for (ii in unique(genesum$pathway)){
  pwayscore[ii] = mean(genesum[which(genesum$pathway == ii),'DF'])
}

metsum = data.frame()
for (ii in 1:length(ss_gene)){
  temp = ss_met[[ii]][,c(3,5:9)]
  if (any(duplicated(temp))){
    temp = temp[-which(duplicated(temp)),]
  }
  metsum = rbind(metsum,temp)
}
metsum$DF = (metsum$numOfItemsUp-metsum$numOfItemsDown)/metsum$numOfItemsMeasured
metsum = metsum[-which(is.na(metsum$DF)),]
metsum$omics = 'Metabolite'

allsum = rbind(metsum,genesum)
allsum$studypway = paste(allsum$dataset,allsum$pathway,sep =':')

# calculate the color of the segments
dfdiff = data.frame(matrix(NA,0,4),stringsAsFactors = FALSE)
for (ii in unique(allsum$studypway)){
  tempg = allsum[which(allsum$studypway == ii & allsum$omics == 'Transcript'),]
  tempm = allsum[which(allsum$studypway == ii & allsum$omics == 'Metabolite'),]
  dfdiff[ii,] = c(strsplit(ii,'\\:')[[1]][1],strsplit(ii,'\\:')[[1]][2],tempg[1,'DF'],tempm[1,'DF'])
}
colnames(dfdiff) = c('dataset','pathway','DFgene','DFmet')
dfdiff$DFgene = as.numeric(dfdiff$DFgene)
dfdiff$DFmet = as.numeric(dfdiff$DFmet)
dfdiff$sign = ifelse(sign(dfdiff$DFgene) == sign(dfdiff$DFmet),1,-1)
dfdiff$x = 'Metabolite'
dfdiff$xend = 'Transcript'

pwayscore2 = c()
for (ii in unique(genesum$pathway)){
  pwayscore2[ii] = sum(dfdiff[which(dfdiff$pathway == ii),'sign'],na.rm = TRUE)
}

pwlist <- (dt_new_cor %>% dplyr::pull(pathway) %>% as.character %>% unique)

# filter pathways
pwayscore2 <- pwayscore2[names(pwayscore2) %in% pwlist]
allsum %<>%
  dplyr::filter(pathway %in% pwlist)

# add score to dataframe
allsum %<>%
  dplyr::left_join(pwayscore2 %>% as.data.frame %>% dplyr::rename(score=".") %>% tibble::rownames_to_column("pathway"),
                   by="pathway") %>%
  # order pathways by score
  dplyr::arrange(score) %>%
  dplyr::mutate(pathway=factor(pathway, levels=pathway %>% unique))

# re-order pathway in the heatmap
pwayscore2 = c()
for (ii in unique(metsum$pathway)){
  pwayscore2[ii] = sum(dfdiff[which(dfdiff$pathway == ii),'sign'],na.rm = TRUE)
}
p2use = unique(allsum$pathway)

# check that dfdiff and allsum have same pathway scores
diffnamestest = setequal(unique(allsum$pathway),unique(dfdiff$pathway))

# order the pathways 
allsum2 = allsum[which(allsum$pathway %in% p2use),]
dfdiff2 = dfdiff[which(dfdiff$pathway %in% p2use),]
allsum2$pathway = factor(allsum2$pathway,levels = names(pwayscore2)[order(pwayscore2)])
dfdiff2$pathway = factor(dfdiff2$pathway,levels = names(pwayscore2)[order(pwayscore2)])

allsum2$dataset<-factor(allsum2$dataset,levels=c("BRCA1","COAD","GBM","PDAC","PRAD","ccRCC3","ccRCC4"))
dfdiff2$dataset<-factor(dfdiff2$dataset,levels=c("BRCA1","COAD","GBM","PDAC","PRAD","ccRCC3","ccRCC4"))

#### Figure 3A ----

pdf("results/Figure3A.pdf" ,height = 7.75,width = 8.5)
print(
  ggplot(allsum2,aes(x=omics, y=pathway,color = DF)) + 
    geom_point(aes(size=log10(numOfItemsMeasured))) + 
    geom_point(aes(size=log10(numOfItemsMeasured)),shape = 21,color = 'black') + 
    scale_size_continuous(range = c(0.1,3)) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    scale_color_gradient2() +
    
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    xlab("") +
    ylab("") +
    geom_segment(data = dfdiff2[which(dfdiff2$sign == 1),],aes(x = x,xend = xend,y = pathway,yend = pathway),color = 'black',linetype = 'solid') +
    geom_segment(data = dfdiff2[which(dfdiff2$sign == -1),],aes(x = x,xend = xend,y = pathway,yend = pathway),color = 'gray',linetype = 'dashed') +
    facet_wrap(~dataset,nrow = 1) +
    theme_minimal() + 
    theme(panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank(),
          text=element_text(size=7), 
          axis.text=element_text(size=7), 
          axis.title=element_text(size=7),
          axis.text.y = element_text(color = "black"),
          axis.text.x = element_text(color="black"),
          legend.title = element_text(size=7), #change legend title font size
          legend.text = element_text(size=7)
    )
) 
dev.off()

##### Figure 3B ----

pdf("results/Figure3B.pdf", width=5, height=3, onefile = F)
print(
  dt_new_cor %>%  
    dplyr::select(pathway,cor,ci_lower,ci_upper,p.value,padj) %>%
    distinct() %>%
    dplyr::left_join(kegg_hierarchy, by=c("pathway"="Name")) %>%
    dplyr::mutate(sign=ifelse(ci_lower>0 | ci_upper<0,"significant","non significant")) %>%
    dplyr::arrange(-cor) %>%
    dplyr::mutate(label=ifelse(sign=="significant", pathway, NA)) %>% 
    mutate_at(vars(pathway), funs(factor(.,levels=unique(.)))) %>% 
    ggplot(aes(x=pathway,y=cor, label=label)) +
    geom_point(aes(color=sign)) +
    geom_hline(yintercept = 0, linetype="solid", color="black") +
    scale_color_manual(values=c(significant="firebrick",`non significant`="gray60")) +
    scale_y_continuous(limits=c(-1.15,1.15),breaks=c(-1,-0.5,0,0.5,1)) +
    geom_label_repel(box.padding = 0.5,
                     size=7 /.pt, 
                     fill="white", 
                     max.overlaps=Inf,
                     nudge_y=0.1,
                     direction="y") +
    xlab("KEGG pathway") +
    ylab("Spearman correlation coefficient") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text=element_text(size=7), 
          axis.text=element_text(size=7), 
          axis.title=element_text(size=7),
          legend.position="none",
          axis.text.x = element_blank(), 
          axis.ticks.x=element_blank(),
          axis.text.y = element_text(color = "black"))
)
dev.off()

##### Figure 3C ----

# plot scores for selected pathways
pw <- "Citrate cycle (TCA cycle)"

pdf("results/Figure3C.pdf", width=2.5, height=2.5, onefile = F)
print(
  dt_new_cor %>%
    dplyr::filter(pathway == pw) %>%
    ggplot(aes(x=metabolite, y=gene, label=dataset)) +
    geom_point(aes(color=dataset),alpha=0.5) +
    geom_smooth(aes(color=ifelse(padj<0.05,"significant","non significant")),method=MASS::rlm, se=F, size=0.2) +
    scale_color_manual(values=setNames(mapping_color$Color, nm=mapping_color$NewName)) +
    geom_text_repel(size=2) +
    theme_bw() +
    xlab("Metabolomics DF score") +
    ylab("Transcriptomics DF score") +
    ggtitle(pw) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=7),
          legend.position="none",
          axis.text.y = element_text(color = "black"),
          axis.text.x = element_text(color="black"))
)
dev.off()

#### Fisher's exact test Enrichment ----

# enrichment dataset
dt_pw_enr <- dt_new_cor %>%
  dplyr::left_join(kegg_hierarchy, by=c("pathway"="Name")) %>%
  dplyr::select(pathway,Group,p.value) %>%
  dplyr::distinct() %>%
  dplyr::mutate(is_sign=ifelse(p.value<0.05,1,0) %>% factor(levels = c(1,0)),
                is_aa=ifelse(Group=="Amino acid metabolism",1,0) %>% factor(levels = c(1,0)),
                is_carb=ifelse(Group=="Carbohydrate metabolism",1,0) %>% factor(levels = c(1,0)))

# Fisher's test for Amino acid metabolism
dt_pw_enr %>%  
  dplyr::select(is_sign,is_aa) %>% 
  table() %>% fisher.test %>% .$p.value

# Fisher's test for Carbohydrate metabolism
dt_pw_enr %>%  
  dplyr::select(is_sign,is_carb) %>% 
  table() %>% fisher.test %>% .$p.value

#### Table S2 ---

# save pathway annotations to file
dt_new %>%
  dplyr::arrange(pathway) %>% 
  write.xlsx2(file=file.path("results/TableS2.xlsx"), sheetName = "Tumor_vs_Normal", 
              col.names = TRUE, row.names = F, append = FALSE)

#### Table S3 ----

TvsN_PathwaySpearman<-dt_new_cor %>%  
  dplyr::select(pathway,cor,ci_lower,ci_upper,p.value,padj) %>%
  distinct() %>%
  dplyr::left_join(kegg_hierarchy, by=c("pathway"="Name")) %>%
  dplyr::mutate(sign=ifelse(ci_lower>0 | ci_upper<0,"significant","non significant")) %>%
  dplyr::arrange(-cor) %>%
  dplyr::mutate(label=ifelse(sign=="significant", pathway, NA)) %>% 
  mutate_at(vars(pathway), funs(factor(.,levels=unique(.))))

write_xlsx(TvsN_PathwaySpearman,
           path="results/TableS3.xlsx")
