#### Initialize ----

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# clear workspace
rm(list = ls())

library(data.table)
library(tidyr)
library(plyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(edgeR)
library(limma)

#########
# differential gene expression test
#########

# define cancer types with both tumor and normal data
selectedDataSetList<-c("BRCA1",
                       "GBM",
                       "COAD",
                       "PDAC",
                       "PRAD",
                       "ccRCC3",
                       "ccRCC4")

setWithRNAseq<-c("BRCA2",
                 "HurthleCC",
                 "GBM",
                 "DLBCL",
                 "ccRCC1",
                 "ccRCC2",
                 "ccRCC3",
                 "ccRCC4")

for(idx in 1:length(selectedDataSetList)){
  
  selectedDataSet<-selectedDataSetList[idx]
  
  message(sprintf("Work on differential gene expression test of %s",selectedDataSet))
  
  withRNAseqData<-FALSE
  withMicroArrayData<-FALSE
  
  if(selectedDataSet %in% setWithRNAseq ){ withRNAseqData<-TRUE }
  if(!(selectedDataSet %in% setWithRNAseq) ){ withMicroArrayData<-TRUE }
  
  ####
  
  filePath<-"data"
  fileName<-"MasterMapping_MetImmune_03_16_2022_release.csv"
  
  fileName<-file.path(filePath,fileName)
  
  masterMappingDat<-read.table(fileName,sep=",",header=TRUE,stringsAsFactors = FALSE)
  
  masterMappingDatSubset<-masterMappingDat[masterMappingDat$Dataset %in% selectedDataSet,]
  
  rownames(masterMappingDatSubset)<-masterMappingDatSubset$CommonID
  
  tumor_commonID<-rownames(masterMappingDatSubset[masterMappingDatSubset$TN %in% "Tumor",])
  tumor_ITHID<-masterMappingDatSubset[tumor_commonID,]$ITHID
  tumor_MetabID<-masterMappingDatSubset[tumor_commonID,]$MetabID
  tumor_RNAID<-masterMappingDatSubset[tumor_commonID,]$RNAID

  ####
  
    normal_commonID<-rownames(masterMappingDatSubset[masterMappingDatSubset$TN %in% "Normal",])
    normal_ITHID<-masterMappingDatSubset[normal_commonID,]$ITHID
    normal_MetabID<-masterMappingDatSubset[normal_commonID,]$MetabID
    normal_RNAID<-masterMappingDatSubset[normal_commonID,]$RNAID
  
  ######
  # load gene expression
  ######
  
  filePath<-"data/transcriptomics_processed"
  fileName<-masterMappingDatSubset$RNAFile[1]
  
  if(withRNAseqData){
    
    # use raw counts for gene expressioin analysis
    fileName<-gsub("tpm","raw",fileName)
    
  }
  
  fileName<-file.path(filePath,fileName)
  
  #tmpFrame<-read.table(fileName,sep=",",header=TRUE,check.names = FALSE,stringsAsFactors = FALSE)
  tmpFrame<-fread(fileName,sep=",",data.table = FALSE)
  
  rownames(tmpFrame)<-tmpFrame[,1]
  expressionData<-as.matrix(tmpFrame[,c(2:ncol(tmpFrame))])
  
  #####
  
  # remove bottom 5% expressed genes 
  rowMean_threshold<-quantile(rowMeans(expressionData),0.05)
  expressionData_trim = expressionData[which(rowMeans(expressionData) > rowMean_threshold),]
  
  ######
  
  expressionData_tumor_normal<-expressionData_trim[,c(tumor_RNAID,normal_RNAID)]
  
  ######
  # Run limma
  
  ################################################################
  # Now compare tumor to normal
  
  genotype2 = rep('Tumor',length(c(tumor_RNAID,normal_RNAID)))
  names(genotype2) = c(tumor_RNAID,normal_RNAID)
  genotype2[normal_RNAID] = 'Normal'
  genotype2 = factor(genotype2,levels=c('Normal','Tumor'))
  
  design2 = model.matrix(~0 + genotype2)
  rownames(design2) = c(tumor_RNAID,normal_RNAID)
  
  if(withRNAseqData){
    
    # Make voom
    dge2 = DGEList(counts = expressionData_tumor_normal)
    dge2 = calcNormFactors(dge2)
    
    v_tvn = voom(dge2,design2,plot = TRUE)
    
    fit_tvn = lmFit(v_tvn,design2)
    head(fit_tvn$coefficient,3)
    
    #fit_tvn = eBayes(fit_tvn)
    
    M_tvn = makeContrasts(status = genotype2Tumor - genotype2Normal, 
                          levels = design2)
    
    fit2 = contrasts.fit(fit_tvn,contrasts = M_tvn)
    head(fit2$coefficient,3)
    
    efit = eBayes(fit2)
    
    plotSA(efit, main="Final model: Mean-variance tren")
    
    volcanoplot(efit)
    
  }
  
  if(withMicroArrayData){ 
    
    fit = lmFit(expressionData_tumor_normal,design2)
    
    head(fit$coefficient,3)
    
    #fit_tvn = eBayes(fit_tvn)
    
    M_tvn = makeContrasts(status = genotype2Tumor - genotype2Normal, 
                          levels = design2)
    
    fit2 = contrasts.fit(fit,contrasts = M_tvn)
    head(fit2$coefficient,3)
    
    efit = eBayes(fit2)
    
    plotSA(efit, main="Final model: Mean-variance tren")
    
    volcanoplot(efit)
  }
  
  top_genes<-topTable(efit,number=nrow(efit))
  
  results<-decideTests(fit2)
  summary(results)
  
  filePath<-file.path("results/Workspace_7_TumorVsNormalAnalysis/diffGeneExpression",selectedDataSet)
  dir.create(filePath,recursive = TRUE)
  fileName<-paste(selectedDataSet,"_DGE_limma_tumor_vs_normal.txt",sep="")
  fileName<-file.path(filePath,fileName)
  
  output<-top_genes
  output$geneName<-rownames(top_genes)
  output<-output[,c("geneName","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
  
  write.table(output,sep="\t",file=fileName,quote=FALSE,row.names = FALSE,col.names = TRUE)
  
}

##########
# differential abundance test
##########

# # select cancer types with both have tumor and normal data
# selectedDataSetList<-c("BRCA1",
#                        "GBM",
#                        "COAD",
#                        "PDAC",
#                        "PRAD",
#                        "ccRCC3",
#                        "ccRCC4")

for(idx in 1:length(selectedDataSetList)){
  
  selectedDataSet<-selectedDataSetList[idx]
  
  message(sprintf("Work on differential abundance of %s",selectedDataSet))
  
  ####
  
  filePath<-"data"
  fileName<-"MasterMapping_MetImmune_03_16_2022_release.csv"
  
  fileName<-file.path(filePath,fileName)
  
  masterMappingDat<-read.table(fileName,sep=",",header=TRUE,stringsAsFactors = FALSE)
  
  masterMappingDatSubset<-masterMappingDat[masterMappingDat$Dataset %in% selectedDataSet,]
  
  rownames(masterMappingDatSubset)<-masterMappingDatSubset$CommonID
  
  tumor_commonID<-rownames(masterMappingDatSubset[masterMappingDatSubset$TN %in% "Tumor",])
  tumor_ITHID<-masterMappingDatSubset[tumor_commonID,]$ITHID
  tumor_MetabID<-masterMappingDatSubset[tumor_commonID,]$MetabID
  
  ####
  
    normal_commonID<-rownames(masterMappingDatSubset[masterMappingDatSubset$TN %in% "Normal",])
    normal_ITHID<-masterMappingDatSubset[normal_commonID,]$ITHID
    normal_MetabID<-masterMappingDatSubset[normal_commonID,]$MetabID
    
  ######
  # load metabolite abundance
  ######
  # handling tumor samples
  ######
    
  filePath<-"data/metabolomics_processed"
  tmpStr<-masterMappingDatSubset[masterMappingDatSubset$TN %in% "Tumor",]$MetabFile[1]
  fileName<-tmpStr
  
  fileName<-file.path(filePath,fileName)
  
  sheetName<-masterMappingDatSubset[masterMappingDatSubset$TN %in% "Tumor",]$MetabFile_sheet[1]
  
  tmpFrame<-read_xlsx(fileName,sheet = sheetName)
  tmpFrame<-as.data.frame(tmpFrame)
  rownames(tmpFrame)<-tmpFrame[,1]
  log2metaboliteData<-as.matrix(tmpFrame[,c(2:ncol(tmpFrame))])
  
  if(selectedDataSet %in% "GBM"){
    log2metaboliteData<-log2metaboliteData[!grepl("Unknown",rownames(log2metaboliteData)),]
  }
  
  log2metaboliteData_tumor<-log2metaboliteData[,tumor_MetabID]
  colnames(log2metaboliteData_tumor)<-tumor_commonID
  
  data_subset<-log2metaboliteData_tumor
  sdValue<-apply(data_subset,1,sd, na.rm=TRUE)
  removedNames<-names(sdValue[ sdValue==0 | is.na(sdValue) ])
  data_subset<-data_subset[!(rownames(data_subset) %in% removedNames),]
  
  log2metaboliteData_tumor<-data_subset
  
  #####
  # handling normal samples
  #####
  
  filePath<-"data/metabolomics_processed"
  tmpStr<-masterMappingDatSubset[masterMappingDatSubset$TN %in% "Normal",]$MetabFile[1]
  fileName<-tmpStr
  
  fileName<-file.path(filePath,fileName)
  
  sheetName<-masterMappingDatSubset[masterMappingDatSubset$TN %in% "Normal",]$MetabFile_sheet[1]
  
  tmpFrame<-read_xlsx(fileName,sheet = sheetName)
  tmpFrame<-as.data.frame(tmpFrame)
  rownames(tmpFrame)<-tmpFrame[,1]
  log2metaboliteData<-as.matrix(tmpFrame[,c(2:ncol(tmpFrame))])
    
    if(selectedDataSet %in% "GBM"){
      log2metaboliteData<-log2metaboliteData[!grepl("Unknown",rownames(log2metaboliteData)),]
    }
    
    log2metaboliteData_normal<-log2metaboliteData[,normal_MetabID]

    colnames(log2metaboliteData_normal)<-normal_commonID
    
    data_subset<-log2metaboliteData_normal
    sdValue<-apply(data_subset,1,sd, na.rm=TRUE)
    removedNames<-names(sdValue[ sdValue==0 | is.na(sdValue) ])
    data_subset<-data_subset[!(rownames(data_subset) %in% removedNames),]
    
    log2metaboliteData_normal<-data_subset
  
  ######
  
  if(TRUE){
    
    common_metabolite<-intersect(rownames(log2metaboliteData_tumor),rownames(log2metaboliteData_normal))
    
    log2metaboliteData_tumor_normal<-cbind(log2metaboliteData_tumor[common_metabolite,],
                                           log2metaboliteData_normal[common_metabolite,])
    
    data_subset<-log2metaboliteData_tumor_normal
    sdValue<-apply(data_subset,1,sd, na.rm=TRUE)
    removedNames<-names(sdValue[ sdValue==0 | is.na(sdValue) ])
    message(sprintf("Removed %s metabolite with no variation",length(removedNames)))
    data_subset<-data_subset[!(rownames(data_subset) %in% removedNames),]
    
    log2metaboliteData_tumor_normal<-data_subset
    
    metaboliteData_tumor_normal<-2^(log2metaboliteData_tumor_normal)
    
    metaboliteDat<-metaboliteData_tumor_normal
    
  }
  
  ######
  # differential metabolite abundance test function
  ######
    
    diffMetaboliteAbundance = function(x,ix1,ix2,outputfilePath = NULL,makeBoxPlot=FALSE,ressort = TRUE){
      # Ed Reznik wrote the initial version
      # do standard differential abundance analysis.
      # x is a matrix, metabolites x samples
      # ix1,ix2 are column names for the two groups to compare. WE COMPARE GROUP 2 TO GROUP 1!
      
      # do differential abundance with mann-whitney tests, calculate p-value, FDR correct, return
      
      res = data.frame()
      
      for (rname in rownames(x)){
        
        cat(sprintf("metabolite name: %s\n",rname))
        
        temp1 = as.matrix( x[rname,(colnames(x) %in% ix1)] )
        temp2 = as.matrix( x[rname,(colnames(x) %in% ix2)] )
        
        meanValueGroup1<-mean(temp1)
        meanValueGroup2<-mean(temp2)
        
        if( !is.na(meanValueGroup1) && !is.na(meanValueGroup2) ){
          
          res[rname,'log2FC'] = log2( meanValueGroup2/meanValueGroup1 )
          
          temptest = wilcox.test(temp1,temp2,alternative = "two.sided", exact=FALSE)
          
          res[rname,'pvalue'] = temptest$p.value
          res[rname,'UStat'] = temptest$statistic
          
          ####
          if(makeBoxPlot){
            
            frame1<-data.frame(temp1,rep("WT",nrow(temp1)),stringsAsFactors = FALSE)
            colnames(frame1)<-c("abundance","group")
            
            frame2<-data.frame(temp2,rep("MUT",nrow(temp2)),stringsAsFactors = FALSE)
            colnames(frame2)<-c("abundance","group")
            
            dataCombined<-rbind(frame1,frame2)
            metabolite<-rname
            pvalue2<-round(temptest$p.value,5)
            
            p <- ggplot(dataCombined,aes(x=group,y=log2(abundance),fill=group))
            p <- p + ggtitle(paste(metabolite,", raw p-value: ",pvalue2,sep=""))
            p <- p + geom_boxplot(outlier.shape = NA)
            #p <- p + geom_violin()
            p <- p + geom_quasirandom()
            p <- p + theme_classic()
            #p <- p + theme_bw()
            p <- p + theme(
              #axis.text.x = element_text(angle=-45, hjust=0, vjust=1), 
              #plot.margin = unit(c(1.5, 1, 1, 1), "cm"), 
              plot.title = element_text(size = 10, face = "bold", colour = "black", vjust = 0, hjust=0.5),
              plot.subtitle = element_text(hjust=0.5)
            )
            
            #print(p)
            
            condensedName<-gsub(" ","",rname)
            condensedName<-gsub("/","#",condensedName)
            #condensedName<-gsub(":","",condensedName)
            
            fileName<-paste(condensedName,".pdf",sep="")
            fileName<-file.path(outputfilePath,fileName)
            
            pdf(fileName,width=6.51,height=6.36)
            print(p)
            dev.off()
            #p<- p + ggsave(fileName)
          }  
          
          ####
          
        }else{
          res[rname,'log2FC'] = NA
          res[rname,'pvalue'] = NA
          res[rname,'UStat'] = NA
        }
        
      }
      
      resFiltered<-res[is.na(res$log2FC),]
      res<-res[!is.na(res$log2FC),]
      
      # remove the situation that pvalue is not available
      res<-res[!is.na(res$pvalue),]
      
      res$padj = p.adjust(res$pvalue,method = 'BH')
      
      if (ressort){
        res = res[order(res$pvalue,decreasing = FALSE),,drop = FALSE]
      }
      
      return(res)
    }  
  
  ######
  
  # log2( mean(ix2) / mean (ix1) ) fold change
  results_tumor_vs_normal_wilcox_test<-diffMetaboliteAbundance(x=metaboliteDat,ix1=normal_commonID,ix2=tumor_commonID)
  
  ####
  res<-results_tumor_vs_normal_wilcox_test
  res$name<-rownames(res)
  res$cancerType<-selectedDataSet
  res<-res[,c("cancerType","name","log2FC","UStat","pvalue","padj")]
  
  ###
  
  filePath<-file.path("results/Workspace_7_TumorVsNormalAnalysis/diffMetaboliteAbundance",selectedDataSet)
  dir.create(filePath,recursive = TRUE)
  fileName<-paste(selectedDataSet,"_differential_metabolite_abundance_wilcox_test_tumor_vs_normal.txt",sep="")
  fileName<-file.path(filePath,fileName)
  
  write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)
  
}

##########
# start of making aggregated summary
##########

# set metabolite and gene threshold, assuming we are using padj column
pthresh_gene = 0.05
pthresh_met = 0.05

######
# load KEGG pathway data - metabolite centric
######

fileName<-paste("data_for_scripts/TumorVsNormalAnalysis/KEGG_metabolic_pathway_metaboliteInfoAll.txt",sep="")

metPathwayList<-fread(file=fileName,header=TRUE,sep="\t",stringsAsFactors = FALSE)

#######
# load KEGG pathway data - gene centric
#######

fileName<-paste("data_for_scripts/TumorVsNormalAnalysis/KEGG_metabolic_pathway_geneInfoAll.txt",sep="")

genePathwayList<-fread(file=fileName,header=TRUE,sep="\t",stringsAsFactors = FALSE)

pathwayList<-unique(metPathwayList$pathwayName)

outDat_gene_list<-list()
outDat_met_list<-list()

#####

for(idx_i in 1:length(pathwayList)){
  
  pathwayName2<-pathwayList[idx_i]
  
  message(sprintf("Work on %s",pathwayName2))
  
  # name for file name saving
  prettyname = gsub("/","#", pathwayName2)
  
  # get metabolites and genes
  
  metix = metPathwayList[metPathwayList$pathwayName %in% pathwayName2,]$keggMetaboliteID
  
  metix<-unique(metix)
  
  metix_base<-data.frame(pathway=pathwayName2,
                         H_KEGG=metix,
                         stringsAsFactors = FALSE)
  
  message(sprintf("make %s",prettyname))
  
  geneix = genePathwayList[genePathwayList$pathwayName %in% pathwayName2,]$geneSymbol
  
  geneix_base<-data.frame(pathway=pathwayName2,
                          geneName=geneix,
                          stringsAsFactors = FALSE)
  
  if(pathwayName2 == "Thiamine metabolism"){
    geneix_base<-geneix_base[!(geneix_base$geneName %in% "adenylate kinase isoenzyme 1-like [KO:K00939] [EC:2.7.4.3]"),]
  }
  
  #####
  
  # cancer types both have tumor and normal data
  cancerTypeList<-c("BRCA1","GBM","COAD","PDAC","PRAD","ccRCC3","ccRCC4")
  
  testCaseName<-"tumor_vs_normal"
  
  for(idx_j in 1:length(cancerTypeList)){
    
    cancerTypeName<-cancerTypeList[idx_j]
    
    message(sprintf("Work on %s",cancerTypeName))
    
    ####
    
    filePath<-"data/metabolomics_processed"
    fileName<-paste("PreprocessedData_",cancerTypeName,".xlsx",sep="")
    fileName<-file.path(filePath,fileName)
    
    metinfo<-read_xlsx(path=fileName,sheet="metanno")
    
    metinfoSubset<-metinfo[,c("H_name","H_KEGG")]
    colnames(metinfoSubset)<-c("name","H_KEGG")
    
    ####
    
    filePath<-file.path("results/Workspace_7_TumorVsNormalAnalysis/diffMetaboliteAbundance",cancerTypeName)
    fileName<-paste(cancerTypeName,"_differential_metabolite_abundance_wilcox_test_tumor_vs_normal.txt",sep="")
    
    fileName<-file.path(filePath,fileName)
    
    mdata<-fread(file=fileName,sep="\t",data.table=FALSE)
    rownames(mdata)<-mdata$name
    
    mdata<-merge(mdata,metinfoSubset,by="name",all.x = TRUE)
    rownames(mdata)<-mdata$name
    
    mdata[which(mdata$padj > pthresh_met),'log2FC'] = "non-sig"
    
    
    if(idx_j == 1){
      
      metabolite_name_combined<-merge(metix_base,mdata[,c("name","H_KEGG")],all.x=TRUE,by="H_KEGG")
      metabolite_name_combined$combined_name<-paste(metabolite_name_combined$H_KEGG,metabolite_name_combined$name,sep="_")
      colnames(metabolite_name_combined)<-c("H_KEGG","pathway",paste(cancerTypeName,"name",sep="_"),"combined_name")
      
      metabolite_name_combined<-metabolite_name_combined[,c("combined_name","pathway",paste(cancerTypeName,"name",sep="_"))]
      
      metabolite_log2FC_combined<-merge(metix_base,mdata[,c("log2FC","H_KEGG")],all.x=TRUE,by="H_KEGG")
      metabolite_log2FC_combined$combined_name<-metabolite_name_combined$combined_name
      colnames(metabolite_log2FC_combined)<-c("H_KEGG","pathway",paste(cancerTypeName,"log2FC",sep="_"),"combined_name")
      
      metabolite_log2FC_combined<-metabolite_log2FC_combined[,c("combined_name","pathway",paste(cancerTypeName,"log2FC",sep="_"))]
      
    }else{
      
      tmp_name_combined<-merge(metix_base,mdata[,c("name","H_KEGG")],all.x=TRUE,by="H_KEGG")
      tmp_name_combined<-tmp_name_combined[,c("H_KEGG","name")]
      tmp_name_combined$combined_name<-paste(tmp_name_combined$H_KEGG,tmp_name_combined$name,sep="_")
      
      colnames(tmp_name_combined)<-c("H_KEGG",paste(cancerTypeName,"name",sep="_"),"combined_name")
      tmp_name_combined<-tmp_name_combined[,c("combined_name",paste(cancerTypeName,"name",sep="_"))]
      
      
      metabolite_name_combined<-merge(metabolite_name_combined,tmp_name_combined,by="combined_name",all.x=TRUE,all.y=TRUE)
      
      tmp_log2FC_combined<-merge(metix_base,mdata[,c("log2FC","H_KEGG")],all.x=TRUE,by="H_KEGG")
      tmp_log2FC_combined<-tmp_log2FC_combined[,c("H_KEGG","log2FC")]
      tmp_log2FC_combined$combined_name<-tmp_name_combined$combined_name
      
      colnames(tmp_log2FC_combined)<-c("H_KEGG",paste(cancerTypeName,"log2FC",sep="_"),"combined_name")
      tmp_log2FC_combined<-tmp_log2FC_combined[,c("combined_name",paste(cancerTypeName,"log2FC",sep="_"))]
      
      metabolite_log2FC_combined<-merge(metabolite_log2FC_combined,tmp_log2FC_combined,by="combined_name",all.x=TRUE, all.y=TRUE)
      
      
    }
    
    ####
    # load differential gene expression result
    ####
    
    if(TRUE){
      
      filePath<-file.path("results/Workspace_7_TumorVsNormalAnalysis/diffGeneExpression",cancerTypeName)
      fileName<-paste(cancerTypeName,"_DGE_limma_tumor_vs_normal.txt",sep="")
      fileName<-file.path(filePath,fileName)
      
      dge = fread(file=fileName,data.table = FALSE)
      rownames(dge)<-dge$geneName
      
      dge[which(dge$adj.P.Val > pthresh_gene),'logFC'] = "non-sig"
      
      dge$name<-dge$geneName
      
      #####
      
      if(idx_j == 1){
        
        gene_name_combined<-merge(geneix_base,dge[,c("name","geneName")],all.x=TRUE,by="geneName")
        colnames(gene_name_combined)<-c("geneName","pathway",paste(cancerTypeName,"name",sep="_"))
        
        gene_log2FC_combined<-merge(geneix_base,dge[,c("logFC","geneName")],all.x=TRUE,by="geneName")
        colnames(gene_log2FC_combined)<-c("geneName","pathway",paste(cancerTypeName,"log2FC",sep="_"))
        
        
        
      }else{
        
        tmp_name_combined<-merge(geneix_base,dge[,c("name","geneName")],all.x=TRUE,by="geneName")
        tmp_name_combined<-tmp_name_combined[,c("geneName","name")]
        colnames(tmp_name_combined)<-c("geneName",paste(cancerTypeName,"name",sep="_"))
        gene_name_combined<-merge(gene_name_combined,tmp_name_combined,by="geneName",all.x=TRUE)
        
        
        tmp_log2FC_combined<-merge(geneix_base,dge[,c("logFC","geneName")],all.x=TRUE,by="geneName")
        tmp_log2FC_combined<-tmp_log2FC_combined[,c("geneName","logFC")]
        colnames(tmp_log2FC_combined)<-c("geneName",paste(cancerTypeName,"log2FC",sep="_"))
        gene_log2FC_combined<-merge(gene_log2FC_combined,tmp_log2FC_combined,by="geneName",all.x=TRUE)
        
        
        
      }
      
      
      
    }
  }    
  
  
  #####
  # side by side plot
  #####
  
  gene_log2FC_subdat<-gene_log2FC_combined[,c(1,3:ncol(gene_log2FC_combined))]
  colnames(gene_log2FC_subdat)[1]<-"name"
  gene_log2FC_subdat$omics<-"gene"
  
  metabolite_log2FC_subdat<-metabolite_log2FC_combined[,c(1,3:ncol(metabolite_log2FC_combined))]
  colnames(metabolite_log2FC_subdat)[1]<-"name"
  metabolite_log2FC_subdat$omics<-"metabolite"
  
  ####
  message(sprintf("Work summary on gene %s pathway",pathwayName2))
  
  gene_log2FC_long<-reshape2::melt(gene_log2FC_subdat, id=c("name","omics"))
  gene_log2FC_long$variable<-gsub("_log2FC","",gene_log2FC_long$variable)
  colnames(gene_log2FC_long)[3]<-"dataset"
  
  if(sum(is.na(gene_log2FC_long$value))>0){
    gene_log2FC_long[is.na(gene_log2FC_long$value),]$value<-"not-measured"
  }
  
  datasetList<-unique(gene_log2FC_long$dataset)
  outDat_gene<-lapply(1:length(datasetList), function(x){
    
    tmpDat<-gene_log2FC_long[gene_log2FC_long$dataset %in% datasetList[x],]
    num_Of_total_items_in_this_pathway<-length(geneix_base$geneName)
    valueVector<-as.numeric(tmpDat[ !tmpDat$value %in% c("non-sig","not-measured"),]$value)
    num_Of_items_up<-sum(valueVector>0)
    num_Of_items_down<-sum(valueVector<0)
    num_Of_items_measured<-length(tmpDat[ !tmpDat$value %in% c("not-measured"),]$value)
    
    tmpDat$pathway<-pathwayName2
    tmpDat$numOfTotalItemsInPathway<-num_Of_total_items_in_this_pathway
    tmpDat$numOfItemsMeasured<-num_Of_items_measured
    tmpDat$numOfItemsUp<-num_Of_items_up
    tmpDat$numOfItemsDown<-num_Of_items_down
    
    return(tmpDat)
    
  })
  
  outDat_gene<-rbind.fill(outDat_gene)
  
  outDat_gene_list[[idx_i]]<-outDat_gene
  
  filePath<-"~/work/Ed_lab/pancancer_metabo_immu/results/tumor_vs_normal_aggregated_summary/11_15_2022"
  dir.create(filePath,recursive=TRUE)
  fileName<-paste(prettyname,"_gene_aggregated_summary.txt")
  fileName<-file.path(filePath,fileName)
  
  write.table(outDat_gene,sep="\t",file=fileName,quote=FALSE,row.names = FALSE,col.names = TRUE)
  
  
  #####
  message(sprintf("Work summary on met %s pathway",pathwayName2))
  
  metabolite_log2FC_long<-reshape2::melt(metabolite_log2FC_subdat, id=c("name","omics"))
  metabolite_log2FC_long$variable<-gsub("_log2FC","",metabolite_log2FC_long$variable)
  colnames(metabolite_log2FC_long)[3]<-"dataset"
  
  if(sum(is.na(metabolite_log2FC_long$value))>0 ){
    metabolite_log2FC_long[is.na(metabolite_log2FC_long$value),]$value<-"not-measured"
  }
  
  datasetList<-unique(metabolite_log2FC_long$dataset)
  outDat_met<-lapply(1:length(datasetList), function(x){
    
    tmpDat<-metabolite_log2FC_long[metabolite_log2FC_long$dataset %in% datasetList[x],]
    num_Of_total_items_in_this_pathway<-length(metix_base$H_KEGG)
    valueVector<-as.numeric(tmpDat[ !tmpDat$value %in% c("non-sig","not-measured"),]$value)
    num_Of_items_up<-sum(valueVector>0)
    num_Of_items_down<-sum(valueVector<0)
    num_Of_items_measured<-length(tmpDat[ !tmpDat$value %in% c("not-measured"),]$value)
    
    tmpDat$pathway<-pathwayName2
    tmpDat$numOfTotalItemsInPathway<-num_Of_total_items_in_this_pathway
    tmpDat$numOfItemsMeasured<-num_Of_items_measured
    tmpDat$numOfItemsUp<-num_Of_items_up
    tmpDat$numOfItemsDown<-num_Of_items_down
    
    return(tmpDat)
    
  })
  
  outDat_met<-rbind.fill(outDat_met)
  
  outDat_met_list[[idx_i]]<-outDat_met
  
  filePath<-"results/Workspace_7_tumorVsNormalAnalysis/tumor_vs_normal_aggregated_summary"
  dir.create(filePath,recursive=TRUE)
  fileName<-paste(prettyname,"_metabolite_aggregated_summary.txt")
  fileName<-file.path(filePath,fileName)
  
  write.table(outDat_met,sep="\t",file=fileName,quote=FALSE,row.names = FALSE,col.names = TRUE)
  
  
  
}    

filePath<-"results/Workspace_7_TumorVsNormalAnalysis/tumor_vs_normal_aggregated_summary"
dir.create(filePath,recursive=TRUE)
fileName<-paste("tumor_vs_normal_gene_aggregated_summary.rds")
fileName<-file.path(filePath,fileName)

saveRDS(outDat_gene_list,file=fileName)

filePath<-"results/Workspace_7_TumorVsNormalAnalysis/tumor_vs_normal_aggregated_summary"
dir.create(filePath,recursive=TRUE)
fileName<-paste("tumor_vs_normal_metabolite_aggregated_summary.rds")
fileName<-file.path(filePath,fileName)

saveRDS(outDat_met_list,file=fileName)
