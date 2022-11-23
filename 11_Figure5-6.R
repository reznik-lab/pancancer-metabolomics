#### Initialize ----

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load libraries
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggpubr)

#### Load Data ----

# load immune deconvolution analysis results
if(file.exists("results/Workspace_8_ImmuneAnalysis.Rdata")){
  load("results/Workspace_8_ImmuneAnalysis.Rdata")
} else {
  load("data_for_scripts/ImmuneDeconvolution/ImmuneAnalysis.Rdata")
}

#### Figure 5 ----

#Panel a; left barplot
fig5a1 <- ggplot(immunescore_percentage_sig, aes(x = cancer, y = percent, fill = n_samples)) + geom_bar(stat="identity") + 
  ylab("Percentage of Significant Metabolites") + xlab("") +
  theme_classic() + theme(legend.position = "bottom") +
  scale_fill_gradient2(low = "#ffb0b0", high = "#ff0000") + coord_flip() + scale_y_reverse() +
  theme(axis.text = element_text(size = 7),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        text = element_text(size=7),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        axis.ticks=element_line(colour="black"))

#Panel a; middle boxplots
fig5a2 <- ggplot(is_expression, aes(x=expression, y=dataset)) + geom_boxplot(lwd=0.25,outlier.size = 0.1) + theme_classic() +
  theme(axis.text = element_text(size = 7),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        text = element_text(size=7),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        axis.ticks=element_line(colour="black")) +
  ylab("")

#Panel a; right scatterplot
imp_metabolites <- c("quinolinate","kynurenine","nicotinamide ribonucleotide (NMN)",
                     "cholate","creatinine","pyruvate","citrate") #metabolites we want to label in the plot
imp_metabolites <- c(imp_metabolites, HCC_metabolites, ICC_metabolites)

fig5a3 <- ggplot(immunescore_concordance_values, aes(x=concordance,y=dataset)) + geom_point(aes(color=concordance_new), size=0.45) + 
  scale_color_gradient2(low = "blue",mid="grey",high="red",midpoint=0) + theme_bw() +
  xlab("Concordance") + ylab("") + scale_x_continuous(breaks=c(-1,-0.5,0,0.5,1)) + expand_limits(x=c(-1,1)) +
  geom_text_repel(data=subset(immunescore_concordance_values,metabolite %in% imp_metabolites),aes(label = metabolite),
                  color = "black", min.segment.length = 0, box.padding = 0.2, max.overlaps = Inf,
                  segment.size  = 0.2, segment.color = "black", size=1.5) + 
  theme(axis.text = element_text(size = 7),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        text = element_text(size=7),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        axis.ticks=element_line(colour="black"),
        panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        plot.margin= margin(10, 10, 5, 5, "pt"),
        plot.title = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
        plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
        legend.key.size = unit(3, 'mm'), #change legend key size
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7), #change legend text font size
        legend.position="bottom")

#Put panel a together using patchwork
fig5a <- (fig5a1|fig5a2|fig5a3) + plot_layout(widths = c(1, 0.45, 2.5))
ggsave("results/Figure5A.pdf", fig5a, height = 4, width = 6)


#Panel b
#Plot each HCC and ICC sample's ImmuneScore expression
HCC_hm <- ggplot(HCC_IS, aes(x=sample, y=sig, fill = IS)) + geom_tile() + scale_fill_gradient(low = "white",high = "red") + theme_classic() +
  xlab("") + ylab("") + theme(axis.text.x=element_blank(),axis.text.y = element_blank(),legend.position="none")
ICC_is_hm <- ggplot(ICC_IS, aes(x=sample, y=sig, fill = IS)) + geom_tile() + scale_fill_gradient(low = "white",high = "red") + theme_classic() +
  xlab("") + ylab("") + theme(axis.text.x=element_blank(),axis.text.y = element_blank(),legend.position="none")

#Plot metabolite abundance of each sample
HCCis_mat <- ggplot(HCC_matrix_melt,aes(x=variable,y=metabolite,fill=as.numeric(value))) + geom_tile() + 
  scale_fill_gradient2(low = "blue",mid="white",high = "red", midpoint = 0) +
  theme_classic() + xlab("") + ylab("") + theme(axis.text.x=element_blank(),axis.text.y = element_blank(), legend.position = "bottom")
ICCis_mat <- ggplot(ICC_matrix_melt,aes(x=variable,y=metabolite,fill=as.numeric(value))) + geom_tile() + 
  scale_fill_gradient2(low = "blue",mid="white",high = "red", midpoint = 0) + 
  theme_classic() + xlab("") + ylab("") + theme(axis.text.x=element_blank(),axis.text.y = element_blank(), legend.position = "bottom") 

#Put together panel using patchwork
fig5b <- (HCC_hm|ICC_is_hm)/(HCCis_mat|ICCis_mat)
ggsave("results/Figure5B.pdf", fig5b, height = 4, width = 6)


#Panel c
fig5c <- ggplot(immunescore_concordance[1:23,], aes(x = -log10(p.adj), y = metabolite)) + geom_bar(stat="identity") + theme_classic() +
  xlab("-Log10 FDR-Corrected P-Value") + ylab("") + ggtitle("ImmuneScore Concordance")
ggsave("results/Figure5C.pdf", fig5c, height = 3.5, width = 6)


#Panel d
fig5d1 <- ggplot(quinolinate_is_exp, aes(x=Quinolinate, y=ImmuneScore, color = CancerType)) + geom_smooth(method=lm, se=FALSE, fullrange=TRUE,size=0.5) + 
  theme_classic() + geom_point(size=0.25) + theme(legend.position = "bottom") +
  scale_color_manual(values=mapping_color$Color[match(levels(quinolinate_is_exp$CancerType),mapping_color$Name)])

fig5d2 <- ggplot(nmn_is_exp, aes(x=NMN, y=ImmuneScore, color = CancerType)) + geom_smooth(method=lm, se=FALSE, fullrange=TRUE,size=0.5) + 
  theme_classic() + geom_point(size=0.25) + theme(legend.position = "bottom") +
  scale_color_manual(values=mapping_color$Color[match(levels(nmn_is_exp$CancerType),mapping_color$Name)])

#Put together panel using ggpubr
fig5d <- ggarrange(fig5d1,fig5d2, nrow=2, ncol=1, common.legend = T)
ggsave("results/Figure5D.pdf", fig5d, width = 4, height = 8, units = "in")


#Panel e
fig5e <- ggplot(is_pathway_analysis, aes(x=-log10(p.adj), y=pathway)) + 
  geom_bar(stat="identity") + theme_classic() + ylab("") + xlab("-log10(Adjusted P-Value)") + 
  ggtitle("ImmuneScore")  + geom_vline(xintercept = -log10(0.05), color = "red", linetype ="dashed")
ggsave("results/Figure5E.pdf", fig5e, width = 5, height = 3.65, units = "in")


#Panel f
#Flow-sorted data boxplots
fig5f1 <- ggplot(kilgour_nad, aes(x = expression,y=cell_type)) +geom_boxplot(outlier.shape = NA) + theme_classic() + ylab("") + xlab("NAD+") + 
  geom_jitter(color="black", size=0.2, alpha=0.9) 

#cAMP ovarian data scatterplot
fig5f2 <- ggplot(ov_nad, aes(x=NAD,y=ImmuneScore)) + geom_point(size=0.5) + geom_smooth(method = "lm", se = FALSE) + theme_classic() +
  xlab("NAD+")

fig5f <- fig5f1/fig5f2
ggsave("results/Figure5F.pdf", fig5f, width = 2.55, height = 3.85, units = "in")


#### Figure 6 ----

#Panel a
fig6a1 <- ggplot(bindea_concordance, aes(x=concordance,y=-log10(p.adj))) + geom_point(aes(color = color), size = 1) + 
  scale_color_manual(values= c("#D9D9D9","#FBB4AE")) + theme_classic() +
  geom_text_repel(data = subset(bindea_concordance,pair %in% label_pairs_bindea), aes(label = pair, size = 4)) + 
  geom_hline(yintercept = -log10(0.05)) + xlim(-0.35,0.35) + theme(legend.position="none")
fig6a2 <- ggplot(bindea_concordance,aes(x=concordance,y=col2)) + geom_point(shape=124,size=3,color=bindea_concordance$col3) + 
  theme_classic() + ylab("") + xlab("") + theme(axis.ticks = element_blank(),axis.text = element_blank()) + xlim(-0.35,0.35)
ggsave("results/Figure6A.pdf", ggarrange(fig6a1,fig6a2,nrow=2, heights = c(5,1)), width = 3, height = 4, units = "in")

#Panel b
fig6b <- ggplot(bindea_concordance_dots, aes(x=immune_signature,y=-log10(p.adj))) + geom_point(aes(color = color), size = 0.55) + 
  scale_color_manual(values= c("#D9D9D9","#FBB4AE")) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") + 
  xlab("") + ylab("-Log10 Adjusted P-Value")
ggsave("results/Figure6B.pdf", fig6b, height = 3.5, width = 2.75)

#Panel c
fig6c <- ggplot(histamine_mc_exp, aes(x=Histamine, y=Mast.Cells, color = CancerType)) + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + 
  theme_classic() + xlab("Histamine") + ylab("Mast Cells") + geom_point(size=0.2) +
  scale_color_manual(values=mapping_color$Color[match(sort(unique(histamine_mc_exp$CancerType)),mapping_color$Name)])
ggsave("results/Figure6C.pdf", ggarrange(fig6c, common.legend = T), width = 4, height = 4, units = "in")

#Panel d
fig6d <- ggplot(mastcell_concordance[1:21,], aes(x = -log10(p.adj), y = metabolite)) + 
  geom_bar(stat="identity") + theme_classic() +
  xlab("-Log10 FDR-Corrected P-Value") + ylab("") + ggtitle("Mast Cell Concordance") + 
  geom_vline(xintercept = -log10(0.05), color = "red", linetype ="dashed")
ggsave("results/Figure6D.pdf", fig6d, height = 2.5, width = 5)

#Panel e
fig6e <- ggplot(histamine_hdc_exp, aes(x=Histamine, y=HDC, color = CancerType)) + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + 
  theme_classic() + xlab("Histamine") + ylab("HDC") + geom_point(size=0.2) +
  scale_color_manual(values=mapping_color$Color[match(sort(unique(histamine_hdc_exp$CancerType)),mapping_color$Name)])
ggsave("results/Figure6E.pdf", ggarrange(fig6e, common.legend = T), width = 4, height = 4, units = "in")

#Panel f
fig6f <- ggplot(adc_kyn_exp, aes(x=Kynurenine, y=aDC, color = CancerType)) + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + 
  theme_classic() + xlab("Kynurenine") + ylab("aDC") + geom_point(size=0.2) +
  scale_color_manual(values=mapping_color$Color[match(sort(unique(adc_kyn_exp$CancerType)),mapping_color$Name)])
ggsave("results/Figure6F.pdf", ggarrange(fig6f, common.legend = T), width = 4, height = 4, units = "in")

#### Figure S2 ----

supp2 <- ggplot(histamine_concordance, aes(x=pair, y=value)) + geom_boxplot() + 
  theme_classic() + xlab("") + ylab("Concordance") + geom_jitter(shape=16, position=position_jitter(0.2), aes(colour=variable)) +
  scale_color_manual(values=mapping_color$Color[match(sort(unique(histamine_concordance$variable)),mapping_color$Name)]) +
  theme(axis.text.x = element_text(angle = 345, hjust = 0))
ggsave("results/FigureS2.pdf", ggarrange(supp2, common.legend = T), width = 4, height = 4, units = "in")
