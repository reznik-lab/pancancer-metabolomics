#Load libraries
library(ggplot2)
library(ggrepel)
library(patchwork)

#Load Rdata file to create figures
load("immune_deconvolution_results.Rdata")

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
ggsave("Figure 5A ImmuneScore Concordance by Cancer Type.pdf", fig5a, height = 4, width = 6)


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
ggsave("Figure 5b ICC HCC ImmuneScore Concordance", fig5b, height = 4, width = 4)


#Panel c
fig5c <- ggplot(immunescore_concordance[1:23,], aes(x = -log10(p.adj), y = metabolite)) + geom_bar(stat="identity") + theme_classic() +
  xlab("-Log10 FDR-Corrected P-Value") + ylab("") + ggtitle("ImmuneScore Concordance")
ggsave("Figure 5c ImmuneScore Concordance by Metabolite.pdf", fig5c, height = 3.5, width = 3.5)


#Panel d
fig5d1 <- ggplot(quinolinate_is_exp, aes(x=Quinolinate, y=ImmuneScore, color = CancerType)) + geom_smooth(method=lm, se=FALSE, fullrange=TRUE,size=0.5) + 
  theme_classic() + geom_point(size=0.25) + theme(legend.position = "bottom") +
  scale_color_manual(values=mapping_color$Color[match(levels(quinolinate_is_exp$CancerType),mapping_color$Name)])

fig5d2 <- ggplot(nmn_is_exp, aes(x=NMN, y=ImmuneScore, color = CancerType)) + geom_smooth(method=lm, se=FALSE, fullrange=TRUE,size=0.5) + 
  theme_classic() + geom_point(size=0.25) + theme(legend.position = "bottom") +
  scale_color_manual(values=mapping_color$Color[match(levels(nmn_is_exp$CancerType),mapping_color$Name)])

#Put together panel using patchwork
fig5d <- fig5d1/fig5d2
ggsave("Figure 5d ImmuneScore vs Quinolinate and NMN Concordance Scatterplots Across Datasets.pdf", fig5d, width = 2.25, height = 4.5, units = "in")


#Panel e
fig5e <- ggplot(is_pathway_analysis, aes(x=-log10(p.adj), y=pathway)) + 
  geom_bar(stat="identity") + theme_classic() + ylab("") + xlab("-log10(Adjusted P-Value)") + 
  ggtitle("ImmuneScore")  + geom_vline(xintercept = -log10(0.05), color = "red", linetype ="dashed")
ggsave("Figure 5e ImmuneScore Pathway Analysis Barplot.pdf", fig5e, width = 3, height = 3.65, units = "in")


#Panel f
#Flow-sorted data boxplots
fig5f1 <- ggplot(kilgour_nad, aes(x = expression,y=cell_type)) +geom_boxplot(outlier.shape = NA) + theme_classic() + ylab("") + xlab("NAD+") + 
  geom_jitter(color="black", size=0.2, alpha=0.9) 

#cAMP ovarian data scatterplot
fig5f2 <- ggplot(ov_nad, aes(x=NAD,y=ImmuneScore)) + geom_point(size=0.5) + geom_smooth(method = "lm", se = FALSE) + theme_classic() +
  xlab("NAD+")

fig5f <- fig5f1/fig5f2
ggsave("Figure 5f Flow sorted and cAMP NAD vs Immune Infiltration Plots.pdf", fig5f, width = 2.55, height = 3.85, units = "in")



#### Figure 6 ----

#Panel a
fig6a1 <- ggplot(bindea_concordance, aes(x=concordance,y=-log10(p.adj))) + geom_point(aes(color = color), size = 1) + 
  scale_color_manual(values= c("#D9D9D9","#FBB4AE")) + theme_classic() +
  geom_text_repel(data = subset(bindea_concordance,pair %in% label_pairs_bindea), aes(label = pair, size = 4)) + 
  geom_hline(yintercept = -log10(0.05)) + xlim(-0.35,0.35)
ggsave("Figure 6a1 Concordance Volcano Plot n > 7 and Bindea Signatures Only.pdf", fig6a1, width = 2.5, height = 2.75, units = "in")

fig6a2 <- ggplot(bindea_concordance,aes(x=concordance,y=col2)) + geom_point(shape=124,size=3,color=bindea_concordance$col3) + 
  theme_classic() + ylab("") + xlab("") + theme(axis.ticks = element_blank(),axis.text = element_blank()) + xlim(-0.35,0.35)
ggsave("Figure 6a2 Concordance Rug Plot n > 7 and Bindea Signatures Only.pdf", fig6a2, width = 2.5, height = 1, units = "in")


#Panel b
fig6b <- ggplot(bindea_concordance_dots, aes(x=immune_signature,y=-log10(p.adj))) + geom_point(aes(color = color), size = 0.55) + 
  scale_color_manual(values= c("#D9D9D9","#FBB4AE")) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("") + ylab("-Log10 Adjusted P-Value")
ggsave("Figure 6b Cell Signatures Concordance of Metabolites Dot Plot.pdf", fig6b, height = 3.5, width = 2.75)


#Panel c
fig6c <- ggplot(histamine_mc_exp, aes(x=Histamine, y=Mast.Cells, color = CancerType)) + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + 
  theme_classic() + xlab("Histamine") + ylab("Mast Cells") + geom_point(size=0.2) +
  scale_color_manual(values=mapping_color$Color[match(sort(unique(histamine_mc_exp$CancerType)),mapping_color$Name)])
ggsave("Figure 6c Mast Cells Histamine Correlation Line Plot in Tumors.pdf", fig6c, width = 2.8, height = 2.8, units = "in")

#Panel d
fig6d <- ggplot(mastcell_concordance[1:21,], aes(x = -log10(p.adj), y = metabolite)) + 
  geom_bar(stat="identity") + theme_classic() +
  xlab("-Log10 FDR-Corrected P-Value") + ylab("") + ggtitle("Mast Cell Concordance") + 
  geom_vline(xintercept = -log10(0.05), color = "red", linetype ="dashed")
ggsave("Figure 6d Mast Cells Concordance of Metabolites.pdf", fig6d, height = 2.5, width = 3)

#Panel e
fig6e <- ggplot(histamine_hdc_exp, aes(x=Histamine, y=HDC, color = CancerType)) + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + 
  theme_classic() + xlab("Histamine") + ylab("HDC") + geom_point(size=0.2) +
  scale_color_manual(values=mapping_color$Color[match(sort(unique(histamine_hdc_exp$CancerType)),mapping_color$Name)])
ggsave("Figure 6e HDC Histamine Correlation Line Plot in Tumors.pdf", fig6e, width = 2.7, height = 2.8, units = "in")

#Panel f
fig6f <- ggplot(adc_kyn_exp, aes(x=Kynurenine, y=aDC, color = CancerType)) + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + 
  theme_classic() + xlab("Kynurenine") + ylab("aDC") + geom_point(size=0.2) +
  scale_color_manual(values=mapping_color$Color[match(sort(unique(adc_kyn_exp$CancerType)),mapping_color$Name)])
ggsave("Figure 6f aDC Kynurenine Correlation Line Plot in Tumors.pdf", fig6f, width = 2.7, height = 2.8, units = "in")

#### Supplemental Figure 2 ----
supp2 <- ggplot(histamine_concordance, aes(x=pair, y=value)) + geom_boxplot() + 
  theme_classic() + xlab("") + ylab("Concordance") + geom_jitter(shape=16, position=position_jitter(0.2), aes(colour=variable)) +
  scale_color_manual(values=mapping_color$Color[match(sort(unique(histamine_concordance$variable)),mapping_color$Name)]) +
  theme(axis.text.x = element_text(angle = 345))
ggsave("Supp Figure 2 Mast Cells Histmamine Concordance Boxplot in Tumors.pdf", supp2, width = 2, height = 3, units = "in")