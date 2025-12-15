source('~/script/circadian/circadian_core.R')

#FIG. 1B UMAP of PBMC atlas
test<-readRDS("~/analysis/16individual.srt.annotated.rds")
test$type<-test$predicted.celltype.l2
celltypes<-test$type %>% unique() %>% sort()
cell_legend_data<-data.frame("type"=celltypes,"pos.x"=rep(1,length(celltypes)),"pos.y"=rev(1:length(celltypes)),"ID"=1:length(celltypes))
mappings<-cell_legend_data[,c("type","ID")]
colnames(mappings)[2]<-"new_label"
test<-addLabels(test,mappings)
dimplot_publication(test,reduction = "umap",reduction.name="umap",group.by = "new_label",label=T,
                    colors=generateColor(30,col.dist.min=0.3,seed = 10))+NoLegend()
ggsave("/lustre/home/acct-medll/medll/figures/circadian/Fig1B.pdf",width=4.5,height=5)
ggplot(cell_legend_data)+geom_point(aes(x=pos.x,y=pos.y,fill=type),shape=21,size=8,color="black")+
  geom_text(aes(x=pos.x,y=pos.y,label=ID),color="white")+
  geom_text(aes(x=pos.x+0.075,y=pos.y,label=type),hjust="left")+
  scale_x_continuous(limits=c(1,2))+scale_fill_manual(values=generateColor(30,col.dist.min=0.3,seed = 10))+
  theme_nothing()
ggsave("/lustre/home/acct-medll/medll/figures/circadian/Fig1B_legend.pdf",width=3,height=6)

#FIG. 1C composition of cells at each CT and markers
test<-subset(test,type!="Eryth"&type!="Platelet")
meta<-test@meta.data
meta$CT<-factor(meta$CT,levels=sort(unique(meta$CT)))
meta$type<-meta$predicted.celltype.l2
celltypes<-sort(unique(meta$type))
meta$type<-factor(meta$type,levels=rev(celltypes))


plotFIG1C1<-ggplot(meta)+geom_bar(aes(x=type,fill=CT),position="fill",color="black")+
  scale_fill_d3(alpha = 0.75)+scale_y_continuous(expand = c(0,0))+
  theme(plot.margin=margin(0.5,0,0,0,"cm"),
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.text.y=element_text(color="black"))+
  ylab("")+xlab("")+NoLegend()

plotFIG1C1_legend<-ggplot(meta)+geom_bar(aes(x=type,fill=CT),position="fill",color="black")+
  scale_fill_d3(alpha = 0.75)+scale_y_continuous(expand = c(0,0))+
  theme(plot.margin=margin(0.5,0,0,0,"cm"),
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.text.y=element_text(color="black"))+
  ylab("")+xlab("")
plotFIG1C1_legend<-get_legend(plotFIG1C1_legend)
grid.draw(plotFIG1C1_legend)
ggsave("/lustre/home/acct-medll/medll/figures/circadian/Fig1C_extraLegend.pdf",width=2,height=2)

plotFIG1C2<-DotPlot(test,rev(c("AXL","SIGLEC6",
                               "MS4A1","CD27",
                               "CD14","FCGR3A",
                               "CD3D","CD4","CD8A","FOXP3",
                               "CD1C","CLEC9A","HBA1",
                               "TRGC1","TRDC","CD34","SLC4A10",
                               "GNLY","GZMK","GZMH","LILRA4",
                               "JCHAIN","PPBP")),group.by="type",
                    cols = c("#507AAF","#BE5C37"))+scale_y_discrete(limits=sort(unique(test$type)))+
  scale_size_continuous(range = c(0,4))+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=9),
        plot.margin=margin(0,0,0,0,"cm"),axis.text.y=element_text(size=9),
        panel.border = element_rect(color = "black",linewidth = 1),
        legend.text = element_text(size=10),legend.title=element_text(size=10))+
  xlab("")+ylab("")+coord_flip()

ggarrange(plotFIG1C1,plotFIG1C2,nrow=2,align="v",heights = c(1,5))
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig1C.pdf',width = 6,height=6)

# FIG. 1D confusion matrix, manual annotation vs Azimuth
test<-readRDS('~/analysis/circadian/R/16individual.srt.annotated.rds')
## method1: by cell proportion
predicted_labels <- test$predicted.celltype.l2
manual_labels <- test$manual_NI
confusion_matrix <- table(predicted = predicted_labels, manual = manual_labels)
normalized_confusion <- prop.table(confusion_matrix, margin = 2)
# Convert confusion matrix to a data frame
confusion_df <- as.data.frame(confusion_matrix)
confusion_df<-confusion_df %>% group_by(predicted) %>% mutate(total=sum(Freq)) %>% mutate(proportion=Freq/total)
# Plot
ggplot(confusion_df, aes(x = manual, y = predicted, fill = proportion)) +
  geom_tile(color="white") +
  #geom_text(aes(label = Freq), color = "black") +
  scale_fill_gradientn(colours = c("#2166AC","#67A9CF","#D1E5F0","#F7F7F7",
                                  "#FDDBC7", "#EF8A62", "#B2182B")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black"),
        axis.text.y = element_text(color="black")) +
  labs(x = "manual annotations", y = "Azimuth annotations")

## method2: by correlation of transcription
mat_azi<-srt2bulkMatrix(test,split.by = "predicted.celltype.l2",layer = "count")
mat_manual<-srt2bulkMatrix(test,split.by = "manual_NI",layer = "count")
mat_manual<-srt2bulkMatrix(test,split.by = "manual.level2",layer = "count")
mat_azi<-mat_azi[,colnames(mat_azi)!="Platelet"&colnames(mat_azi)!="Eryth"]
mat_manual<-mat_manual[,colnames(mat_manual)!="Platelet"]
mat_azi<-mat_azi[rowSums(mat_azi)!=0,]
mat_manual<-mat_manual[rowSums(mat_manual)!=0,]
shared.features<-intersect(rownames(mat_azi),rownames(mat_manual))
mat_azi<-mat_azi[shared.features,]
mat_manual<-mat_manual[shared.features,]

mat_azi_z <- t(scale(t(mat_azi)))
mat_manual_z <- t(scale(t(mat_manual)))
# Calculate correlation matrix between cell types of the two dataframes
cor_matrix<-cor(mat_azi_z, mat_manual_z,method = "spearman")

library(reshape2)
cor_long <- melt(cor_matrix)
colnames(cor_long) <- c("Azimuth", "manual", "correlation")

ggplot(cor_long, aes(x = Azimuth, y = manual, fill = correlation)) +
  geom_tile(color="black") +
  #geom_text(aes(label = round(correlation, 2)), color = "white", size = 3) +
  scale_fill_gradientn(colours = c("#2166AC","#67A9CF","#D1E5F0","#F7F7F7",
                                            "#FDDBC7", "#EF8A62", "#B2182B"),limits = c(-1, 1)) +
                                              theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Expression Correlation", 
       y = "manual annotated cell types", 
       x = "Azimuth annotation")+
  #scale_y_discrete(limits=c("B_naive","B_memory","Plasma cell","CD14⁺Monocytes","CD16⁺Monocytes",
  #                          "mDC","pDC","CD4_naive_CCR7","CD4_TCM_AQP3","CD4_TEM_GNLY","CD4_TEM_ANXA1","CD4_TEM_GZMK","CD4_Treg_FOXP3",
  #                          "CD8_naive_LEF1","CD8_TEM_CMC1","CD8_TEM_GNLY","CD8_TEM_ZNF683",
  #                          "CD8_MAIT_SLC4A10","TNK_proliferatig_MKI67","NK_CD56ᵈⁱᵐ","NK_CD56ᵇʳⁱᵍʰᵗ","γδT","dnT_LYST",
  #                          "HSC_CD34","MS4A2_Granulocyte"))+
  scale_y_discrete(limits=c("cB01_TCL1A+B_naive","cB02_EGR1+B_intermediate","cB03_IgG+IgA+B_memory","cB04_Plasma",
                            "cM01_CD14+FCGR3A-Monocyte","cM02_CD14-FCGR3A+Monocyte","cM03_CLEC9A+cDC1","cM04_FCER1A+cDC2",
                            "cM05_LILRA4+pDC","cTNK01_CD4+CCR7+Tn","cTNK03_CD4+PTGER2+Tcm","cTNK12_CD4+cTRBV+Tcm-like","cTNK08_CD4+NKG7+CTL",
                            "cTNK09_CD4+FOXP3+Treg","cTNK02_CD8+CCR7+Tn","cTNK04_CD8+GZMK+Tem","cTNK15_CD8+GZMH+TEMRA/TEFF","cTNK13_CD8+cTRBV+Tcm-like",
                            "cTNK14_CD8+IKZF2+IEL","cTNK07_CD8+SLC4A10+MAIT","cTNK10_MKI67+proliferating","cTNK05_GNLY+NK","cTNK06_gdT","cTNK11_LYST+dnT",
                            "cS01_CD34+HSPC","cM05_MS4A2+Granulocyte"))+
  scale_x_discrete(limits=c("B naive","B intermediate","B memory","Plasmablast","CD14 Mono","CD16 Mono",
                            "cDC1","cDC2","pDC","ASDC","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","Treg",
                            "CD8 Naive","CD8 TEM","CD8 TCM","MAIT","CD4 Proliferating","CD8 Proliferating",
                            "NK Proliferating","NK","NK_CD56bright","gdT","dnT","HSPC","ILC"))


# FIG. 1E confusion matrix, CYTOF annotation vs Azimuth
## method 1, by expression correlation between protein and RNA
test<-readRDS("~/analysis/circadian/R/4multiplexedindividual.srt.annotated.rds")
mat<-srt2bulkMatrix(test,split.by = "predicted.celltype.l2",layer = "count")
mat<-mat[CYTOF2RNA_FEATURES %>% as.vector(),]
mat<-mat[rownames(mat)!="CD151",]
mat<-mat[rownames(mat)!="NA",]
mat<-mat[,colnames(mat)!="Platelet"]


test.cytof<-readRDS('~/analysis/circadian/R/CYTOF.300k.reannotation.srt.rds')
test.cytof@meta.data$type<-test.cytof@meta.data$manual.new
mat.cytof<-srt2bulkMatrix(subset(test.cytof,!is.na(type)),split.by = "type",layer = "count")
rownames(mat.cytof)<-CYTOF2RNA_FEATURES[rownames(mat.cytof)]
mat.cytof<-mat.cytof[!rownames(mat.cytof)%in%c("CD15","CD45RA"),]

#mat<-mat[,colnames(mat) %in% colnames(mat.cytof)]
mat_z <- t(scale(t(mat)))
mat.cytof_z <- t(scale(t(mat.cytof)))

# Calculate correlation matrix between cell types of the two dataframes
cor_matrix<-cor(mat_z, mat.cytof_z,method = "spearman")

# Convert to long format for ggplot
library(reshape2)
cor_long <- melt(cor_matrix)
colnames(cor_long) <- c("scRNAseq", "CYTOF", "correlation")

ggplot(cor_long, aes(x = scRNAseq, y = CYTOF, fill = correlation)) +
  geom_tile(color="black") +
  #geom_text(aes(label = round(correlation, 2)), color = "white", size = 3) +
  scale_fill_gradientn(colours = c("#2166AC","#67A9CF","#D1E5F0","#F7F7F7",
                  "#FDDBC7", "#EF8A62", "#B2182B"),limits = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Expression Correlation", 
       x = "single cell RNA seq cell types (Azimuth)", 
       y = "CYTOF cell types")+
  scale_x_discrete(limits=c("B naive","B intermediate","B memory","Plasmablast","CD14 Mono","CD16 Mono",
                            "cDC1","cDC2","pDC","ASDC","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","Treg",
                            "CD8 Naive","CD8 TEM","CD8 TCM","MAIT","CD4 Proliferating","CD8 Proliferating",
                            "NK Proliferating","NK","NK_CD56bright","gdT","dnT","HSPC","ILC"))+
  scale_y_discrete(limits=c("B naive","B memory","Plasmablast","CD14 Mono","CD16 Mono","cDC2","pDC",
                            "CD4 Naive","CD4 TCM","CD4 CTL","Treg","CD8 Naive","CD8 TEM","CD8 TCM","MAIT",
                            "NK","NK_CD56ᵇʳⁱᵍʰᵗ","gdT","dnT"))

## method 2, by correlation between cell type proportion in each sample
test<-readRDS("~/analysis/circadian/R/4multiplexedindividual.srt.annotated.rds")
test.cytof<-readRDS('~/analysis/circadian/R/CYTOF.300k.reannotation.srt.rds')

frequency_rna<-test@meta.data[,c("CT","individual","predicted.celltype.l2")] %>% group_by(CT,individual,predicted.celltype.l2) %>% summarise(cell.count=n())
frequency_CYTOF<-test.cytof@meta.data[,c("CT","individual","manual.new")] %>% group_by(CT,individual,manual.new) %>% summarise(cell.count=n())

frequency_rna<-frequency_rna %>% group_by(CT,individual) %>% mutate(total.cell=sum(cell.count)) %>% mutate(prop.cell=cell.count/total.cell)
frequency_CYTOF<-frequency_CYTOF %>% group_by(CT,individual) %>% mutate(total.cell=sum(cell.count)) %>% mutate(prop.cell=cell.count/total.cell)
#share.type<-intersect(frequency_rna$predicted.celltype.l2,frequency_CYTOF$manual.new)
#frequency_rna<-frequency_rna[frequency_rna$predicted.celltype.l2 %in% share.type,]
#frequency_CYTOF<-frequency_CYTOF[frequency_CYTOF$manual.new %in% share.type,]
frequency_CYTOF<-frequency_CYTOF[frequency_CYTOF$individual!="WMY",]

frequency_rna$sample<-paste0(frequency_rna$individual,"_",frequency_rna$CT)
frequency_CYTOF$sample<-paste0(frequency_CYTOF$individual,"_",frequency_CYTOF$CT)

frequency_rna<-spread(frequency_rna[c("sample","predicted.celltype.l2","prop.cell")],
                      key=sample,value=prop.cell) %>% column_to_rownames(.,"predicted.celltype.l2")
frequency_CYTOF<-spread(frequency_CYTOF[c("sample","manual.new","prop.cell")],
                      key=sample,value=prop.cell) %>% column_to_rownames(.,"manual.new")
frequency_rna[is.na(frequency_rna)]<-0
frequency_CYTOF[is.na(frequency_CYTOF)]<-0
frequency_rna<-as.data.frame(t(frequency_rna))
frequency_CYTOF<-as.data.frame(t(frequency_CYTOF))
cor_matrix<-cor(frequency_rna, frequency_CYTOF,method = "spearman")

# Convert to long format for ggplot
library(reshape2)
cor_long <- melt(cor_matrix)
colnames(cor_long) <- c("scRNAseq", "CYTOF", "correlation")

ggplot(cor_long, aes(x = scRNAseq, y = CYTOF, fill = correlation)) +
  geom_tile(color="black") +
  #geom_text(aes(label = round(correlation, 2)), color = "white", size = 3) +
  scale_fill_gradientn(colours = c("#2166AC","#67A9CF","#D1E5F0","#F7F7F7",
                                            "#FDDBC7", "#EF8A62", "#B2182B"),limits = c(-1, 1)) +
                                              theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Proportion Correlation", 
       x = "single cell RNA seq cell types (Azimuth)", 
       y = "CYTOF cell types")+
  scale_x_discrete(limits=c("B naive","B intermediate","B memory","Plasmablast","CD14 Mono","CD16 Mono",
                            "cDC1","cDC2","pDC","ASDC","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","Treg",
                            "CD8 Naive","CD8 TEM","CD8 TCM","MAIT","CD4 Proliferating","CD8 Proliferating",
                            "NK Proliferating","NK","NK_CD56bright","gdT","dnT","HSPC","ILC"))+
  scale_y_discrete(limits=c("B naive","B memory","Plasmablast","CD14 Mono","CD16 Mono","cDC2","pDC",
                            "CD4 Naive","CD4 TCM","CD4 CTL","Treg","CD8 Naive","CD8 TEM","CD8 TCM","MAIT",
                            "NK","NK_CD56ᵇʳⁱᵍʰᵗ","gdT","dnT"))

# FIG. S1XX, plot showing low expression level of genes 
count.mat<-LayerData(test,features=CIRCADIAN_GENES_MAIN,layer = "count") %>% as.data.frame()
count.mat<-rownames_to_column(count.mat,var = "features")
count.mat.long<-gather(count.mat,key="observations",value="counts",-features) 
meta.mat<-test@meta.data[,c("nCount_RNA","type")]
meta.mat<-rownames_to_column(meta.mat,var="observations")
count.mat.long<-left_join(count.mat.long,meta.mat,by="observations")
count.mat.long.summar<-count.mat.long %>% group_by(features,type) %>% 
  summarise(ratio.gt.1=mean(counts>1))
count.mat.long.summar$ratio.gt.10<-ifelse(count.mat.long.summar$ratio.gt.1>0.1,"yes","no")
ggplot(count.mat.long.summar) +
  geom_point(aes(x = type, y = features, size = ratio.gt.1,color=ratio.gt.10),
             shape=21,fill=pal_viridis()(10)[5],stroke=1) +
  scale_size_continuous(name = "Expression Ratio",range = c(0, 8),breaks = seq(0, 1, 0.3),limits = c(0,1)) +
  scale_color_manual(name = "Ratio > 10%",values=c("white","black")) +
  scale_y_discrete(limits=rev(c("BMAL1","CLOCK","PER1","PER2","PER3","CRY1","CRY2","DBP","NR1D1","NR1D2","TEF","HLF","CIART")))+
  labs(x = "", y = "",title = "Ratio of cells with UMI count > 1") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle=60,hjust=1),
    legend.position = "top",legend.direction = "horizontal",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    title=element_text(size = 11, face = "bold")
  )


# FIG. 2A
srt.metacell<-readRDS("/tmpdata/LyuLin/analysis/circadian/R/seacell.0.08.bytype.rds")
# part1
VlnPlot(srt.metacell,CIRCADIAN_GENES_MAIN,group.by = "type",stack = T,flip = T,cols = generateColor(15,seed = 2025))+NoLegend()+xlab("")+theme(strip.text=element_text(size=10,face="plain"))
# part2
plotdata<-srt.metacell$type %>% table() %>% as.data.frame()

ggplot(plotdata)+
  geom_text(data=plotdata[plotdata$Freq<2500,],aes(x=get("."),y=Freq,label=Freq),hjust=-0.25,size=4,angle=90)+
  geom_text(data=plotdata[plotdata$Freq>=2500,],aes(x=get("."),y=Freq/2,label=Freq),size=4,angle=90)+
  geom_bar(aes(x=get("."),y=Freq),stat="identity",fill="darkblue",alpha=0.5,color="black")+
  theme_minimal_hgrid(font_size = 12)+theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
  xlab("")+ylab("number of metacells")+scale_y_continuous(limits = c(0,10000),expand = c(0,0))+xlim(plotdata$.)

# FIG.2B
JTK.batch1<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/JTK.result.filtered.addp2bkg.addp2t.bytype.0.04.12Healthy.rds')
JTK.batch2<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/JTK.result.filtered.addp2bkg.addcor2batch.addp2t.bytype.0.08.rds')
JTK.ALL<-rbind(JTK.batch1,JTK.batch2[-10])
JTK.ALL<-JTK.ALL %>% group_by(CycID,celltype) %>% mutate(n.individual=n())
JTK.ALL.filtered<-dplyr::filter(JTK.ALL,n.individual>=2)
plotCountOscilattingByMetaCelltype(JTK.ALL)
plotCountOscilattingByMetaCelltype(JTK.ALL.filtered,flip = T)

JTK.raw1<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/JTK.result.filtered.bytype.0.04.12Healthy.rds')
JTK.raw2<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/JTK.result.filtered.bytype.0.08.rds')
JTK.raw.ALL<-rbind(JTK.raw1,JTK.raw2)
JTK.raw.ALL.core<-JTK.raw.ALL[JTK.raw.ALL$CycID %in% CIRCADIAN_GENES_MAIN,]
JTK.raw.ALL.core<-JTK.raw.ALL.core %>% group_by(CycID,celltype) %>% mutate(log.p.adj=mean(-log10(ADJ.P)),n.individual=n())
JTK.raw.ALL.core$significant<-ifelse(JTK.raw.ALL.core$log.p.adj>=2,"YES","NO")
ggplot(JTK.raw.ALL.core)+geom_point(aes(x=CycID,y=fct_rev(celltype),size=n.individual,fill=log.p.adj,color=significant),shape=21,width=3)+theme_classic()+
  scale_fill_gradient2(low = "lightblue",high = "darkred",name="-log10(p.adj)")+scale_color_manual(values=c("white","black"),name="p<=0.01")+
  theme(axis.text.x=element_text(angle=60,hjust=1,color="black"),legend.box.margin=margin(),
        axis.text.y=element_text(color="black"),legend.margin=margin(),
        legend.position="left")+xlab("")+ylab('')
plotCountOscilattingByMetaCelltype(JTK.ALL.filtered,flip = T,show.zero = F)

#FIG. 2C
plot.list<-list()
for (feature in c("BMAL1","CLOCK","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2")) {
  plot.list[[feature]]=plotMetaCellByIndividual(srt=srt.metacell,cell.type = "CD14 Mono",feature = feature,layer = "data",merge = T)
}
ggarrange(plotlist = plot.list)


# FIG. 2D
srt.metacell<-readRDS('~/analysis/circadian/R/seacell.16individual.rds')
srt.metacell.TF<-readRDS('~/analysis/circadian/R/seacell.16individual.TF.rds')
DefaultAssay(srt.metacell.TF)<-"TF"
srt.metacell.TF<-subset(srt.metacell.TF,!is.na(type))
VlnPlot(srt.metacell.TF,"CLOCK",group.by = "type",pt.size = 0,cols=generateColor(31,col.dist.min = 0.3,alpha = 0.5))+NoLegend()+xlab("")+
  ylab("transcription factor\nactivity")

# FIG. 2E
JTK_result_all<-readRDS('~/analysis/circadian/R/JTK.16individual.cor2batch0.2.fold2bkg1.5.fold2trough1.5.rds')
plotdata<-JTK_result_all[JTK_result_all$CycID %in% CIRCADIAN_GENES_MAIN,]
plotdata<-plotdata %>% group_by(CycID,individual) %>% summarise(phase=median(LAG),mean.log10.padj=mean(-log10(ADJ.P))) %>% as.data.frame()
ggplot(plotdata)+geom_vline(aes(xintercept=CycID))+
  geom_boxplot(aes(x=CycID,y=phase),fill="white",outlier.size = 0)+geom_jitter(aes(x=CycID,y=phase,fill=mean.log10.padj),color="black",shape=21,size=2)+
  ylim(c(0,23))+scale_x_discrete(limits = levels(plotdata$CycID)) +
  labs(x = "", y = "peaking time") + scale_fill_gradient2(low="lightblue",mid="blue",high="darkblue",midpoint = 12.5,name="-log10(p.adj)")+
  theme_minimal()+theme(axis.text.x = element_text(angle = 60, hjust = 1, color="black"),
                        axis.text.y = element_text(color="black"),
                        plot.title = element_text(hjust=0.5),legend.direction = "horizontal",legend.position = "top")

# FIG. S3A
plotdata<-srt.metacell.TF$type %>% table %>% as.data.frame()
ggplot(plotdata)+geom_bar(aes(x=get("."),y=Freq),stat="identity",fill="darkblue")+
  geom_text(data=plotdata[plotdata$Freq>10,],aes(x=get("."),y=Freq/8,label=Freq),angle=90,color="white")+
  geom_text(data=plotdata[plotdata$Freq<=10,],aes(x=get("."),y=Freq,label=Freq),angle=90,hjust=-0.5)+
  scale_y_log10()+theme_minimal()+
  theme(axis.text.x = element_text(angle=60,hjust=1))+xlab("")+ylab("number of metacells")

# FIG. S3B
plotdata<-srt.metacell.TF@meta.data[,c("type","nFeature_RNA")] %>% group_by(type) %>% summarise(median_nFeature=median(nFeature_RNA))
ggplot(plotdata)+geom_bar(aes(x=type,y=median_nFeature),stat="identity",fill="darkred")+
  geom_text(aes(x=type,y=median_nFeature/2,label=round(median_nFeature)),angle=90,color="white")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=60,hjust=1))+xlab("")+ylab("median number of\nexpressing genes")

# FIG. 3A
celltypes<-unique(srt.metacell$type)
JTK.result.filtered<-JTK_result_all[JTK_result_all$ADJ.P<=0.01,]
JTK.result.filtered<-JTK.result.filtered[,c("CycID","celltype")] %>% unique
share.data=NULL
for(type.this in celltypes){
  for(type.that in celltypes){
    n.gene.type.this=nrow(JTK.result.filtered[JTK.result.filtered$celltype==type.this,])
    n.share.gene=length(intersect(JTK.result.filtered[JTK.result.filtered$celltype==type.this,"CycID"],
                                  JTK.result.filtered[JTK.result.filtered$celltype==type.that,"CycID"]))
    if(n.gene.type.this==0){
      this.share.data=data.frame(from.type=type.this,to.type=type.that,shared.ratio=NA)
    }else{
      this.share.data=data.frame(from.type=type.this,to.type=type.that,shared.ratio=n.share.gene/n.gene.type.this)
    }
    if(is.null(share.data)){
      share.data=this.share.data
    }else{
      share.data=rbind(share.data,this.share.data)
    }
  }
}
invalid.type=share.data[is.na(share.data$shared.ratio),"from.type"]

ggplot(share.data[!is.na(share.data$shared.ratio)&!(share.data$from.type%in%invalid.type)&!(share.data$to.type%in%invalid.type),])+
  geom_raster(aes(x=from.type,y=to.type,fill=shared.ratio))+
  scale_fill_gradientn(colours = c("#283168","#009fc0","#f0f0cc","white"),values = c(0,0.5,0.9,1))+
  theme(axis.text.x = element_text(angle=60,hjust=1),plot.margin=margin(0,0,0,0,"cm"))+ylab("")+xlab("")

# FIG. 3B
result_list<-readRDS('~/analysis/circadian/R/enrichGO.16individual.rds')
plotdata=NULL
for(celltype in names(result_list)){
  this.data=result_list[[celltype]]@result['GO:0048511',]
  this.data$type=celltype
  if(is.null(plotdata)){
    plotdata=this.data
  }else{
    plotdata=rbind(plotdata,this.data)
  }
}
plotdata$GeneRatio.num<-plotdata$GeneRatio %>% sapply(., function(x) eval(parse(text = x)))
plotdata<-arrange(plotdata,by=desc(p.adjust))
plotdata$maintype<-CELL_TYPES[plotdata$type]
ggplot(plotdata)+geom_point(aes(x=fct_inorder(type),y=-log10(p.adjust),size=GeneRatio.num,color=maintype))+
  geom_hline(yintercept=-log10(0.05),linetype="dashed")+
  coord_flip()+theme_linedraw()+xlab("gene enrichment in\nrhythmic process (GO:0048511)")+
  scale_size_continuous(name="gene ratio")+theme(axis.title=element_text(color="black",size=12))+
  scale_color_manual(values = generateColor(n = length(unique(plotdata$maintype)),alpha = 0.75),name="main group")

# FIG. 3C
plotdata=NULL
for(celltype in names(result_list)){
  this.data=result_list[[celltype]]@result
  this.data=arrange(this.data,by=desc(p.adjust))
  this.data=this.data[this.data$p.adjust<0.05,]
  if(nrow(this.data)==0){
    next
  }
  this.data$type=celltype
  if(is.null(plotdata)){
    plotdata=this.data
  }else{
    plotdata=rbind(plotdata,this.data)
  }
}
IDs<-plotdata$ID %>% table()
IDs<-IDs[IDs==1]
plotdata<-plotdata[plotdata$ID %in% names(IDs),]
plotdata$GeneRatio.num<-plotdata$GeneRatio %>% sapply(., function(x) eval(parse(text = x)))
plotdata$maintype<-CELL_TYPES[plotdata$type]
plotdata<-plotdata[c("ID","type","p.adjust","qvalue","GeneRatio.num","Description")]
plotdata.label<-plotdata %>% group_by(type) %>% slice_max(.,order_by=-log10(p.adjust),n = 3)
ggplot(plotdata)+geom_point(aes(x=qvalue,y=-log10(p.adjust),color=type,size=GeneRatio.num))+
  geom_text_repel(data=plotdata.label,aes(x=qvalue,y=-log10(p.adjust),label=ID),size=3)+
  theme_linedraw()+theme(axis.text.x = element_text(angle=60,hjust=1),legend.position="top")+
  scale_color_manual(values = generateColor(n = length(unique(plotdata$type)),alpha = 0.75),name="cell type")+
  scale_size_continuous(name="gene ratio")+
  guides(color="none")+facet_wrap(~type)


# FIG. 3D
JTK_result_all<-readRDS('~/analysis/circadian/R/JTK.16individual.cor2batch0.2.fold2bkg1.5.fold2trough1.5.rds')
JTK_result_all<-JTK_result_all %>% group_by(CycID,celltype) %>% mutate(n.occurance=n())
JTK_result_filtered<-JTK_result_all[JTK_result_all$n.occurance>=2&JTK_result_all$CycID!="HBB",]
JTK_result_filtered<-JTK_result_filtered[!grepl("ENSG",JTK_result_filtered$CycID),]
plot1<-plotOscillatingGeneSummarise(JTK_result_filtered,c("CD14 Mono"),x.axis="phase.consist",
      y.axis="fold.peak2trough",method = "JTK",xlim=c(0,0.5),ylim=c(1.5,15),
      int.color.index = 2)+xlab("")
plot2<-plotOscillatingGeneSummarise(JTK_result_filtered,c("CD4 Naive"),x.axis="phase.consist",
                                    y.axis="fold.peak2trough",method = "JTK",xlim=c(0,0.5),ylim=c(1.5,5),
                                    int.color.index = 4)+xlab("")+ylab("")
plot3<-plotOscillatingGeneSummarise(JTK_result_filtered,c("CD4 TCM"),x.axis="phase.consist",
                                    y.axis="fold.peak2trough",method = "JTK",xlim=c(0,0.5),ylim=c(1.5,5.5),
                                    int.color.index = 5)+xlab("")+ylab("")
plot4<-plotOscillatingGeneSummarise(JTK_result_filtered,c("CD8 Naive"),x.axis="phase.consist",
                                    y.axis="fold.peak2trough",method = "JTK",xlim=c(0,0.5),ylim=c(1.5,5.5),
                                    int.color.index = 7)+xlab("")
plot5<-plotOscillatingGeneSummarise(JTK_result_filtered,c("CD8 TEM"),x.axis="phase.consist",
                             y.axis="fold.peak2trough",method = "JTK",xlim=c(0,0.5),ylim=c(1.5,6),
                             int.color.index = 9)+ylab("")
plot6<-plotOscillatingGeneSummarise(JTK_result_filtered,c("NK"),x.axis="phase.consist",
                             y.axis="fold.peak2trough",method = "JTK",xlim=c(0,0.5),ylim=c(1.5,6),
                             int.color.index = 12)+xlab("")+ylab("")
ggarrange(plotlist = list(plot1,plot2,plot3,plot4,plot5,plot6),nrow = 2,ncol=3,common.legend = T)

# FIG. 3E
JTK_result_filtered<-JTK_result_filtered %>% group_by(CycID) %>% mutate(n.occurance=n())
top100<-JTK_result_filtered[c("CycID","n.occurance")] %>% unique() %>% arrange(by=desc(n.occurance)) %>% head(100)
plotdata<-JTK_result_filtered[JTK_result_filtered$CycID %in% top100$CycID,]
plotdata<-plotdata %>% group_by(CycID) %>% mutate(median.phase=median(LAG),mean.log10.p.adj=mean(-log10(ADJ.P)))
plotdata2<-plotdata[c("CycID","n.occurance","mean.log10.p.adj")] %>% unique() %>% arrange(desc(n.occurance))
plot1<-ggplot(plotdata)+geom_bar(aes(x=fct_infreq(CycID),fill=median.phase),color="black")+
  geom_text(data=plotdata2,aes(x=CycID,y=n.occurance,label=CycID),angle=60,vjust=0.5,hjust=-0.11,size=3)+
  scale_fill_gradientn(name="peaking time",colours=c("black","darkgrey","lightgrey","white","lightgrey","darkgrey","black"))+
  theme_half_open()+scale_y_continuous(expand = c(0,0),limits = c(0,115))+scale_x_discrete(expand=c(0,2.5))+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.margin = margin(0,0,0,0,"cm"),
        legend.direction="horizontal",
        legend.title.position = "top")+ylab("total occurence")+xlab("")
legend1<-get_legend(plot1)
plot2<-ggplot(plotdata2)+geom_tile(aes(x=fct_inorder(CycID),y="significance",fill=mean.log10.p.adj),color="black",size=0.5)+
  scale_x_discrete(expand=c(0,2.5))+theme_half_open()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line = element_blank(),
        axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        panel.background = element_rect(fill="white"),plot.margin = margin(0,0,0,0),
        legend.direction="horizontal",legend.title.position = "top")+
  scale_fill_gradient2(low="darkblue",high = "darkred",name="mean -log10(p.adj)")+
  xlab("")+ylab("")
legend2<-get_legend(plot2)
legend.plot<-ggarrange(as_ggplot(legend1),as_ggplot(legend2),ncol=2,nrow=1)
plot1<-plot1+NoLegend()
plot2<-plot2+NoLegend()
ggarrange(plot1,plot2,nrow=2,ncol=1,align="v",heights = c(6,1))+
  inset_element(legend.plot, left = 0.7, bottom = 0.7, right = 0.95, top = 0.85)

# FIG. S3 XX
JTK_result_all<-readRDS('~/analysis/circadian/R/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
JTK_result_all<-JTK_result_all %>% group_by(CycID,celltype) %>% mutate(n.occurance=n())
JTK_result_filtered<-JTK_result_all[JTK_result_all$fold_peak2trough>=1.5&JTK_result_all$median_expression>10&JTK_result_all$fold_to_background>2,]
used.celltype<-"cDC2"
JTK_result_used<-JTK_result_filtered[!is.na(JTK_result_filtered$CycID)&JTK_result_filtered$celltype==used.celltype&JTK_result_filtered$PER>=20,]
JTK_result_used<-JTK_result_used %>% group_by(CycID,celltype) %>% mutate(n.individual=n())
JTK_result_used<-JTK_result_used %>% group_by(individual) %>% mutate(percentile_rank=percent_rank(ADJ.P))
JTK_result_used<-JTK_result_used %>% group_by(CycID) %>% mutate(median.percentile=median(percentile_rank),median.Fpeak2trough=median(fold_peak2trough))
top100.genes<-(JTK_result_used[c("CycID","n.individual")] %>% unique() %>% arrange(desc(n.individual)) %>% head(.,100))$"CycID"
ggplot(JTK_result_used[JTK_result_used$CycID%in%top100.genes,])+geom_bar(aes(x=fct_infreq(CycID),fill=median.percentile),color="black")+
  geom_text(aes(x=CycID,y=n.individual,label=CycID),angle=60,vjust=0.5,hjust=-0.11,size=3)+
  scale_fill_gradientn(name="median percentile of\nadjusted p",colours=rev(c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"))) +
  theme_half_open()+scale_y_continuous(expand = c(0,0),limits = c(0,4))+scale_x_discrete(expand=c(0,2.5))+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.margin = margin(1,0,0,0,"cm"),
        legend.direction="horizontal",
        legend.title.position = "top")+ylab("number of individual")+xlab("")+
  ggtitle(paste0("rank and prevalence of significant genes for ",used.celltype))

# FIG. 3F
data<-read.delim('~/analysis/pidata/tauFisher_test/GSE113883_hs_full_adj_meta2d/JTKresult_GSE113883_hs_full_adj.txt')
data<-data[data$ADJ.P<0.05,]
plot1<-ggplot(data)+geom_point(aes(x=log10(AMP),y=-log10(ADJ.P),fill=LAG),shape=21,size=3)+
  geom_text_repel(aes(x=log10(AMP),y=-log10(ADJ.P),label=CycID),hjust=0,vjust=-1)+
  scale_fill_gradientn(colours = c("white","lightgrey","darkgrey","black"),name="relative\npeaking time")+theme_half_open()+
  ylab("-log10(p.adj)")+xlab("log10(amplitude)")+
  theme(legend.position = c(0.75,0.9),legend.direction = "horizontal",
        legend.title.position = "top",plot.margin = margin(0,0,0,0))

plot2<-ggplot(data)+geom_point(aes(x=1,y=-log10(ADJ.P),fill=LAG),shape=21,size=3)+
  geom_text(data=data[data$CycID %in% c("DDIT4","FKBP5","TSC22D3","CXCR4"),],aes(x=1.2,y=-log10(ADJ.P),label=CycID),hjust=0)+
  scale_fill_gradientn(colours = c("white","lightgrey","darkgrey","black"))+theme_nothing()+
  ylab("-log10(p.adj)")+NoLegend()+theme(plot.margin = margin(0,1,0,0,"cm"))+scale_x_continuous(expand=c(0,0),limits=c(0.8,2))

ggarrange(plot1,plot2,ncol = 2,nrow=1,align="hv",widths = c(3,1))

# FIG. 3G
plot.list<-list()
for (celltype in c("CD14 Mono","CD4 Naive", "CD4 TCM", "CD8 Naive", "CD8 TEM", "NK")) {
  plot.list[[celltype]]=plotMetaCellByIndividual(srt.metacell,celltype,"DDIT4","data",merge = T,phase.adjust = JTK_result_all)
}
ggarrange(plotlist = plot.list,nrow=2,ncol=3)


#ggplot(plotdata)+geom_boxplot(aes(x=time,y=relative_expression))+
#  geom_point(aes(x=time,y=relative_expression),alpha=0.5,color="blue",size=3)+
#  scale_x_discrete(limits=CT_TIME_ORDER)+theme_bw()+
#  facet_wrap(~feature,scales="free_y")+ylab("z-score")+xlab("")+
#  ggtitle("CD14 Mono\np.adj: 0.016, peak: CT18")
#ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig2A1.pdf',width = 3,height=3)

#ggplot(plotdata2)+geom_boxplot(aes(x=time,y=relative_expression))+
#  geom_point(aes(x=time,y=relative_expression),alpha=0.5,color="blue",size=3)+
#  scale_x_discrete(limits=CT_TIME_ORDER)+theme_bw()+
#  facet_wrap(~feature,scales="free_y")+ylab("z-score")+xlab("")+
#  ggtitle("NK\np.adj: 0.001, peak: CT10")
#ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig2A2.pdf',width = 3,height=3)

#ggplot(plotdata3)+geom_boxplot(aes(x=time,y=relative_expression))+
#  geom_point(aes(x=time,y=relative_expression),alpha=0.5,color="blue",size=3)+
#  scale_x_discrete(limits=CT_TIME_ORDER)+theme_bw()+
#  facet_wrap(~feature,scales="free_y")+ylab("z-score")+xlab("")+
#  ggtitle("Plasma\np.adj: 0.037, peak: CT6")
#ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig2A3.pdf',width = 3,height=3)

#ggplot(plotdata4)+geom_boxplot(aes(x=time,y=relative_expression))+
#  geom_point(aes(x=time,y=relative_expression),alpha=0.5,color="blue",size=3)+
#  scale_x_discrete(limits=CT_TIME_ORDER)+theme_bw()+
#  facet_wrap(~feature,scales="free_y")+ylab("z-score")+xlab("")+
#  ggtitle("CD16 Mono\np.adj: 0.032, peak: CT4")
#ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig2A4.pdf',width = 3,height=3)

#plotdataFIG2B1<-test@meta.data[,c("nFeature_RNA","type")] %>% group_by(type) %>% summarise(median_genes=max(nFeature_RNA))
#plotdataFIG2B1<-left_join(plotdataFIG2B1,(JTK.individual %>% dplyr::filter(.,ADJ.P<=0.01))$celltype %>% table() %>% as.data.frame() %>% `colnames<-`(.,c("type","rhythmic")))
#plotdataFIG2B1$percent_rhythmic<-plotdataFIG2B1$rhythmic*100/plotdataFIG2B1$median_genes
#plotdataFIG2B1$percent_rhythmic<-round(plotdataFIG2B1$percent_rhythmic,2)
#plotdataFIG2B1[is.na(plotdataFIG2B1)]<-0

#FIG2B1<-ggplot(plotdataFIG2B1)+geom_bar(aes(y=type,x=percent_rhythmic),stat="identity",fill="black",alpha=0.5,color="black")+
#  geom_text(data=plotdataFIG2B1[plotdataFIG2B1$percent_rhythmic>=2,],aes(y=type,x=percent_rhythmic/2,label=percent_rhythmic),color="white")+
#  geom_text(data=plotdataFIG2B1[plotdataFIG2B1$percent_rhythmic<2,],aes(y=type,x=percent_rhythmic+1,label=percent_rhythmic),color="black")+
#  scale_y_discrete(limits=rev(celltypes %>% sort()))+scale_x_continuous(expand=c(0,0))+
#  theme(panel.grid = element_blank(),panel.background=element_rect(fill="white"))+
#  ylab("")+xlab("percentage of\nrhythmic genes(%)")

#plotdataFIG2B<-JTK.individual[JTK.individual$CycID %in% CIRCADIAN_GENES_MAIN,]
#plotdataFIG2B$celltype<-factor(plotdataFIG2B$celltype,levels=rev(celltypes %>% sort()))
#FIG2B2<-ggplot(plotdataFIG2B)+geom_point(aes(x=celltype,y=CycID,size=-log10(ADJ.P),color=log2(AMP)))+
#  scale_color_gradientn(colours = c("#507AAF","#BE5C37"))+
#  theme(axis.text.x = element_text(angle=45,hjust=1),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
#        panel.background = element_rect(fill="white"))+coord_flip()+
#  xlab("")+ylab("rhythmicity of\ncore circadian genes")

#ggarrange(FIG2B1,FIG2B2,ncol=2,widths=c(1.5,2.5),align="h")
#ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig2B.pdf',width=6.5,height=4.5)

#FIG. 3A
#plotdataFIG3A1<-test$predicted.celltype.l1.5 %>% table() %>% as.data.frame()
#test$type<-test$predicted.celltype.l1.5
#celltypes<-sort(unique(test$predicted.celltype.l1.5))
#colnames(plotdataFIG3A1)<-c("cell_type","number_of_cell")
#FIG3A1<-ggplot(plotdataFIG3A1)+geom_bar(aes(y=cell_type,x=log10(number_of_cell)),stat="identity",fill="blue",alpha=0.5,color="black")+
#  geom_text(aes(y=cell_type,x=log10(number_of_cell)/2,label=number_of_cell),color="white")+
#  scale_y_discrete(limits=rev(plotdataFIG3A1$cell_type))+scale_x_continuous(expand=c(0,0))+
#  theme(panel.grid = element_blank(),panel.background=element_rect(fill="white"))+ylab("")+xlab("log10(# of cell)")

#JTK_CYCL<-NULL
#for(block_data in list.files("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5/")){
#  block_path=paste0("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5/",block_data,"/")
#  block=getJTK_CYCLEouts(block_path,celltypes=celltypes,filter.data = F)
#  if(is.null(JTK_CYCL)){
#    JTK_CYCL=block
#  }else{
#    JTK_CYCL=rbind(JTK_CYCL,block)
#  }
#}

#JTK_CYCL_bulk<-read.delim('/lustre/home/acct-medll/medll/data/analysis/tauFisher_test/GSE113883_hs_full_adj_meta2d/meta2d_GSE113883_hs_full_adj.txt')
#JTK_CYCL_bulk<-dplyr::filter(JTK_CYCL_bulk,JTK_pvalue<0.01)
#plotdataFIG3A2<-test@meta.data[,c("nFeature_RNA","type")] %>% group_by(type) %>% summarise(median_genes=max(nFeature_RNA))
#plotdataFIG3A2<-left_join(plotdataFIG3A2,(JTK_CYCL %>% dplyr::filter(.,ADJ.P<=0.01,CycID %in% JTK_CYCL_bulk$CycID))$celltype %>% table() %>% as.data.frame() %>% `colnames<-`(.,c("type","rhythmic")))
#plotdataFIG3A2$percent_rhythmic<-plotdataFIG3A2$rhythmic*100/plotdataFIG3A2$median_genes
#plotdataFIG3A2$percent_rhythmic<-round(plotdataFIG3A2$percent_rhythmic,2)
#plotdataFIG3A2[is.na(plotdataFIG3A2)]<-0

#FIG3A2<-ggplot(plotdataFIG3A2)+geom_bar(aes(y=type,x=percent_rhythmic),stat="identity",fill="red",alpha=0.5,color="black")+
#  geom_text(data=plotdataFIG3A2[plotdataFIG3A2$percent_rhythmic>=2,],aes(y=type,x=percent_rhythmic/2,label=percent_rhythmic),color="white")+
#  geom_text(data=plotdataFIG3A2[plotdataFIG3A2$percent_rhythmic<2,],aes(y=type,x=percent_rhythmic+1,label=percent_rhythmic),color="black")+
#  scale_y_discrete(limits=rev(plotdataFIG3A1$cell_type))+scale_x_continuous(expand=c(0,0))+
#  theme(panel.grid = element_blank(),panel.background=element_rect(fill="white"),
#        axis.text.y=element_blank(),axis.ticks.y=element_blank())+
#  ylab("")+xlab("percentage of\nrhythmic genes(%)")

#plotdataFIG3A3<-JTK_CYCL[JTK_CYCL$CycID %in% CIRCADIAN_GENES_MAIN,]
#plotdataFIG3A3$celltype<-factor(plotdataFIG3A3$celltype,levels=rev(celltypes))
#FIG3A3<-ggplot(plotdataFIG3A3)+geom_point(aes(x=celltype,y=CycID,size=-log10(ADJ.P),color=log2(AMP)))+
#  scale_color_gradientn(colours = c("#507AAF","#BE5C37"))+
#  theme(axis.text.x = element_text(angle=45,hjust=1),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
#        panel.background = element_rect(fill="white"))+coord_flip()+
#  xlab("")+ylab("rhythmicity of\ncore circadian genes")

#ggarrange(FIG3A1,FIG3A2,FIG3A3,ncol=3,widths=c(1.5,1.1,2.7),align="h")
#ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig3A.pdf',width=8,height=4.5)

#FIG. 3B
#plotdataFIG3B<-NULL
#for(celltype in c("CD14 Mono","CD16 Mono","CD4 T")){
#  plotdataFIG3B.part1=plotPseudobulkByCT(test,celltype,features=c("ARNTL","CRY1","NR1D1","PER1","TEF"),return.data = T,rep=2,proportion=0.5)
#  plotdataFIG3B.part1$type=celltype
#  plotdataFIG3B.part2=plotPseudobulkByCT(test,celltype,features=c("ARNTL","CRY1","NR1D1","PER1","TEF"),return.data = T,rep=2,proportion=0.5)
#  plotdataFIG3B.part2$type=celltype
#  plotdataFIG3B.part2$time=case_when(plotdataFIG3B.part2$time=="CT1"~"CT25",
#                                     plotdataFIG3B.part2$time=="CT5"~"CT29",
#                                     plotdataFIG3B.part2$time=="CT9"~"CT33",
#                                     plotdataFIG3B.part2$time=="CT13"~"CT37",
#                                     plotdataFIG3B.part2$time=="CT17"~"CT41",
#                                     plotdataFIG3B.part2$time=="CT21"~"CT45")
#  plotdataFIG3B.part=rbind(plotdataFIG3B.part1,plotdataFIG3B.part2)
#  if(is.null(plotdataFIG3B)){
#    plotdataFIG3B=plotdataFIG3B.part
#  }else{
#    plotdataFIG3B=rbind(plotdataFIG3B,plotdataFIG3B.part)
#  }
#}
#summ.data<-plotdataFIG3B %>% group_by(feature,type) %>% summarise(min=min(count)+1)
#plotdataFIG3B[plotdataFIG3B$count!=0,]
#plotdataFIG3B<-left_join(plotdataFIG3B,summ.data,by=c("feature","type"))
#plotdataFIG3B$relative_expression<-plotdataFIG3B$count/plotdataFIG3B$min
#plotdataFIG3B.sup<-plotdataFIG3B %>% group_by(feature,type,time) %>% summarise(median_expression=median(relative_expression))

#use_color<-generateColor(30,col.dist.min=0.3,seed = 10)
#ggplot(plotdataFIG3B)+geom_boxplot(aes(x=time,y=relative_expression,color=type))+
#  geom_point(data=plotdataFIG3B.sup,aes(x=time,y=median_expression,color=type))+
#  geom_line(data=plotdataFIG3B.sup,aes(x=time,y=median_expression,group=type,color=type))+
#  scale_color_manual(values=c(use_color[4],use_color[5],use_color[6],use_color[9]))+
#  scale_x_discrete(limits=c(CT_TIME_ORDER,"CT25","CT29","CT33","CT37","CT41","CT45"))+
#  facet_wrap(~feature,scales = "free_y",nrow=2)+xlab("")+
#  theme(axis.text.x=element_text(angle=60,hjust=1,color="black",size=9),
#        axis.text.y=element_text(color="black",size=9),legend.text=element_text(color="black",size=9),
#        panel.background=element_rect(fill="white"),
#        legend.position=c(0.85,0.1))
#ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig3B.pdf',width=6,height=3)


#FIG. 3C
#JTK.result.filtered<-JTK_CYCL[JTK_CYCL$ADJ.P<=0.01&JTK_CYCL$CycID %in% JTK_CYCL_bulk$CycID,]
#share.data=NULL
#for(type.this in celltypes){
#  for(type.that in celltypes){
#    n.gene.type.this=nrow(JTK.result.filtered[JTK.result.filtered$celltype==type.this,])
#    n.share.gene=length(intersect(JTK.result.filtered[JTK.result.filtered$celltype==type.this,"CycID"],
#                                  JTK.result.filtered[JTK.result.filtered$celltype==type.that,"CycID"]))
#    if(n.gene.type.this==0){
#      this.share.data=data.frame(from.type=type.this,to.type=type.that,shared.ratio=NA)
#    }else{
#      this.share.data=data.frame(from.type=type.this,to.type=type.that,shared.ratio=n.share.gene/n.gene.type.this)
#    }
#    if(is.null(share.data)){
#      share.data=this.share.data
#    }else{
#      share.data=rbind(share.data,this.share.data)
#    }
#  }
#}
#ggplot(share.data[share.data$from.type!="Eryth"&share.data$to.type!="Eryth",])+geom_raster(aes(x=from.type,y=to.type,fill=shared.ratio))+
#  scale_fill_gradientn(colours = c("#283168","#009fc0","#f0f0cc","white"),values = c(0,0.5,0.9,1))+
#  theme(axis.text.x = element_text(angle=60,hjust=1),plot.margin=margin(0,0,0,0,"cm"))+ylab("")+xlab("")
#ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig3C.pdf',width=5,height=4)

#FIG. 3D
#unique_in_sc<-JTK_CYCL[JTK_CYCL$ADJ.P<=0.01&!(JTK_CYCL$CycID %in% JTK_CYCL_bulk$CycID)&JTK_CYCL$PER==24,"CycID"] %>% unique()
#shared_between_sc_bulk<-JTK_CYCL[JTK_CYCL$ADJ.P<=0.01&(JTK_CYCL$CycID %in% JTK_CYCL_bulk$CycID)&JTK_CYCL$PER==24,"CycID"] %>% unique()
#(LayerData(test,layer="count")[unique_in_sc,]>0) %>% rowSums() %>% median()
#(LayerData(test,layer="count")[shared_between_sc_bulk,]>0) %>% rowSums() %>% median()
#NK_genes_by_lag<-(JTK.result.filtered[JTK.result.filtered$celltype=="NK",] %>% arrange(.,by=LAG))$CycID
#plotGeneExpressionHeatmapByCT(test,"NK",NK_genes_by_lag)
#ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig3D.pdf',width=4,height=2.5)

#FIG. 3E
#gdT_genes_by_lag<-(JTK.result.filtered[JTK.result.filtered$celltype=="gdT",] %>% arrange(.,by=LAG))$CycID
#plotGeneExpressionHeatmapByCT(test,"gdT",gdT_genes_by_lag)
#ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig3E.pdf',width=4,height=2.5)

#FIG. 3F
#use_color<-generateColor(30,col.dist.min=0.3,seed = 10)[-11]
#ggplot(JTK.result.filtered[JTK.result.filtered$celltype!="Eryth",])+
#  geom_bar(aes(x=LAG,fill=celltype),alpha=0.5,color="black")+
#  facet_wrap(~celltype,scales = "free_y",ncol=6)+
#  scale_fill_manual(values=use_color)+NoLegend()+
#  theme(axis.text=element_text(color="black"))+ylab("count of peaking genes")+
#  xlab("ciradian time")
#ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig3F.pdf',width=8,height=4)


#FIG. 3G
#plotPhaseofCellType(JTK.result.filtered[JTK.result.filtered$celltype!="Eryth",],p.cutoff = 0.01,
#                    colors=generateColor(30,col.dist.min=0.3,seed = 10)[-11])+NoLegend()
#ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig3G.pdf',width=4,height=4)

#FIG. 4A
celltypes<-names(CELL_TYPES)
celltypes<-celltypes[!celltypes%in%c("Eryth","CD4 Proliferating","CD8 Proliferating","NK Proliferating",
                                     "NK_CD56bright","ASDC","cDC1","pDC","NK_CD56ᵇʳⁱᵍʰᵗ","Platelet","HSPC","ILC")]
adh<-c("SELPLG","SELL","CD44","ITGAL","ITGAM","ITGAX","ITGB2","ITGB1","ITGA2","ITGA4","ITGA5","ITGA6")
adhAchemo<-c("SELPLG","SELL","CD44","ITGAL","ITGAM","ITGAX","ITGB2","ITGB1","ITGA2","ITGA4","ITGA5","ITGA6",
             "CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CCR8","CCR9","CCR10","CXCR2","CXCR3","CXCR4","CXCR5","CXCR6","CX3CR1")
JTK_result_all<-readRDS('~/analysis/circadian/R/JTK.16individual.cor2batch0.2.fold2bkg1.5.fold2trough1.5.rds')
JTK_result_raw.part1<-readRDS('~/analysis/circadian/R/JTK.12individual.raw.rds')
JTK_result_raw.part2<-readRDS('~/analysis/circadian/R/JTK.4individual.raw.rds')
JTK_result_raw<-rbind(JTK_result_raw.part1,JTK_result_raw.part2)
adhesionAndChemokine<-JTK_result_raw[JTK_result_raw$CycID %in% adhAchemo,]
adhesionAndChemokine2<-JTK_result_all[JTK_result_all$CycID %in% adhAchemo,]
significance<-adhesionAndChemokine %>% group_by(CycID,celltype) %>% summarise(mean.log.adj.p=mean(-log10(ADJ.P)))
fold2background<-adhesionAndChemokine2 %>% group_by(CycID,celltype) %>% summarise(fold_to_background=median(fold_to_background))
significance$significant<-ifelse(significance$mean.log.adj.p<(-log10(0.05)),"not oscillating","oscillating")
plotdata<-left_join(significance,fold2background,by=c("CycID","celltype"))
plotdata$significant<-ifelse(is.na(plotdata$fold_to_background)&plotdata$significant=="not oscillating","no/low expression",plotdata$significant)

#adhesionAndChemokine<-JTK.individual[JTK.individual$CycID %in% adhAchemo&JTK.individual$celltype!="Eryth" ,]
#adhesionAndChemokine$`p.adj<0.05`<-case_when(adhesionAndChemokine$ADJ.P<0.05 ~ "oscillating",
#                                             adhesionAndChemokine$fold_to_background<2 ~ "no/low expression",
#                                             TRUE~"not oscillating")
for(celltype in celltypes){
  this.data=data.frame("CycID"=c("ITGA2","CCR3","CCR8","CCR9"),"celltype"=rep(celltype,4),
                       mean.log.adj.p=rep(NA,4),significant=rep("no/low expression",4),
                       fold_to_background=rep(NA,4))
  plotdata=rbind(plotdata,this.data)
}


bottom<-ggplot(plotdata)+geom_tile(aes(x=CycID,y=celltype,fill=significant),color="black",linewidth=0.5)+
  scale_y_discrete(limits=rev(celltypes))+scale_fill_manual(values=c("white","grey","red"))+
  scale_x_discrete(limits=adhAchemo)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,color="black",size=12),
        axis.text.y=element_text(color="black",size=10),legend.text=element_text(size=12),
        panel.border = element_rect(fill=NA,color="black"),legend.title = element_text(size=14),
        plot.margin = margin(0,0.5,0,0,"cm"))+
  xlab("")+ylab("")+guides(fill=guide_legend(title=NULL,position = "bottom"))


adhesionAndChemokine.f<-plotdata[plotdata$mean.log.adj.p>(-log10(0.05)),]
#mid<-ggplot(adhesionAndChemokine.f)+geom_boxplot(aes(x=CycID,y=LAG),outliers = FALSE)+
#  geom_jitter(aes(x=CycID,y=LAG,color=celltype))+
#  scale_x_discrete(limits=adhAchemo)+xlab("")+ylab("peaking time by\ncell type (CT)")+
#  scale_color_manual(values=generateColor(30,col.dist.min=0.3,seed = 10)[-11])+
#  theme(axis.text.x=element_blank(),axis.ticks = element_blank(),plot.margin = margin(0,0.5,0,0,"cm"),
#        axis.text.y=element_text(color="black",size=10),panel.grid.major.y=element_line(colour = "grey"),
#        panel.background = element_rect(fill=NA))+NoLegend()

adhesionAndChemokine.f$gene_type<-ifelse(adhesionAndChemokine.f$CycID %in% adh,"adhesion molecules","chemokine receptors")
top<-ggplot(adhesionAndChemokine.f)+geom_bar(aes(x=CycID,fill=gene_type),color="black")+
  scale_x_discrete(limits=adhAchemo)+xlab("")+ylab("\n\n\ncount of\ncell types")+scale_y_continuous(n.breaks = 4)+
  theme(axis.text.x=element_blank(),axis.ticks = element_blank(),plot.margin = margin(0,0,0,0,"cm"),
        axis.text.y=element_text(color="black",size=10),panel.grid.major.y=element_line(colour = "grey"),
        axis.title.y = element_text(margin = margin(r = -20)),legend.text=element_text(size=12),
        legend.title = element_text(size=14),
        panel.background = element_rect(fill=NA),legend.position="top")+
        guides(fill=guide_legend(title=NULL))

ggarrange(top,bottom,nrow=2,ncol=1,align = "v",heights = c(1,3))

# Fig. 4XX CYTOF result
srt.cytof<-readRDS('~/analysis/circadian/R/CYTOF.300k.reannotation.srt.rds')
JTK.cytof<-readRDS('~/analysis/circadian/R/CYTOF.JTK.newresult.filtered.addmedianexp.addp2t.bytype.rds')
JTK.cytof$CycID_RNA<-CYTOF2RNA_FEATURES[JTK.cytof$CycID]
JTK.cytof<-JTK.cytof %>% group_by(celltype,CycID_RNA) %>% summarise(ADJ.P=median(ADJ.P),median_expression=median(median_expression))
JTK.cytof$significant<-case_when(JTK.cytof$ADJ.P<0.05&JTK.cytof$median_expression>0.01~"oscillating",
                                 JTK.cytof$ADJ.P>0.05&JTK.cytof$median_expression>0.01~"not oscillating",
                                 JTK.cytof$median_expression<0.01~"no/low expression")
ggplot(JTK.cytof)+geom_tile(aes(x=CycID_RNA,y=celltype,fill=significant),color="black",linewidth=0.5)+
  scale_fill_manual(values=c("white","grey","red"))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,color="black",size=12),
        axis.text.y=element_text(color="black",size=10),legend.text=element_text(size=12),
        panel.border = element_rect(fill=NA,color="black"),legend.title = element_text(size=14),
        plot.margin = margin(0,0.5,0,0,"cm"))+
  xlab("")+ylab("")+guides(fill=guide_legend(title=NULL,position = "bottom"))

# Fig. 4B
plotPeakAmongIndividuals(JTK_result_all,gene = "CXCR4",float.barlength = 0.3)

# Fig. 4C
plotPeakAmongIndividuals(JTK_result_all,gene = "ITGB2",float.barlength = 0.3)

# Fig. 4D
JTK_result_all<-readRDS('~/analysis/circadian/R/JTK.16individual.cor2batch0.2.fold2bkg1.5.fold2trough1.5.rds')
ggplot(JTK_result_all)+geom_bar(aes(x=LAG,fill=celltype))+
  facet_wrap(~celltype,scales = "free_y",ncol = 6,nrow=3)+
  scale_fill_manual(values=generateColor(20,alpha = 0.75))+NoLegend()+
  xlab("circadian time")+ylab("count of peaking genes")

# Fig. 4E
plotPhaseofCellType(JTK_result_all,colors =generateColor(20,alpha = 0.75))+NoLegend()

# Fig. XX
n.count.osc.by.age<-JTK_result_all %>% group_by(individual,celltype) %>% summarise(n.count.oscillating=n())
n.count.osc.by.age<-n.count.osc.by.age %>% group_by(individual) %>% mutate(.,n_total=sum(n.count.oscillating))
n.count.osc.by.age$n.count.oscillating.norm<-n.count.osc.by.age$n.count.oscillating/n.count.osc.by.age$n_total
n.count.osc.by.age$age<-AGES[n.count.osc.by.age$individual]
n.count.osc.by.age$group<-ifelse(n.count.osc.by.age$age>30,"middle","young")
n.group<-n.count.osc.by.age[c("celltype","group")] %>% unique() %>% group_by(celltype) %>% summarise(n.group=n())
valid.celltype<-n.group[n.group$n.group>1,"celltype"]$celltype
n.count.osc.by.age<-n.count.osc.by.age[n.count.osc.by.age$celltype %in% valid.celltype,]
# Perform Wilcoxon test for each celltype
library(rstatix)
stat.test <- n.count.osc.by.age %>%
  group_by(celltype) %>%
  pairwise_wilcox_test(
    n.count.oscillating.norm ~ group,  # Formula: Value ~ Group
    p.adjust.method = "bonferroni",   # Adjust for multiple testing
    paired = FALSE                    # Assuming your groups are independent (different individuals)
  ) %>%
  add_xy_position(x = "celltype", group = "group", dodge = 0.8) # Crucial for positioning!
p<-ggplot(n.count.osc.by.age)+geom_boxplot(aes(x=celltype,y=n.count.oscillating.norm,fill=group),position="dodge",outlier.colour=NA)+
  theme_half_open()+theme(axis.text.x=element_text(angle=60,hjust=1),plot.margin = margin(1,0.2,0.2,0.2,"cm"))+xlab("")+ylab("number of oscillating genes")+
  scale_y_continuous(expand = c(0,0))
final_plot <- p +
  stat_pvalue_manual(
    stat.test,
    label = "p.adj.signif",    # Show asterisks (***, **, *, ns)
    tip.length = 0.01,         # Short bracket tips
    hide.ns = FALSE,           # Set to TRUE to hide 'ns' labels
    bracket.nudge.y = 0.5      # Nudge brackets higher if they overlap
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) # Add space on top for brackets

final_plot
