source('~/script/circadian/circadian_core.R')

# FIG. 1A study design
# FIG. 1B individual metadata

# FIG. 1C UMAP of PBMC atlas
test<-readRDS("~/analysis/core_data/16individuals_annotated_integrated.rds")
test$type<-test$predicted.celltype.l2
test<-subset(test,type!="Eryth"&type!="Platelet")
celltypes<-test$type %>% unique() %>% sort()
cell_legend_data<-data.frame("type"=celltypes,"pos.x"=rep(1,length(celltypes)),
                             "pos.y"=rev(1:length(celltypes)),"ID"=1:length(celltypes))
mappings<-cell_legend_data[,c("type","ID")]
colnames(mappings)[2]<-"new_label"
test<-addLabels(test,mappings)
dimplot_publication(test,reduction = "umap",reduction.name="umap.harmony",group.by = "new_label",label=T,
                    colors=generateColor(30,col.dist.min=0.3,seed = 10))+NoLegend()
ggsave("/lustre/home/acct-medll/medll/figures/circadian/Fig1C.png",width=6,height=6)
ggplot(cell_legend_data)+geom_point(aes(x=pos.x,y=pos.y,fill=type),shape=21,size=8,color="black")+
  geom_text(aes(x=pos.x,y=pos.y,label=ID),color="white")+
  geom_text(aes(x=pos.x+0.075,y=pos.y,label=type),hjust="left")+
  scale_x_continuous(limits=c(1,2,3,4))+scale_fill_manual(values=generateColor(30,col.dist.min=0.3,seed = 10))+
  theme_nothing()
ggsave("/lustre/home/acct-medll/medll/figures/circadian/Fig1C_legend.pdf",width=3,height=6)
rm(test)
gc()

#FIG. 1D composition of cells at each CT and their markers
test<-readRDS("~/analysis/core_data/16individual.srt.annotated.rds")
test$type<-test$predicted.celltype.l2
test<-subset(test,type!="Eryth"&type!="Platelet")
meta<-test@meta.data
meta$CT<-factor(meta$CT,levels=sort(unique(meta$CT)))
meta$type<-meta$predicted.celltype.l2
celltypes<-sort(unique(meta$type))
meta$type<-factor(meta$type,levels=rev(celltypes))

plot1<-ggplot(meta)+geom_bar(aes(x=type,fill=CT),position="fill",color="black")+
  scale_fill_d3(alpha=0.75,palette = "category20")+scale_y_continuous(expand = c(0,0))+
  theme(plot.margin=margin(0.5,0,0,0,"cm"),
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.text.y=element_text(color="black"))+
  ylab("")+xlab("")+NoLegend()

plot1_legend<-ggplot(meta)+geom_bar(aes(x=type,fill=CT),position="fill",color="black")+
  scale_fill_d3(alpha=0.75,palette = "category20")+scale_y_continuous(expand = c(0,0))+
  theme(plot.margin=margin(0.5,0,0,0,"cm"),
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.text.y=element_text(color="black"))+
  ylab("")+xlab("")
plot1_legend<-get_legend(plot1_legend)
grid.draw(plot1_legend)
test<-NormalizeData(test)
plot2<-DotPlot(test,rev(c("AXL","SIGLEC6","MS4A1","CD27",
                          "CD14","FCGR3A","CD3D","CD4","CD8A","FOXP3",
                          "CD1C","CLEC9A","HBA1","TRGC1","TRDC","CD34","SLC4A10",
                          "GNLY","GZMK","GZMH","LILRA4","JCHAIN","PPBP")),group.by="type",
                    cols = c("#507AAF","#BE5C37"))+scale_y_discrete(limits=sort(unique(test$type)))+
  scale_size_continuous(range = c(0,4))+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=9),
        plot.margin=margin(0,0,0,0,"cm"),axis.text.y=element_text(size=9),
        panel.border = element_rect(color = "black",linewidth = 1),
        legend.text = element_text(size=10),legend.title=element_text(size=10))+
  xlab("")+ylab("")+coord_flip()

ggarrange(plot1,plot2,nrow=2,align="v",heights = c(1,5))
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig1D.pdf',width = 6,height=6)


# FIG. 1E confusion matrix, CYTOF annotation vs Azimuth
## method 1, by expression correlation between protein and RNA
test<-readRDS("~/analysis/core_data/4multiplexedindividual.srt.annotated.rds")
mat<-srt2bulkMatrix(test,split.by = "predicted.celltype.l2",layer = "count")
mat<-mat[CYTOF2RNA_FEATURES %>% as.vector(),]
mat<-mat[rownames(mat)!="CD151",]
mat<-mat[rownames(mat)!="NA",]
mat<-mat[,colnames(mat)!="Platelet"]


test.cytof<-readRDS('~/analysis/core_data/CYTOF.300k.reannotation.srt.rds')
test.cytof@meta.data$type<-test.cytof@meta.data$manual.new
mat.cytof<-srt2bulkMatrix(test.cytof,split.by = "type",layer = "count")
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
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig1E.pdf',width = 7,height=4)

# FIG. S2A, plot showing low expression level of genes
test.annotated<-readRDS("~/analysis/core_data/16individual.srt.annotated.rds")
test.annotated$type<-test.annotated$predicted.celltype.l2
count.mat<-LayerData(test.annotated,features=CIRCADIAN_GENES_MAIN,layer = "count") %>% as.data.frame()
count.mat<-rownames_to_column(count.mat,var = "features")
count.mat.long<-gather(count.mat,key="observations",value="counts",-features) 
meta.mat<-test.annotated@meta.data[,c("nCount_RNA","type")]
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
ggsave('/lustre/home/acct-medll/medll/figures/circadian/FigS2A.pdf',width = 8,height=5)


# FIG. S2B
# circadian metacell size benchmark
# compare metacell proportion 0.08 0.04 0.03
plotdata=NULL
for (proportion in c(0.03,0.04,0.08)) {
  AllJTKresult.filtered=readRDS(paste0("~/analysis/core_data/JTK.result.filtered.addp2bkg.addcor2batch.bytype.",proportion,".rds"))
  AllJTKresult.filtered=AllJTKresult.filtered[AllJTKresult.filtered$fold_to_background>=2&AllJTKresult.filtered$cor_to_batch<0.2,]
  #AllJTKresult.filtered<-AllJTKresult.filtered[AllJTKresult.filtered$cor_to_batch<0.2,]
  AllJTKresult.filtered<-AllJTKresult.filtered[!is.na(AllJTKresult.filtered$celltype) & !is.na(AllJTKresult.filtered$individual), ]
  this.plotdata=AllJTKresult.filtered %>% group_by(celltype) %>% summarise(n.gene.found=n())
  this.plotdata$ratio.shrink=proportion
  if(is.null(plotdata)){
    plotdata=this.plotdata
  }else{
    plotdata=rbind(plotdata,this.plotdata)
  }
}
plotdata$ratio.shrink<-as.factor(plotdata$ratio.shrink)
plotdata$meta.cell.size<-case_when(plotdata$ratio.shrink==0.03~30,
                                   plotdata$ratio.shrink==0.04~20,
                                   plotdata$ratio.shrink==0.08~10)
plotdata$meta.cell.size<-as.factor(plotdata$meta.cell.size)
plotdata$celltype<-factor(plotdata$celltype,levels=unique(plotdata$celltype))
ggplot(plotdata)+geom_bar(aes(x=celltype,y=n.gene.found,fill=meta.cell.size),color="black",
                          position=position_dodge2(width = 2,preserve = "single"),stat="identity")+
  scale_y_continuous(expand = c(0,0))+scale_fill_aaas()+xlab("")+
  ggtitle("number of circadian genes found\nunder different metacell size")+
  theme_minimal()+theme(axis.text.x=element_text(angle = 60,hjust=1),plot.title=element_text(hjust = 0.5))
ggsave('/lustre/home/acct-medll/medll/figures/circadian/FigS2B.pdf',width = 6,height=3)


# FIG. S2C
srt.metacell<-readRDS("~/analysis/core_data/seacell.16individual.TF.rds")
# part1
part1<-VlnPlot(srt.metacell,CIRCADIAN_GENES_MAIN,group.by = "type",stack = T,flip = T,cols = generateColor(15,seed = 2025))+NoLegend()+xlab("")+theme(strip.text=element_text(size=10,face="plain"))
# part2
plotdata<-srt.metacell$type %>% table() %>% as.data.frame()

part2<-ggplot(plotdata)+
  geom_text(data=plotdata[plotdata$Freq<2500,],aes(x=get("."),y=Freq,label=Freq),hjust=-0.25,size=4,angle=90)+
  geom_text(data=plotdata[plotdata$Freq>=2500,],aes(x=get("."),y=Freq/2,label=Freq),size=4,angle=90)+
  geom_bar(aes(x=get("."),y=Freq),stat="identity",fill="darkblue",alpha=0.5,color="black")+
  theme_minimal_hgrid(font_size = 12)+theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
  xlab("")+ylab("number of metacells")+scale_y_continuous(limits = c(0,10000),expand = c(0,0))+xlim(plotdata$.)
part2/part1
ggsave('/lustre/home/acct-medll/medll/figures/circadian/FigS2C.pdf',width = 7,height=5)

# FIG. S2DEFG ?
test.table<-readRDS('~/analysis/core_data/16individual.pseudobulk.bytype.CPM.rds')
srt.metacell<-readRDS("~/analysis/core_data/seacell.0.04.16individual.RawCount.srt.rds")
imposeMetaCellOnPseudoBulk(srt.metacell = srt.metacell,char.cell.type="CD14 Mono",char.feature="BMAL1",
                           df.pseudobulk = test.table)
imposeMetaCellOnPseudoBulk(srt.metacell = srt.metacell,char.cell.type="CD14 Mono",char.feature="PER2",
                           df.pseudobulk = test.table)
imposeMetaCellOnPseudoBulk(srt.metacell = srt.metacell,char.cell.type="CD14 Mono",char.feature="CRY1",
                           df.pseudobulk = test.table)
imposeMetaCellOnPseudoBulk(srt.metacell = srt.metacell,char.cell.type="CD14 Mono",char.feature="NR1D2",
                           df.pseudobulk = test.table)

# FIG.2A figure elaborating principle of metacells

# FIG.2B
meta2d.result0<-readRDS('~/analysis/core_data/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')

plot1<-plotCountOscilattingByCelltypePseudobulk(meta2d.result0,cut.off = 0.05,flip = T,min.occur = 3,metacell.mode = T)

meta2d.result<-meta2d.result0[(meta2d.result0$CycID %in% CIRCADIAN_GENES_MAIN) & meta2d.result0$ADJ.P<0.05,]
meta2d.result<-meta2d.result %>% group_by(CycID,celltype) %>% mutate(log.p.adj=median(-log10(ADJ.P)),n.individual=n())
meta2d.result$significant<-ifelse(meta2d.result$log.p.adj>=2,"YES","NO")
meta2d.result$celltype<-gsub("_"," ",meta2d.result$celltype)
plot2<-ggplot(meta2d.result)+geom_point(aes(x=CycID,y=celltype,size=n.individual,fill=log.p.adj,color=significant),shape=21,stroke=0.75)+theme_classic()+
  scale_fill_gradient2(low = "lightblue",high = "darkred",name="-log10(p.adj)")+scale_color_manual(values=c("white","black"),name="p<=0.01")+
  scale_y_discrete(limits=rev(CELL_TYPES_METACELL))+
  theme(axis.text.x=element_text(angle=60,hjust=1,color="black"),legend.box.margin=margin(),
        axis.text.y=element_text(color="black"),legend.margin=margin(),
        legend.position="left")+xlab("")+ylab('')
plot2|plot1
ggsave("~/figures/circadian/FIG2B.pdf",width=6,height=4)

# FIG.2C
meta2d.result0<-readRDS('~/analysis/core_data/meta2d.pseudobulkByCelltypeByIndividual.rds')
plot.list<-NULL
for(feature in c("BMAL1","CLOCK","PER1","PER2","CRY1","CRY2","NR1D1","NR1D2","DBP")){
  plot.list[[feature]]=plotPseudobulk(char.cell.type = "CD14 Mono",vector.features=feature,
                                      char.data.from.mat=test.table,df.meta2d=meta2d.result0,
                                      bool.merge = T)+xlab("")
}
ggarrange(plotlist = plot.list,ncol = 3,nrow = 3)
ggsave("~/figures/circadian/FIG2C.pdf",width=6.5,height=6.5)

# FIG.2D
test.metacell.TF<-readRDS('~/analysis/core_data/seacell.16individual.TF.rds')
# Get the cells that don't have NA in the type column
cells_to_keep <- !is.na(test.metacell.TF$type)
test.metacell.TF <- test.metacell.TF[, cells_to_keep]
DefaultAssay(test.metacell.TF)<-"TF"
TF.activity<-LayerData(test.metacell.TF,layer = "data")
TF.activity<-rownames_to_column(as.data.frame(TF.activity),"TF_name")
TF.activity<-gather(TF.activity,key="metacell",value="TF_activity",-TF_name)
#TF.activity$individual<-getField(TF.activity$metacell,"_",2)
#TF.activity$CT<-getField(TF.activity$metacell,"_",3)
TF.activity$celltype<-getField(TF.activity$metacell,"_",4)
TF.activity<-TF.activity %>% group_by(TF_name,celltype) %>% mutate(n.cell=n(),median_activity=median(TF_activity))
TF.activity<-TF.activity[TF.activity$n.cell>2,]
TF.activity<-arrange(TF.activity,by=desc(median_activity))
TF.activity$celltype<-gsub("-"," ",TF.activity$celltype)
TF.activity$rough_type<-CELL_TYPES2[TF.activity$celltype]
ggplot(TF.activity[TF.activity$TF_name=="CLOCK",])+
  geom_boxplot(aes(x=fct_inorder(celltype),y=TF_activity,fill=rough_type),outliers = F,alpha=0.5)+
  theme_half_open()+ggtitle(label = "",subtitle = "TF activity of CLOCK")+
  theme(axis.text.x = element_text(angle=60,hjust=1))+ylab("intensity of activity")+
  xlab("")+scale_fill_manual(name="",values = generateColor(n = length(unique(TF.activity$rough_type))))
ggsave("~/figures/circadian/FIG2E.pdf",width=6.5,height=3.5)
rm(test.metacell.TF)
gc()

# FIG.2E
meta2d.result0<-readRDS('~/analysis/core_data/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
meta2d.result<-meta2d.result0[meta2d.result0$ADJ.P<0.05&meta2d.result0$fold_peak2trough>=1.5,]
meta2d.result<-meta2d.result[!grepl("ENSG",meta2d.result$CycID),]
meta2d.result<-meta2d.result %>% group_by(CycID) %>% mutate(n.occurance=n())
top100<-meta2d.result[c("CycID","n.occurance")] %>% unique() %>% arrange(by=desc(n.occurance)) %>% head(100)
plotdata<-meta2d.result[meta2d.result$CycID %in% top100$CycID,]
plotdata<-plotdata %>% group_by(CycID) %>% mutate(median.phase=median(LAG),mean.log10.p.adj=mean(-log10(ADJ.P)))
plotdata2<-plotdata[c("CycID","n.occurance","mean.log10.p.adj")] %>% unique() %>% arrange(desc(n.occurance))
plot1<-ggplot(plotdata)+geom_bar(aes(x=fct_infreq(CycID),fill=median.phase),color="black")+
  geom_text(data=plotdata2,aes(x=CycID,y=n.occurance,label=CycID),angle=60,vjust=0.5,hjust=-0.11,size=3)+
  scale_fill_gradientn(name="peaking time",colours=c("black","darkgrey","lightgrey","white","lightgrey","darkgrey","black"))+
  theme_half_open()+scale_y_continuous(expand = c(0,0),limits = c(0,130))+scale_x_discrete(expand=c(0,2.5))+
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
  scale_fill_gradient2(low="darkblue",high = "darkred",name="mean -log10(p.adj)",midpoint = -log10(0.05))+
  xlab("")+ylab("")
legend2<-get_legend(plot2)
legend.plot<-ggarrange(as_ggplot(legend1),as_ggplot(legend2),ncol=2,nrow=1)
plot1<-plot1+NoLegend()
plot2<-plot2+NoLegend()
ggarrange(plot1,plot2,nrow=2,ncol=1,align="v",heights = c(6,1))+
  patchwork::inset_element(legend.plot, left = 0.7, bottom = 0.8, right = 0.95, top = 0.95)
ggsave("~/figures/circadian/FIG2E.pdf",width=17,height=3.5)

# FIG.2F dotplot of bulk RNA seq, x=log10(amplitude), y=-log10(p.adj)
data<-read.delim('~/analysis/core_data/tauFisher_test/GSE113883_hs_full_adj_meta2d/JTKresult_GSE113883_hs_full_adj.txt')
data<-data[data$ADJ.P<0.05,]
plot1<-ggplot(data)+geom_point(aes(x=log10(AMP),y=-log10(ADJ.P),fill=LAG),shape=21,size=3)+
  geom_text_repel(aes(x=log10(AMP),y=-log10(ADJ.P),label=CycID),hjust=0,vjust=-1)+
  scale_fill_gradientn(colours = c("white","lightgrey","darkgrey","black"),name="relative\npeaking time")+theme_half_open()+
  ylab("-log10(p.adj)")+xlab("log10(amplitude)")+
  theme(legend.position = c(0.75,0.9),legend.direction = "horizontal",
        legend.title.position = "top",plot.margin = margin(0,0,0,0))

plot2<-ggplot(data)+geom_point(aes(x=1,y=-log10(ADJ.P),fill=LAG),shape=21,size=3)+
  geom_text(data=data[data$CycID %in% c("DDIT4","FKBP5","RBFOX2","RBM3","CXCR4"),],aes(x=1.2,y=-log10(ADJ.P),label=CycID),hjust=0)+
  scale_fill_gradientn(colours = c("white","lightgrey","darkgrey","black"))+theme_nothing()+
  ylab("-log10(p.adj)")+NoLegend()+theme(plot.margin = margin(0,1,0,0,"cm"))+scale_x_continuous(expand=c(0,0),limits=c(0.8,2))

ggarrange(plot1,plot2,ncol = 2,nrow=1,align="hv",widths = c(3,1))
ggsave("~/figures/circadian/FIG2F.pdf",width=8,height=5)

# FIG.2G correlation of hormone and RNA expression

# FIG.2H RNA velosity of circadian genes
spliced.meta2d<-readRDS('~/analysis/core_data/JTK.spliced.filtered.16individual.addmedianexp.addp2t.rds')
spliced.meta2d<-spliced.meta2d[spliced.meta2d$ADJ.P<0.05&spliced.meta2d$fold_peak2trough>1.5,]

unspliced.meta2d<-readRDS('~/analysis/core_data/JTK.unspliced.filtered.16individual.addmedianexp.addp2t.rds')
unspliced.meta2d<-unspliced.meta2d[unspliced.meta2d$ADJ.P<0.05&unspliced.meta2d$fold_peak2trough>1.5,]

holistic.meta2d<-readRDS('~/analysis/core_data/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
holistic.meta2d<-holistic.meta2d[holistic.meta2d$ADJ.P<0.05&holistic.meta2d$fold_peak2trough>1.5,]

spliced.genes<-unique(spliced.meta2d$CycID)
unspliced.genes<-unique(unspliced.meta2d$CycID)
holistic.genes<-unique(holistic.meta2d$CycID)
my_list <- list(spliced = spliced.genes, unspliced = unspliced.genes, holistic = holistic.genes)

plotVenn(my_list)
ggsave("~/figures/circadian/FIG2H.pdf",width=5,height=5)

# FIG.2I compare peaking time
spliced.summa<-spliced.meta2d %>% group_by(CycID,individual) %>% summarise(.,phase=median(LAG))
colnames(spliced.summa)[3]<-"spliced"
unspliced.summa<-unspliced.meta2d %>% group_by(CycID,individual) %>% summarise(.,phase=median(LAG))
colnames(unspliced.summa)[3]<-"unspliced"
phases<-left_join(spliced.summa,unspliced.summa,by=c("CycID","individual"))
phases<-phases[!is.na(phases$spliced)&!is.na(phases$unspliced),]
phases$diff<-phases$unspliced-phases$spliced
phases$spliced<-ifelse(phases$diff>12,phases$spliced+24,phases$spliced)
phases$spliced<-ifelse(phases$diff< -12,phases$spliced-24,phases$spliced)
phases$diff<-NULL
phases<-gather(phases,key="status",value="phase",-CycID,-individual)


ggplot(phases,aes(x=status,y=phase,fill=status))+geom_boxplot(outliers = F)+
  scale_x_discrete(limits=c("unspliced","spliced"))+
  stat_compare_means(method = "wilcox.test",comparisons = list(c("unspliced", "spliced")),label = "p.signif")+
  theme_half_open()+theme(axis.text.x=element_text(angle=60,hjust=1))+NoLegend()+ylab("peaking time")+xlab("")
ggsave("~/figures/circadian/FIG2I.pdf",width=2.5,height=4)

# FIG.3A dot plot showing genes in different cell type
meta2d.result0<-readRDS('~/analysis/core_data/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
#meta2d.result<-meta2d.result0[meta2d.result0$meta2d_pvalue<=0.1&!grepl("ENSG|LINC",meta2d.result0$CycID),]
meta2d.result<-filterMeta2dTable(meta2d.result0,float.pvalue.s1 = 0.05,float.pvalue.s2 = 0.05,int.min.individual = 2,mode = "JTK")

plot1<-plotOscillatingGeneSummarise(meta2d.result,vec.celltype="CD14 Mono",xlim = c(0,0.6),ylim = c(1,8),
          x.axis="phase.consist",y.axis="fold.peak2trough",method = "JTK",int.color.index = 5)+NoLegend()
plot2<-plotOscillatingGeneSummarise(meta2d.result,vec.celltype="CD4 Naive",xlim = c(0,0.6),ylim = c(1,8),
          x.axis="phase.consist",y.axis="fold.peak2trough",method = "JTK",int.color.index = 8)+NoLegend()
plot3<-plotOscillatingGeneSummarise(meta2d.result,vec.celltype="CD8 Naive",xlim = c(0,0.6),ylim = c(1,8),
          x.axis="phase.consist",y.axis="fold.peak2trough",method = "JTK",int.color.index = 4)
plot1|plot2|plot3
ggsave("~/figures/circadian/FIG3A.pdf",width=12,height=4)

# FIG.3B heatmap showing shared genes between celltypes
#meta2d.result<-meta2d.result0[meta2d.result0$meta2d_pvalue<=0.1&!grepl("ENSG|LINC",meta2d.result0$CycID),]
meta2d.result0<-readRDS('~/analysis/core_data/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
#meta2d.result<-meta2d.result0[meta2d.result0$meta2d_pvalue<=0.1&!grepl("ENSG|LINC",meta2d.result0$CycID),]
meta2d.result<-filterMeta2dTable(meta2d.result0,float.pvalue.s1 = 0.05,float.pvalue.s2 = 0.05,int.min.individual = 2,mode = "JTK")
meta2d.result<-meta2d.result %>% group_by(CycID,celltype) %>% mutate(n.individual=n())
#meta2d.result<-meta2d.result[meta2d.result$n.individual>=2,]
meta2d.result$celltype<-gsub("_"," ",meta2d.result$celltype)
celltypes<-unique(meta2d.result$celltype)
meta2d.result<-meta2d.result[,c("CycID","celltype")] %>% unique
meta2d.result<-as.data.frame(meta2d.result)
share.data=NULL
for(type.this in celltypes){
  for(type.that in celltypes){
    n.gene.type.this=nrow(meta2d.result[meta2d.result$celltype==type.this,])
    n.share.gene=length(intersect(meta2d.result[meta2d.result$celltype==type.this,"CycID"],
                                  meta2d.result[meta2d.result$celltype==type.that,"CycID"]))
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
invalid.type=share.data[is.na(share.data$shared.ratio)|share.data$from.type=="Platelet","from.type"]
share.data$from.type<-gsub("_"," ",share.data$from.type)
share.data$to.type<-gsub("_"," ",share.data$to.type)
ggplot(share.data[!is.na(share.data$shared.ratio)&!(share.data$from.type%in%invalid.type)&!(share.data$to.type%in%invalid.type),])+
  geom_raster(aes(x=from.type,y=to.type,fill=shared.ratio))+scale_x_discrete(limits=CELL_TYPES_METACELL)+scale_y_discrete(limits=CELL_TYPES_METACELL)+
  scale_fill_gradientn(colours = c("#283168","#009fc0","#f0f0cc","white"),values = c(0,0.5,0.9,1))+
  theme(axis.text.x = element_text(angle=60,hjust=1),plot.margin=margin(0,0,0,0,"cm"))+ylab("")+xlab("")
ggsave("~/figures/circadian/FIG3B.pdf",width=7,height=6)

# FIG.3C gene enrichment of top100 shared genes
meta2d.result0<-readRDS('~/analysis/core_data/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
meta2d.result<-filterMeta2dTable(meta2d.result0,float.pvalue.s1 = 0.05,float.pvalue.s2 = 0.05,int.min.individual = 2,mode = "JTK")
meta2d.result<-meta2d.result %>% group_by(CycID,celltype) %>% mutate(n.individual=n())
meta2d.result<-meta2d.result[meta2d.result$n.individual>=2,]
meta2d.result<-meta2d.result[c("CycID",'celltype')] %>% unique()
genes.shared<-meta2d.result$CycID %>% table() %>% sort() %>% tail(100) %>% names()
result<-enrichGObyHGNC(genes.shared)
plotdata<-result@result[result@result$pvalue<0.01,] %>% arrange(desc(pvalue)) %>% tail(50)
plotdata<-plotdata[!duplicated(plotdata$geneID),]
ggplot(plotdata)+geom_point(aes(x=fct_inorder(Description),color=-log10(pvalue),
                                y=-log10(pvalue),size=Count))+coord_flip()+
  xlab("")+ylim(c(1.95,4.5))+theme_classic2(base_size = 13)+
  theme(legend.position = "bottom",legend.box = "vertical",plot.margin = margin(0.5,2,0.5,0,"cm"))
ggsave("~/figures/circadian/FIG3C.pdf",width=8,height=7)

# FIG.3D phase of FKBP5 and DDIT4
meta2d.result0<-readRDS('~/analysis/core_data/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
meta2d.result<-filterMeta2dTable(meta2d.result0,float.pvalue.s1 = 0.05,float.pvalue.s2 = 0.05,int.min.individual = 2,mode = "JTK")
plot1<-plotPeakAmongIndividuals(meta2d.result,mode = "JTK",gene = "DDIT4")
plot2<-plotPeakAmongIndividuals(meta2d.result,mode = "JTK",gene = "FKBP5")
plot1/plot2
ggsave("~/figures/circadian/FIG3D.pdf",width=5,height=5)

# FIG.3E count of CR genes ant different CT
plotCountOscilattingByCelltypeByCT(meta2d.result,mode = "JTK",float.pvalue = 0.1,
                  vec.celltype.to.plot=c("CD14 Mono", "CD4 Naive","CD8 Naive" ,"NK"))
ggsave("~/figures/circadian/FIG3E.pdf",width=5,height=5)


# FIG. 4A
pseudobulkdata<-readRDS('~/analysis/core_data/16individual.pseudobulk.bytype.CPM.rds')
metacell<-readRDS("~/analysis/core_data/seacell.0.04.16individual.RawCount.srt.rds")
imposeMetaCellOnPseudoBulk(srt.metacell = metacell,char.cell.type = "CD14 Mono",char.feature = "NR1D2",
                           vector.individual = c("LJQ","YXQ"),df.pseudobulk = pseudobulkdata)
plotPseudobulk(char.cell.type = "CD14 Mono",vector.features = c("BMAL1","NR1D2"),
      vector.individual = c("SLY","LJQ","YXQ"),char.data.from.mat = pseudobulkdata,bool.relative = T)+
  ggtitle("BMAL1 and NR1D2 expression in CD14+ Monocytes")+xlab("")
ggsave("~/figures/circadian/FIG4A.pdf",width=7,height=2.5)

# FIG.4B prediction with core circadian genes
bulk_by_type<-readRDS("~/analysis/core_data/16individual.pseudobulk.bytype.CPM.rds")
top.genes<-CIRCADIAN_GENES_MAIN
predicted.data<-predictIndividualwithTauFisher(data=bulk_by_type,individual = c("KD","JJC","HZD","ZYX","WLG","XSP"),celltype = "CD14_Mono",
                                               predictor.genes = top.genes,bool.include.trained = T,bool.return.pred = T)
predicted.data<-predicted.data %>% group_by(individual) %>% mutate(median_diff=median(diff),RMSE=mean(diff^2)^(1/2))
individual_order<-unique((predicted.data %>% arrange(median_diff))$individual)
ggplot(predicted.data)+geom_boxplot(aes(x=individual,y=diff,fill=RMSE),outliers = F)+
  geom_hline(yintercept = 0,linetype=2)+ggtitle("difference between\nactual time and predicted time")+
  theme_half_open()+scale_fill_gradientn(colors=c("#2E86AB","white","#A23B72"))+
  scale_x_discrete(limits=individual_order,labels=INDIVIDUAL_MASKING[individual_order])+
  theme(axis.text.x=element_text(angle=60,hjust=1),plot.title = element_text(hjust=0.5,size=12))+ylab("")
ggsave("~/figures/circadian/FIG4B.pdf",width=5,height=3)

# FIG.4C dot plot of differential gene analysis
test.deseq<-readRDS('~/analysis/R/DEseq2/CD14_Mono_normal_vs_disordered.DESeq2.obj.rds')
DESeq2DotPlot(test.deseq,char.title = "CD14 Monocyte\nnormal vs disordered")
ggsave("~/figures/circadian/FIG4C.pdf",width=4,height=5)

# FIG.4D ADM, disordered vs normal
plotDESeq2Expression(DESeq2.obj = test.deseq,"ADM")
ggsave("~/figures/circadian/FIG4D.pdf",width=3,height=4)

# FIG.S4A dotplot of circaidan gene significance in different individual and celltype
plotMainCircadianGeneByIndividualByType(meta2d.result0,mode = "JTK",cells.to.plot=CELL_TYPES_METACELL)
ggsave("~/figures/circadian/FIGS4A.pdf",width=12,height=10)

# FIG.S4B tauFisher benchmarking with all individual
plotdata=NULL
for (n.individual in 1:8) {
  for (rep.id in 1:5) {
    #message(n.top.genes)
    all.individual=unique(meta2d.result0$individual)
    test.individuals=sample(all.individual,n.individual,replace=F)
    #top.genes=(JTK_result_used[c("CycID","n.individual")] %>% unique() %>% arrange(desc(n.individual)) %>% head(.,10))$"CycID"
    #top.genes=(JTK_result_used[c("CycID","median.percentile")] %>% unique() %>% arrange(median.percentile) %>% head(.,n.top.genes))$"CycID"
    top.genes=c("BMAL1","CLOCK","PER1","PER2","CRY1","CRY2","NR1D1","NR1D2","DBP")
    #top.genes=c(top.genes,c("BMAL1","CLOCK","PER1","PER2","CRY1","CRY2","NR1D1","NR1D2","DBP")) %>% unique()
    res=predictIndividualwithTauFisher(data=bulk_by_type,individual = all.individual[!(all.individual %in% test.individuals)],celltype = "CD14_Mono",
                                       predictor.genes = top.genes,return.res = T)
    this.plotdata=data.frame("n.individual"=n.individual,"accuracy"=res$Accuracy,"RMSE"=res$RMSE,"replicate"=rep.id)
    if(is.null(plotdata)){
      plotdata=this.plotdata
    }else{
      plotdata=rbind(plotdata,this.plotdata)
    }
  }
}
df_long <- plotdata %>%
  pivot_longer(cols = c(accuracy, RMSE), 
               names_to = "metric", 
               values_to = "value")

# For separate panels
ggplot(df_long, aes(x = factor(n.individual), y = value, fill = metric)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    width = 0.6
  ) +
  geom_point(
    position = position_jitter(width = 0.2),
    size = 1.5,
    alpha = 0.6,
    aes(color = metric)
  ) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("accuracy" = "#2E86AB", "RMSE" = "#A23B72")) +
  scale_color_manual(values = c("accuracy" = "#2E86AB", "RMSE" = "#A23B72")) +
  labs(
    x = "Number of Individual used",
    y = "Value"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(),
    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2)
  )
ggsave("~/figures/circadian/FIGS4B.pdf",width=3,height=4)

# FIG.S4C tauFisher benchmarking with normal circaidan pattern individual 
plotdata=NULL
for (n.individual in 1:8) {
  for (rep.id in 1:5) {
    #message(n.top.genes)
    all.individual=names(INDIVIDUAL2CIRCADIAN_GROUP2[INDIVIDUAL2CIRCADIAN_GROUP2=="normal"])
    test.individuals=sample(all.individual,n.individual,replace=F)
    #top.genes=(JTK_result_used[c("CycID","n.individual")] %>% unique() %>% arrange(desc(n.individual)) %>% head(.,10))$"CycID"
    #top.genes=(JTK_result_used[c("CycID","median.percentile")] %>% unique() %>% arrange(median.percentile) %>% head(.,n.top.genes))$"CycID"
    top.genes=c("BMAL1","CLOCK","PER1","PER2","CRY1","CRY2","NR1D1","NR1D2","DBP")
    #top.genes=c(top.genes,c("BMAL1","CLOCK","PER1","PER2","CRY1","CRY2","NR1D1","NR1D2","DBP")) %>% unique()
    res=predictIndividualwithTauFisher(data=bulk_by_type,individual = all.individual[!(all.individual %in% test.individuals)],celltype = "CD14_Mono",
                                       predictor.genes = top.genes,return.res = T)
    this.plotdata=data.frame("n.individual"=n.individual,"accuracy"=res$Accuracy,"RMSE"=res$RMSE,"replicate"=rep.id)
    if(is.null(plotdata)){
      plotdata=this.plotdata
    }else{
      plotdata=rbind(plotdata,this.plotdata)
    }
  }
}
df_long <- plotdata %>%
  pivot_longer(cols = c(accuracy, RMSE), 
               names_to = "metric", 
               values_to = "value")

# For separate panels
ggplot(df_long, aes(x = factor(n.individual), y = value, fill = metric)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    width = 0.6
  ) +
  geom_point(
    position = position_jitter(width = 0.2),
    size = 1.5,
    alpha = 0.6,
    aes(color = metric)
  ) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("accuracy" = "#2E86AB", "RMSE" = "#A23B72")) +
  scale_color_manual(values = c("accuracy" = "#2E86AB", "RMSE" = "#A23B72")) +
  labs(
    x = "Number of Individual used",
    y = "Value"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(),
    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2)
  )
ggsave("~/figures/circadian/FIGS4C.pdf",width=3,height=4)

# FIG.S4D age and circadian rhythm
meta2d.result<-filterMeta2dTable(meta2d.result0,float.pvalue.s1=0.05,float.pvalue.s2=0.05,mode = "JTK")
meta2d.result$age<-AGES[meta2d.result$individual]
meta2d.result$age_group<-ifelse(meta2d.result$age<40,"young","aged")
meta2d.result<-meta2d.result[!is.na(meta2d.result$age_group),]
meta2d.result<-meta2d.result[meta2d.result$fold_peak2trough<3,]
# Create the boxplot with statistical annotations
ggplot(meta2d.result, aes(x = age_group, y = fold_peak2trough)) +
  geom_boxplot(aes(fill = age_group), width = 0.6, alpha = 0.8, outliers = F) +
  stat_compare_means(
    method = "wilcox.test", # or "wilcox.test" for non-parametric
    comparisons = list(c("young", "aged")), # adjust based on your groups
    label = "p.signif", # shows asterisks
    tip.length = 0.01,
    vjust = 0.5) + scale_fill_manual(values = c("#1F77B4", "#FF7F0E"))+
  ylim(c(1,3.5))+
  labs(x = "",y = "relative amplitude",title = "amplitude of oscillating genes\nbetween age groups") +
  theme_classic() + 
  theme(
    text = element_text(size = 12),
    axis.text = element_text(color = "black"),
    legend.position = "none",plot.title = element_text(hjust=0.5)
  )
ggsave("~/figures/circadian/FIGS4D.pdf",width=3,height=4)

# FIG.S4E amplitude of core circadian genes
plt.list<-list()
for(gene in c("BMAL1","PER1","PER2","CRY1","NR1D2","DBP")){
  plt.list[[gene]]=ggplot(meta2d.result[meta2d.result$CycID==gene,], aes(x = age_group, y = fold_peak2trough)) +
    geom_boxplot(aes(fill = age_group), width = 0.6, alpha = 0.8, outliers = F) +
    stat_compare_means(
      method = "wilcox.test", # or "wilcox.test" for non-parametric
      comparisons = list(c("young", "aged")), # adjust based on your groups
      label = "p.signif", # shows asterisks
      tip.length = 0.01,
      vjust = 0.1) + scale_fill_manual(values = c("#1F77B4", "#FF7F0E"))+scale_y_continuous(expand = c(0,0.5))+
    labs(x = "",y = "",title = gene) +
    theme_classic() +
    theme(
      text = element_text(size = 11),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "none"
    )
}
ggarrange(plotlist = plt.list)
ggsave("~/figures/circadian/FIGS4E.pdf",width=5,height=4)

# FIG.S4F compare number of circadian genes
plotdata<-meta2d.result %>% group_by(individual,celltype) %>% summarise(n.genes=n())
plotdata$age<-AGES[plotdata$individual]
plotdata$age_group<-ifelse(plotdata$age<40,"young","aged")

plt.list<-list()
for(this.celltype in CELL_TYPES_METACELL[!CELL_TYPES_METACELL%in%c("B naive","Plasmablast","MAIT","B memory")]){
  plt.list[[this.celltype]]=ggplot(plotdata[plotdata$celltype==this.celltype,], aes(x = age_group, y = n.genes)) +
    geom_boxplot(aes(fill = age_group), width = 0.6, alpha = 0.8, outliers = F) +
    stat_compare_means(
      method = "wilcox.test", # or "wilcox.test" for non-parametric
      comparisons = list(c("young", "aged")), # adjust based on your groups
      label = "p.signif", # shows asterisks
      tip.length = 0.01,
      vjust = -1) + scale_fill_manual(values = c("#1F77B4", "#FF7F0E"))+scale_y_continuous(expand = c(0.5,1))+
    labs(x = "",y = "",title = this.celltype) +
    theme_classic() +
    theme(
      text = element_text(size = 11),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "none"
    )
}
ggarrange(plotlist = plt.list)
ggsave("~/figures/circadian/FIGS4F.pdf",width=7,height=6)

# FIG.5A cellmarkers and some adhesin in RNA
celltypes<-names(CELL_TYPES)
meta2d.result.RNA0<-readRDS('~/analysis/core_data/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
meta2d.result.RNA<-filterMeta2dTable(meta2d.result.RNA0,float.pvalue.s1 = 0.05,float.pvalue.s2 = 0.05,int.min.individual = 1,mode = "JTK")
meta2d.result.cytof0<-readRDS('~/analysis/core_data/CYTOF.JTK.newresult.filtered.addmedianexp.addp2t.bytype.rds')
meta2d.result.cytof0$celltype<-meta2d.result.cytof0$celltype %>% gsub("_"," ",.)
meta2d.result.cytof0$CycID<-CYTOF2RNA_FEATURES[meta2d.result.cytof0$CycID]
celltypes<-CELL_TYPES_METACELL[CELL_TYPES_METACELL %in% unique(meta2d.result.cytof0$celltype)]
adhesionAndChemokine<-meta2d.result.RNA[meta2d.result.RNA$CycID %in% meta2d.result.cytof0$CycID,]
adhesionAndChemokine<-adhesionAndChemokine %>% group_by(CycID,celltype) %>% mutate(n.individual=count(ADJ.P<0.05))
adhesionAndChemokine$celltype<-gsub("_"," ",adhesionAndChemokine$celltype)
adhesionAndChemokine<-adhesionAndChemokine[c("CycID","celltype","n.individual")]
adhesionAndChemokine<-as.data.frame(adhesionAndChemokine)
shared_features<-intersect(unique(meta2d.result.RNA0$CycID),unique(meta2d.result.cytof0$CycID))
for (this.feature in shared_features) {
  for (this.celltype in unique(adhesionAndChemokine$celltype)) {
    if(nrow(adhesionAndChemokine[adhesionAndChemokine$CycID==this.feature&adhesionAndChemokine$celltype==this.celltype,])==0){
      adhesionAndChemokine=rbind(adhesionAndChemokine,data.frame("CycID"=this.feature,"celltype"=this.celltype,"n.individual"=0))
    }
  }
}

ggplot(adhesionAndChemokine)+geom_tile(aes(x=CycID,y=celltype,fill=n.individual),color="black",linewidth=0.5)+
  scale_fill_gradient2(low = "white",high = "red",midpoint = 1)+
  ggtitle("mRNA level")+
  theme(axis.text.x=element_text(hjust=1,angle=60,color="black",size=10),
        axis.text.y=element_text(color="black",size=10),
        plot.margin = margin(0,0,0,0,"cm"))+
  scale_y_discrete(limits=rev(celltypes))+ylab("")+xlab("")
ggsave("~/figures/circadian/FIG5A.pdf",width=9,height=4)

# FIG.5B cellmarkers and some adhesin in CYTOF
celltypes<-CELL_TYPES_METACELL[CELL_TYPES_METACELL %in% unique(meta2d.result.cytof0$celltype)]
#adhesionAndChemokine<-meta2d.result0[meta2d.result0$CycID %in% adhAchemo,]
meta2d.result.cytof<-meta2d.result.cytof0 %>% group_by(CycID,celltype) %>% mutate(n.individual=count(ADJ.P<0.05&fold_peak2trough>1.5))

ggplot(meta2d.result.cytof[!meta2d.result.cytof$CycID%in%c("CD15","CD45RA"),])+
  geom_tile(aes(x=CycID,y=celltype,fill=n.individual),color="black",linewidth=0.5)+
  scale_fill_gradient2(low = "white",high = "red",midpoint = 1)+
  ggtitle("protein level")+
  theme(axis.text.x=element_text(hjust=1,angle=60,color="black",size=10),
        axis.text.y=element_text(color="black",size=10),
        plot.margin = margin(0,0,0,0,"cm"))+
  scale_y_discrete(limits=rev(celltypes))+ylab("")+xlab("")
ggsave("~/figures/circadian/FIG5B.pdf",width=9,height=4)

# FIG.5C
meta2d.result.RNA<-filterMeta2dTable(meta2d.result.RNA0,float.pvalue.s1 = 0.05,float.pvalue.s2 = 0.05,mode = "JTK",int.min.individual = 1)
meta2d.result.cytof<-dplyr::filter(meta2d.result.cytof0,ADJ.P<0.05&PER>=20&PER<=28)
plotPhaseBetweenRNAandProtein(meta2d.result.RNA,meta2d.result.cytof,by.individual = T,feature = "CXCR4")
ggsave("~/figures/circadian/FIG5C.pdf",width=4,height=3)


# FIG.S5A
CBC<-read.delim("~/analysis/core_data/circadian_CBC.csv",header = T,sep = ";")
CBC<-gather(CBC,key="category",value="count",-individual,-CT)
CBC<-CBC %>% group_by(individual,category) %>% mutate(z_score=(count-mean(count))/sd(count))
CBC<-CBC %>% group_by(category,CT) %>% mutate(median_z_score=median(z_score),
                                              Q1 = quantile(z_score, 0.25, na.rm = TRUE),
                                              Q3 = quantile(z_score, 0.75, na.rm = TRUE))
CBC$category<-gsub("_lymphocyte_absolute","",CBC$category) %>% gsub("_absolute","",.)
ggplot(CBC[CBC$category%in%c("B","T","NK","Monocyte"),])+
  geom_errorbar(aes(x=CT,y=z_score,ymin=Q1,ymax=Q3),width=0.2)+facet_wrap(~category)+
  geom_line(aes(x=CT,y=median_z_score,group=category))+
  geom_point(aes(x=CT,y=median_z_score))+theme_bw()+
  theme(axis.text.x = element_text(angle=60,hjust=1))+xlab("")
ggsave("~/figures/circadian/FIGS5A.pdf",width=3,height=3)

# FIG.6A
meta2d.result0<-readRDS('~/analysis/core_data/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
plot1<-plotTopOscillatingGenesByCellType(meta2d.result0,"CD14 Mono",char.mode = "JTK")
plot2<-plotTopOscillatingGenesByCellType(meta2d.result0,"CD4 TCM",char.mode = "JTK")
plot1/plot2
ggsave('~/figures/circadian/Fig6A.pdf',width=15,height=8)

# FIG.S6A
plot1<-plotTopOscillatingGenesByCellType(meta2d.result0,"cDC2",char.mode = "JTK")
plot2<-plotTopOscillatingGenesByCellType(meta2d.result0,"CD4 Naive",char.mode = "JTK")
plot3<-plotTopOscillatingGenesByCellType(meta2d.result0,"CD8 Naive",char.mode = "JTK")
plot4<-plotTopOscillatingGenesByCellType(meta2d.result0,"CD8 TEM",char.mode = "JTK")
plot5<-plotTopOscillatingGenesByCellType(meta2d.result0,"NK",char.mode = "JTK")
plot1/plot2/plot3/plot4/plot5
ggsave('~/figures/circadian/FigS6A.pdf',width=15,height=20)


