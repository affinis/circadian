source('~/script/Rscript/circadian_core.R')

#FIG. 1B
test<-readRDS("/tmpdata/LyuLin/analysis/circadian/R/7mixed.integrated.annotated.clean.sct.Azimuth.rds")
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

#FIG. 1C
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

#FIG. 1D
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
  scale_fill_gradient(low = "darkblue", high = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black"),
        axis.text.y = element_text(color="black")) +
  labs(x = "manual annotations", y = "Azimuth annotations")

# FIG. 1E


#FIG. 2A
srt.metacell<-readRDS("/tmpdata/LyuLin/analysis/circadian/R/seacell.0.08.bytype.rds")
VlnPlot(srt.metacell,CIRCADIAN_GENES_MAIN,group.by = "type",stack = T,flip = T,cols = generateColor(15,seed = 2025))+NoLegend()+xlab("")+theme(strip.text=element_text(size=10,face="plain"))

test<-readRDS("/tmpdata/LyuLin/analysis/circadian/R/seacell.0.08.bytype.rds")
JTK.individual<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/JTK.result.filtered.addp2bkg.addcor2batch.rds')
#plotdata<-plotPseudobulkByCTByIndividual(test,"CD14 Mono","ARNTL",return.data = T)
#plotdata2<-plotPseudobulkByCTByIndividual(test,"NK","CRY1",return.data = T)
#plotdata3<-plotPseudobulkByCTByIndividual(test,"Plasma","PER3",return.data = T)
#plotdata4<-plotPseudobulkByCTByIndividual(test,"CD16 Mono","NR1D1",return.data = T)
plotMetaCellByIndividual(test,cell.type = "CD14 Mono",feature = "NR1D2",layer = "data")

ggplot(plotdata)+geom_boxplot(aes(x=time,y=relative_expression))+
  geom_point(aes(x=time,y=relative_expression),alpha=0.5,color="blue",size=3)+
  scale_x_discrete(limits=CT_TIME_ORDER)+theme_bw()+
  facet_wrap(~feature,scales="free_y")+ylab("z-score")+xlab("")+
  ggtitle("CD14 Mono\np.adj: 0.016, peak: CT18")
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig2A1.pdf',width = 3,height=3)

ggplot(plotdata2)+geom_boxplot(aes(x=time,y=relative_expression))+
  geom_point(aes(x=time,y=relative_expression),alpha=0.5,color="blue",size=3)+
  scale_x_discrete(limits=CT_TIME_ORDER)+theme_bw()+
  facet_wrap(~feature,scales="free_y")+ylab("z-score")+xlab("")+
  ggtitle("NK\np.adj: 0.001, peak: CT10")
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig2A2.pdf',width = 3,height=3)

ggplot(plotdata3)+geom_boxplot(aes(x=time,y=relative_expression))+
  geom_point(aes(x=time,y=relative_expression),alpha=0.5,color="blue",size=3)+
  scale_x_discrete(limits=CT_TIME_ORDER)+theme_bw()+
  facet_wrap(~feature,scales="free_y")+ylab("z-score")+xlab("")+
  ggtitle("Plasma\np.adj: 0.037, peak: CT6")
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig2A3.pdf',width = 3,height=3)

ggplot(plotdata4)+geom_boxplot(aes(x=time,y=relative_expression))+
  geom_point(aes(x=time,y=relative_expression),alpha=0.5,color="blue",size=3)+
  scale_x_discrete(limits=CT_TIME_ORDER)+theme_bw()+
  facet_wrap(~feature,scales="free_y")+ylab("z-score")+xlab("")+
  ggtitle("CD16 Mono\np.adj: 0.032, peak: CT4")
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig2A4.pdf',width = 3,height=3)

plotdataFIG2B1<-test@meta.data[,c("nFeature_RNA","type")] %>% group_by(type) %>% summarise(median_genes=max(nFeature_RNA))
plotdataFIG2B1<-left_join(plotdataFIG2B1,(JTK.individual %>% dplyr::filter(.,ADJ.P<=0.01))$celltype %>% table() %>% as.data.frame() %>% `colnames<-`(.,c("type","rhythmic")))
plotdataFIG2B1$percent_rhythmic<-plotdataFIG2B1$rhythmic*100/plotdataFIG2B1$median_genes
plotdataFIG2B1$percent_rhythmic<-round(plotdataFIG2B1$percent_rhythmic,2)
plotdataFIG2B1[is.na(plotdataFIG2B1)]<-0

FIG2B1<-ggplot(plotdataFIG2B1)+geom_bar(aes(y=type,x=percent_rhythmic),stat="identity",fill="black",alpha=0.5,color="black")+
  geom_text(data=plotdataFIG2B1[plotdataFIG2B1$percent_rhythmic>=2,],aes(y=type,x=percent_rhythmic/2,label=percent_rhythmic),color="white")+
  geom_text(data=plotdataFIG2B1[plotdataFIG2B1$percent_rhythmic<2,],aes(y=type,x=percent_rhythmic+1,label=percent_rhythmic),color="black")+
  scale_y_discrete(limits=rev(celltypes %>% sort()))+scale_x_continuous(expand=c(0,0))+
  theme(panel.grid = element_blank(),panel.background=element_rect(fill="white"))+
  ylab("")+xlab("percentage of\nrhythmic genes(%)")

plotdataFIG2B<-JTK.individual[JTK.individual$CycID %in% CIRCADIAN_GENES_MAIN,]
plotdataFIG2B$celltype<-factor(plotdataFIG2B$celltype,levels=rev(celltypes %>% sort()))
FIG2B2<-ggplot(plotdataFIG2B)+geom_point(aes(x=celltype,y=CycID,size=-log10(ADJ.P),color=log2(AMP)))+
  scale_color_gradientn(colours = c("#507AAF","#BE5C37"))+
  theme(axis.text.x = element_text(angle=45,hjust=1),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        panel.background = element_rect(fill="white"))+coord_flip()+
  xlab("")+ylab("rhythmicity of\ncore circadian genes")

ggarrange(FIG2B1,FIG2B2,ncol=2,widths=c(1.5,2.5),align="h")
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig2B.pdf',width=6.5,height=4.5)

#FIG. 3A
plotdataFIG3A1<-test$predicted.celltype.l1.5 %>% table() %>% as.data.frame()
test$type<-test$predicted.celltype.l1.5
celltypes<-sort(unique(test$predicted.celltype.l1.5))
colnames(plotdataFIG3A1)<-c("cell_type","number_of_cell")
FIG3A1<-ggplot(plotdataFIG3A1)+geom_bar(aes(y=cell_type,x=log10(number_of_cell)),stat="identity",fill="blue",alpha=0.5,color="black")+
  geom_text(aes(y=cell_type,x=log10(number_of_cell)/2,label=number_of_cell),color="white")+
  scale_y_discrete(limits=rev(plotdataFIG3A1$cell_type))+scale_x_continuous(expand=c(0,0))+
  theme(panel.grid = element_blank(),panel.background=element_rect(fill="white"))+ylab("")+xlab("log10(# of cell)")

JTK_CYCL<-NULL
for(block_data in list.files("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5/")){
  block_path=paste0("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5/",block_data,"/")
  block=getJTK_CYCLEouts(block_path,celltypes=celltypes,filter.data = F)
  if(is.null(JTK_CYCL)){
    JTK_CYCL=block
  }else{
    JTK_CYCL=rbind(JTK_CYCL,block)
  }
}

JTK_CYCL_bulk<-read.delim('/lustre/home/acct-medll/medll/data/analysis/tauFisher_test/GSE113883_hs_full_adj_meta2d/meta2d_GSE113883_hs_full_adj.txt')
JTK_CYCL_bulk<-dplyr::filter(JTK_CYCL_bulk,JTK_pvalue<0.01)
plotdataFIG3A2<-test@meta.data[,c("nFeature_RNA","type")] %>% group_by(type) %>% summarise(median_genes=max(nFeature_RNA))
plotdataFIG3A2<-left_join(plotdataFIG3A2,(JTK_CYCL %>% dplyr::filter(.,ADJ.P<=0.01,CycID %in% JTK_CYCL_bulk$CycID))$celltype %>% table() %>% as.data.frame() %>% `colnames<-`(.,c("type","rhythmic")))
plotdataFIG3A2$percent_rhythmic<-plotdataFIG3A2$rhythmic*100/plotdataFIG3A2$median_genes
plotdataFIG3A2$percent_rhythmic<-round(plotdataFIG3A2$percent_rhythmic,2)
plotdataFIG3A2[is.na(plotdataFIG3A2)]<-0

FIG3A2<-ggplot(plotdataFIG3A2)+geom_bar(aes(y=type,x=percent_rhythmic),stat="identity",fill="red",alpha=0.5,color="black")+
  geom_text(data=plotdataFIG3A2[plotdataFIG3A2$percent_rhythmic>=2,],aes(y=type,x=percent_rhythmic/2,label=percent_rhythmic),color="white")+
  geom_text(data=plotdataFIG3A2[plotdataFIG3A2$percent_rhythmic<2,],aes(y=type,x=percent_rhythmic+1,label=percent_rhythmic),color="black")+
  scale_y_discrete(limits=rev(plotdataFIG3A1$cell_type))+scale_x_continuous(expand=c(0,0))+
  theme(panel.grid = element_blank(),panel.background=element_rect(fill="white"),
        axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  ylab("")+xlab("percentage of\nrhythmic genes(%)")

plotdataFIG3A3<-JTK_CYCL[JTK_CYCL$CycID %in% CIRCADIAN_GENES_MAIN,]
plotdataFIG3A3$celltype<-factor(plotdataFIG3A3$celltype,levels=rev(celltypes))
FIG3A3<-ggplot(plotdataFIG3A3)+geom_point(aes(x=celltype,y=CycID,size=-log10(ADJ.P),color=log2(AMP)))+
  scale_color_gradientn(colours = c("#507AAF","#BE5C37"))+
  theme(axis.text.x = element_text(angle=45,hjust=1),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        panel.background = element_rect(fill="white"))+coord_flip()+
  xlab("")+ylab("rhythmicity of\ncore circadian genes")

ggarrange(FIG3A1,FIG3A2,FIG3A3,ncol=3,widths=c(1.5,1.1,2.7),align="h")
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig3A.pdf',width=8,height=4.5)

#FIG. 3B
plotdataFIG3B<-NULL
for(celltype in c("CD14 Mono","CD16 Mono","CD4 T")){
  plotdataFIG3B.part1=plotPseudobulkByCT(test,celltype,features=c("ARNTL","CRY1","NR1D1","PER1","TEF"),return.data = T,rep=2,proportion=0.5)
  plotdataFIG3B.part1$type=celltype
  plotdataFIG3B.part2=plotPseudobulkByCT(test,celltype,features=c("ARNTL","CRY1","NR1D1","PER1","TEF"),return.data = T,rep=2,proportion=0.5)
  plotdataFIG3B.part2$type=celltype
  plotdataFIG3B.part2$time=case_when(plotdataFIG3B.part2$time=="CT1"~"CT25",
                                     plotdataFIG3B.part2$time=="CT5"~"CT29",
                                     plotdataFIG3B.part2$time=="CT9"~"CT33",
                                     plotdataFIG3B.part2$time=="CT13"~"CT37",
                                     plotdataFIG3B.part2$time=="CT17"~"CT41",
                                     plotdataFIG3B.part2$time=="CT21"~"CT45")
  plotdataFIG3B.part=rbind(plotdataFIG3B.part1,plotdataFIG3B.part2)
  if(is.null(plotdataFIG3B)){
    plotdataFIG3B=plotdataFIG3B.part
  }else{
    plotdataFIG3B=rbind(plotdataFIG3B,plotdataFIG3B.part)
  }
}
summ.data<-plotdataFIG3B %>% group_by(feature,type) %>% summarise(min=min(count)+1)
plotdataFIG3B[plotdataFIG3B$count!=0,]
plotdataFIG3B<-left_join(plotdataFIG3B,summ.data,by=c("feature","type"))
plotdataFIG3B$relative_expression<-plotdataFIG3B$count/plotdataFIG3B$min
plotdataFIG3B.sup<-plotdataFIG3B %>% group_by(feature,type,time) %>% summarise(median_expression=median(relative_expression))

use_color<-generateColor(30,col.dist.min=0.3,seed = 10)
ggplot(plotdataFIG3B)+geom_boxplot(aes(x=time,y=relative_expression,color=type))+
  geom_point(data=plotdataFIG3B.sup,aes(x=time,y=median_expression,color=type))+
  geom_line(data=plotdataFIG3B.sup,aes(x=time,y=median_expression,group=type,color=type))+
  scale_color_manual(values=c(use_color[4],use_color[5],use_color[6],use_color[9]))+
  scale_x_discrete(limits=c(CT_TIME_ORDER,"CT25","CT29","CT33","CT37","CT41","CT45"))+
  facet_wrap(~feature,scales = "free_y",nrow=2)+xlab("")+
  theme(axis.text.x=element_text(angle=60,hjust=1,color="black",size=9),
        axis.text.y=element_text(color="black",size=9),legend.text=element_text(color="black",size=9),
        panel.background=element_rect(fill="white"),
        legend.position=c(0.85,0.1))
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig3B.pdf',width=6,height=3)


#FIG. 3C
JTK.result.filtered<-JTK_CYCL[JTK_CYCL$ADJ.P<=0.01&JTK_CYCL$CycID %in% JTK_CYCL_bulk$CycID,]
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
ggplot(share.data[share.data$from.type!="Eryth"&share.data$to.type!="Eryth",])+geom_raster(aes(x=from.type,y=to.type,fill=shared.ratio))+
  scale_fill_gradientn(colours = c("#283168","#009fc0","#f0f0cc","white"),values = c(0,0.5,0.9,1))+
  theme(axis.text.x = element_text(angle=60,hjust=1),plot.margin=margin(0,0,0,0,"cm"))+ylab("")+xlab("")
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig3C.pdf',width=5,height=4)

#FIG. 3D
unique_in_sc<-JTK_CYCL[JTK_CYCL$ADJ.P<=0.01&!(JTK_CYCL$CycID %in% JTK_CYCL_bulk$CycID)&JTK_CYCL$PER==24,"CycID"] %>% unique()
shared_between_sc_bulk<-JTK_CYCL[JTK_CYCL$ADJ.P<=0.01&(JTK_CYCL$CycID %in% JTK_CYCL_bulk$CycID)&JTK_CYCL$PER==24,"CycID"] %>% unique()
#(LayerData(test,layer="count")[unique_in_sc,]>0) %>% rowSums() %>% median()
#(LayerData(test,layer="count")[shared_between_sc_bulk,]>0) %>% rowSums() %>% median()
NK_genes_by_lag<-(JTK.result.filtered[JTK.result.filtered$celltype=="NK",] %>% arrange(.,by=LAG))$CycID
plotGeneExpressionHeatmapByCT(test,"NK",NK_genes_by_lag)
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig3D.pdf',width=4,height=2.5)

#FIG. 3E
gdT_genes_by_lag<-(JTK.result.filtered[JTK.result.filtered$celltype=="gdT",] %>% arrange(.,by=LAG))$CycID
plotGeneExpressionHeatmapByCT(test,"gdT",gdT_genes_by_lag)
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig3E.pdf',width=4,height=2.5)

#FIG. 3F
use_color<-generateColor(30,col.dist.min=0.3,seed = 10)[-11]
ggplot(JTK.result.filtered[JTK.result.filtered$celltype!="Eryth",])+
  geom_bar(aes(x=LAG,fill=celltype),alpha=0.5,color="black")+
  facet_wrap(~celltype,scales = "free_y",ncol=6)+
  scale_fill_manual(values=use_color)+NoLegend()+
  theme(axis.text=element_text(color="black"))+ylab("count of peaking genes")+
  xlab("ciradian time")
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig3F.pdf',width=8,height=4)


#FIG. 3G
plotPhaseofCellType(JTK.result.filtered[JTK.result.filtered$celltype!="Eryth",],p.cutoff = 0.01,
                    colors=generateColor(30,col.dist.min=0.3,seed = 10)[-11])+NoLegend()
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig3G.pdf',width=4,height=4)

#FIG. 4A
celltypes<-test$predicted.celltype.l1.5 %>% sort() %>% unique()
celltypes<-celltypes[celltypes!="Eryth"]
adh<-c("SELPLG","SELL","CD44","ITGAL","ITGAM","ITGAX","ITGB2","ITGB1","ITGA2","ITGA4","ITGA5","ITGA6")
adhAchemo<-c("SELPLG","SELL","CD44","ITGAL","ITGAM","ITGAX","ITGB2","ITGB1","ITGA2","ITGA4","ITGA5","ITGA6",
             "CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CCR8","CCR9","CCR10","CXCR2","CXCR3","CXCR4","CXCR5","CXCR6","CX3CR1")
adhesionAndChemokine<-JTK_CYCL[JTK_CYCL$CycID %in% adhAchemo&JTK_CYCL$celltype!="Eryth" ,]
#adhesionAndChemokine<-JTK.individual[JTK.individual$CycID %in% adhAchemo&JTK.individual$celltype!="Eryth" ,]
adhesionAndChemokine$`p.adj<0.05`<-case_when(adhesionAndChemokine$ADJ.P<0.05 ~ "oscillating",
                                             adhesionAndChemokine$AMP<=0.1 ~ "no/low expression",
                                             TRUE~"not oscillating")

bottom<-ggplot(adhesionAndChemokine)+geom_tile(aes(x=CycID,y=celltype,fill=get("p.adj<0.05")),color="black",linewidth=0.5)+
  scale_y_discrete(limits=rev(celltypes))+scale_fill_manual(values=c("white","grey","red"))+scale_x_discrete(limits=adhAchemo)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,color="black",size=12),
        axis.text.y=element_text(color="black",size=10),
        panel.border = element_rect(fill=NA,color="black"),
        plot.margin = margin(0,0.5,0,0,"cm"))+
  xlab("")+ylab("")+guides(fill=guide_legend(title="Oscillating during 24 hours",position = "bottom"))


adhesionAndChemokine.f<-adhesionAndChemokine[adhesionAndChemokine$ADJ.P<0.05,]
mid<-ggplot(adhesionAndChemokine.f)+geom_boxplot(aes(x=CycID,y=LAG),outliers = FALSE)+
  geom_jitter(aes(x=CycID,y=LAG,color=celltype))+
  scale_x_discrete(limits=adhAchemo)+xlab("")+ylab("peaking time by\ncell type (CT)")+
  scale_color_manual(values=generateColor(30,col.dist.min=0.3,seed = 10)[-11])+
  theme(axis.text.x=element_blank(),axis.ticks = element_blank(),plot.margin = margin(0,0.5,0,0,"cm"),
        axis.text.y=element_text(color="black",size=10),panel.grid.major.y=element_line(colour = "grey"),
        panel.background = element_rect(fill=NA))+NoLegend()

adhesionAndChemokine.f$gene_type<-ifelse(adhesionAndChemokine.f$CycID %in% adh,"adhesion molecules","chemokine receptors")
top<-ggplot(adhesionAndChemokine.f)+geom_bar(aes(x=CycID,fill=gene_type),color="black")+
  scale_x_discrete(limits=adhAchemo)+xlab("")+ylab("count of\ncell types")+scale_y_continuous(n.breaks = 4)+
  theme(axis.text.x=element_blank(),axis.ticks = element_blank(),plot.margin = margin(0.5,0.5,0,0,"cm"),
        axis.text.y=element_text(color="black",size=10),panel.grid.major.y=element_line(colour = "grey"),
        panel.background = element_rect(fill=NA),legend.position="top")

ggarrange(top,mid,bottom,nrow=3,align = "v",heights = c(1,1,4),common.legend = T)
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig4A.pdf',width=6,height=7)



