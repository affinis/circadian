source('~/script/Rscript/circadian_core.R')


velosity_raw<-readRDS('../analysis/all.samples.velosity.rds')
velosity_merged<-mergeSplicedData(velosity_raw)
velosity_merged$spliced %>% ncol()
saveRDS(velosity_merged,"../analysis/all.samples.velosity.merged.rds")
velosity_merged<-readRDS("../analysis/all.samples.velosity.merged.rds")

test<-readRDS('../analysis/all.samples.oddsremoved.tsneran.Azimuth.rds')
cells_used<-Cells(test)
features_used<-Features(test)
shared_cells<-intersect(velosity_merged$spliced %>% colnames(),cells_used)
shared_features<-intersect(velosity_merged$spliced %>% rownames(),features_used)
velosity_merged$spliced<-velosity_merged$spliced[shared_features,shared_cells]
velosity_merged$unspliced<-velosity_merged$unspliced[shared_features,shared_cells]
velosity_merged$ambiguous<-velosity_merged$ambiguous[shared_features,shared_cells]

test@assays$spliced<-CreateAssay5Object(velosity_merged$spliced)
test@assays$unspliced<-CreateAssay5Object(velosity_merged$unspliced)
test@assays$ambiguous<-CreateAssay5Object(velosity_merged$ambiguous)

saveRDS(test,"../analysis/all.samples.oddsremoved.tsneran.Azimuth.velosity.rds")

test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.Azimuth.velosity.rds")
test$type<-test$predicted.celltype.l1.5
DefaultAssay(test)<-"spliced"
DefaultAssay(test)<-"unspliced"
DefaultAssay(test)<-"ambiguous"
DefaultAssay(test)<-"RNA"
plotPseudobulkByCT(test,"CD16 Mono")
plotPseudobulkByCTByIndividual(test,cell.type = "CD14 Mono",features = "ARNTL",split.by = "feature")


celltypes<-test$predicted.celltype.l1.5 %>% unique() %>% sort()
JTK_CYCLEouts<-getAllJTK_CYCLEouts("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5/",celltypes)
JTK_CYCLEouts<-dplyr::filter(JTK_CYCLEouts,ADJ.P<0.01,PER==24)

JTK_CYCL_bulk<-read.delim('/lustre/home/acct-medll/medll/data/analysis/tauFisher_test/GSE113883_hs_full_adj_meta2d/meta2d_GSE113883_hs_full_adj.txt')
JTK_CYCL_bulk<-dplyr::filter(JTK_CYCL_bulk,JTK_pvalue<0.01)

cycBylag<-(JTK_CYCLEouts[JTK_CYCLEouts$celltype=="NK",] %>% arrange(.,by=LAG))$CycID
cycBylag<-cycBylag[cycBylag %in% Features(test)]
cycBylag<-cycBylag[cycBylag %in% JTK_CYCL_bulk$CycID]


test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.Azimuth.velosity.rds")
test$type<-test$predicted.celltype.l1.5
celltypes<-unique(test$predicted.celltype.l1.5)
spliced_osi<-getAllJTK_CYCLEouts('../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5_spliced',celltypes = celltypes)
spliced_osi<-dplyr::filter(spliced_osi,ADJ.P<0.01)

unspliced_osi<-getAllJTK_CYCLEouts('../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5_unspliced',celltypes = celltypes)
unspliced_osi<-dplyr::filter(unspliced_osi,ADJ.P<0.01)

total_osi<-getAllJTK_CYCLEouts('../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5',celltypes = celltypes)

shared_spliceunsplice<-intersect(unique(spliced_osi$CycID),unique(unspliced_osi$CycID))
spliced_osi.f<-spliced_osi[spliced_osi$CycID %in% shared_spliceunsplice,c("CycID","LAG","celltype")]
unspliced_osi.f<-unspliced_osi[unspliced_osi$CycID %in% shared_spliceunsplice,c("CycID","LAG","celltype")]
spliced_osi.f$stat<-"spliced"
unspliced_osi.f$stat<-"unspliced"
velosity.f<-rbind(spliced_osi.f,unspliced_osi.f)
ggplot(velosity.f)+geom_boxplot(aes(x=stat,y=LAG,fill=stat))+
  stat_compare_means(aes(x=stat,y=LAG,group=stat))+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.background=element_rect(fill = "white"),
        axis.line = element_line())+
  xlab("")+ylab("peaking time (CT)")

ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig4C.pdf',width=3.5,height=3)

colnames(spliced_osi.f)[2]<-"LAG_spliced"
colnames(unspliced_osi.f)[2]<-"LAG_unspliced"
osi.f<-left_join(spliced_osi.f[1:3],unspliced_osi.f[1:3],by=c("CycID","celltype"))
osi.f<-osi.f[!is.na(osi.f$LAG_unspliced)&!is.na(osi.f$LAG_spliced),]
osi.f$LAG_diff<-abs(osi.f$LAG_spliced-osi.f$LAG_unspliced)
osi.f<-dplyr::filter(osi.f,LAG_diff<=4)

celltype.used="CD8 T"
celltype.used="NK"
cycBylag<-intersect(spliced_osi.f[spliced_osi.f$celltype==celltype.used,"CycID"],unspliced_osi.f[unspliced_osi.f$celltype==celltype.used,"CycID"]) %>% unique()
cycBylag1<-unspliced_osi.f[unspliced_osi.f$CycID %in% cycBylag & unspliced_osi.f$celltype==celltype.used,c("CycID","LAG_unspliced")]
colnames(cycBylag1)<-c("CycID", "LAG.unspliced")
cycBylag2<-spliced_osi.f[spliced_osi.f$CycID %in% cycBylag & spliced_osi.f$celltype==celltype.used,c("CycID","LAG_spliced")]
colnames(cycBylag2)<-c("CycID", "LAG.spliced")
cycBylag<-left_join(cycBylag1,cycBylag2)
cycBylag$lag.diff<-abs(cycBylag$LAG.unspliced-cycBylag$LAG.spliced)
cycBylag$mean.lag<-(cycBylag$LAG.unspliced+cycBylag$LAG.spliced)/2
cycBylag<-cycBylag[cycBylag$lag.diff<=4|cycBylag$lag.diff>=22,]
cycBylag<-(cycBylag %>% arrange(.,by=mean.lag))$CycID

DefaultAssay(test)<-"spliced"
data_spliced<-plotGeneExpressionHeatmapByCT(test,celltype.used,gene.list=cycBylag,return.data = T)
data_spliced$relative_expression<-as.vector(data_spliced$relative_expression)
data_spliced$CT<-paste0("spliced",data_spliced$CT)

DefaultAssay(test)<-"unspliced"
data_unspliced<-plotGeneExpressionHeatmapByCT(test,celltype.used,gene.list=cycBylag,return.data = T)
data_unspliced$relative_expression<-as.vector(data_unspliced$relative_expression)
data_unspliced$CT<-paste0("unspliced",data_unspliced$CT)

data<-rbind(data_spliced,data_unspliced)
ggplot(data)+geom_bin_2d(aes(x=genes,y=CT,fill=relative_expression))+
  scale_fill_gradientn(colors = rev(brewer.pal(11,"RdBu")))+
  scale_x_discrete(limits=cycBylag)+
  scale_y_discrete(limits=rev(c("unsplicedCT1","splicedCT1","unsplicedCT5","splicedCT5",
                                "unsplicedCT9","splicedCT9","unsplicedCT13","splicedCT13",
                                "unsplicedCT17","splicedCT17","unsplicedCT21","splicedCT21")))+
  ylab("")+xlab("")+theme(axis.text.x=element_blank(),axis.ticks=element_blank(),legend.position = "top")+NoLegend()+
  ggtitle(celltype.used)

ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig4B1.pdf',width=8,height=3)
ggsave('/lustre/home/acct-medll/medll/figures/circadian/Fig4B2.pdf',width=5,height=3)

internal_controls<-(total_osi %>% dplyr::filter(.,ADJ.P>=0.5,AMP>=50) %>% arrange(desc(AMP)))$CycID %>% table() %>% sort() %>% tail(50)
internal_controls<-names(internal_controls)
saveRDS(internal_controls,'../analysis/internal_controls.rds')

