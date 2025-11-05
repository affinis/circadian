source('~/script/circadian/circadian_core.R')

# merge 4 individual in batch3
i=1
srts<-list()
for (file in list.files('~/analysis/circadian/R/',pattern = "X5_")){
  this.srt=readRDS(file.path('~/analysis/circadian/R',file))
  srts[[i]]=subset(this.srt,individual!="unknown")
  i=i+1
}
srts=merge(srts[[1]],srts[2:length(srts)])
saveRDS(srts,'~/analysis/circadian/R/4multiplexedindividual.srt.raw.rds')
anno<-read.delim('~/analysis/circadian/R/cell.annotation.clean.sct.Azimuth.tsv')
rownames(anno)<-anno$cell_id
srts<-AddMetaData(srts,anno)
srts<-subset(srts,!is.na(predicted.celltype.l2))
srts$sample<-BATCH[paste0(srts$individual,"-",srts$CT)] %>% as.vector()
saveRDS(srts,'~/analysis/circadian/R/4multiplexedindividual.srt.annotated.rds')

srts<-readRDS('~/analysis/circadian/R/4multiplexedindividual.srt.annotated.rds')
srts2<-readRDS('~/analysis/circadian/R/12individual.srt.annotated.rds')
srts<-merge(srts,srts2)
srts[["RNA"]]<-JoinLayers(srts[["RNA"]])
saveRDS(srts,'~/analysis/circadian/R/16individual.srt.annotated.rds')