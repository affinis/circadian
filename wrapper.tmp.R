source('/tmpdata/LyuLin/script/circadian/circadian_core.R')

calculateCorBatch<-function(arg.type,arg.feature,mat){
  this.data=dplyr::filter(mat,feature==arg.feature,type==arg.type)[c("batch","individual","count")] %>% spread(key="individual",value="count")
  this.data=column_to_rownames(this.data,"batch")
  #print(cor(this.data))
  this.cor=cor(this.data)[upper.tri(cor(this.data), diag = FALSE)]
  this.batch.data=data.frame("feature"=rep(arg.feature,length(this.cor)),"type"=rep(arg.type,length(this.cor)),"cor"=this.cor)
  return(this.batch.data)
}

AllJTKresult<-readJTKFromMetaCells('/tmpdata/LyuLin/analysis/circadian/R/seacell.meta2d.0.04')
AllJTKresult.filtered<-filter(AllJTKresult,ADJ.P<0.05,PER>=20,PER<=28)
#filter for genes oscillate in at least 2 people
filtered.genes<-AllJTKresult.filtered[c("CycID","celltype")] %>% table() %>% as.data.frame() %>% dplyr::filter(.,Freq>=2)
filtered.genes<-filtered.genes$CycID %>% unique()
AllJTKresult.filtered<-AllJTKresult.filtered[AllJTKresult.filtered$CycID %in% filtered.genes,]

matBySample.nr<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/matBySample.rds')

AllJTKresult.filtered$cor_to_batch=NA
for (i in 1:nrow(AllJTKresult.filtered)) {
  if(is.na(AllJTKresult.filtered$cor_to_batch[i])){
    feature=AllJTKresult.filtered$CycID[i]
    type=AllJTKresult.filtered$celltype[i] %>% gsub(" ","_",.)
    message(paste0(i," ",feature," ",type))
    cor=median(calculateCorBatch(type,feature,matBySample.nr)$cor)
    AllJTKresult.filtered$cor[i]=cor
  }else{
    next
  }
  if(i%%100==0){
    message(i)
  }
}

saveRDS(AllJTKresult.filtered,"/tmpdata/LyuLin/analysis/circadian/R/JTK.result.filtered.rds")
