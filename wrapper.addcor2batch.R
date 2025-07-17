source('/tmpdata/LyuLin/script/circadian/circadian_core.R')

# upstream: wrapper.compareExpressionWithDroplet.R
# downstream: 
# dependency: circadian_core.R
# caller: NSF
##
# speed: 36000 rows elapsed about 3.5 hours

calculateCorBatch<-function(arg.type,arg.feature,mat){
  this.data=dplyr::filter(mat,feature==arg.feature,type==arg.type)[c("batch","individual","count")] %>% spread(key="individual",value="count")
  this.data=column_to_rownames(this.data,"batch")
  #print(cor(this.data))
  this.cor=cor(this.data)[upper.tri(cor(this.data), diag = FALSE)]
  this.batch.data=data.frame("feature"=rep(arg.feature,length(this.cor)),"type"=rep(arg.type,length(this.cor)),"cor"=this.cor)
  return(this.batch.data)
}
AllJTKresult.filtered<-readRDS("/tmpdata/LyuLin/analysis/circadian/R/JTK.result.filtered.addp2bkg.bytype.0.04.rds")

matBySample.nr<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/matBySample.rds')

AllJTKresult.filtered$cor_to_batch=NA
for (i in 1:nrow(AllJTKresult.filtered)) {
  if(is.na(AllJTKresult.filtered$cor_to_batch[i])){
    feature=AllJTKresult.filtered$CycID[i]
    type=AllJTKresult.filtered$celltype[i] %>% gsub(" ","_",.)
    #message(paste0(i," ",feature," ",type))
    cor=median(calculateCorBatch(type,feature,matBySample.nr)$cor)
    AllJTKresult.filtered$cor_to_batch[i]=cor
  }else{
    next
  }
  if(i%%100==0){
    message(i)
  }
}

saveRDS(AllJTKresult.filtered,"/tmpdata/LyuLin/analysis/circadian/R/JTK.result.filtered.addp2bkg.addcor2batch.bytype.0.04.rds")
