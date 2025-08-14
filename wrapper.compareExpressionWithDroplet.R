source('/tmpdata/LyuLin/script/circadian/circadian_core.R')

# upstream: wrapper.testRhythmicity.R, metacell2srt
# downstream: wrapper.addcor2batch.R, wrapper.addpeak2trough.R
# dependency: circadian_core.R
# caller: NSF
##
# speed: 36000 rows elapsed about 1 hours

droplets<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/seacell.droplet.0.08.rds')
metacells<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/seacell.0.08.bytype.rds')

#Allresult<-readJTKFromMetaCells('/tmpdata/LyuLin/analysis/circadian/R/seacell.meta2d.0.04.12Healthy.bytype/')
Allresult<-readLSFromMetaCells('/tmpdata/LyuLin/analysis/circadian/R/seacell.meta2d.0.08.bytype/')
Allresult.filtered<-dplyr::filter(Allresult,p<0.05,BH.Q<0.05,Period>=20,Period<=28,!is.na(Period))
saveRDS(Allresult.filtered,'/tmpdata/LyuLin/analysis/circadian/R/LS.result.filtered.bytype.0.08.rds')

droplet.mat<-LayerData(droplets,layer="count")
metacell.mat<-LayerData(metacells,layer="count")
batches<-BATCH %>% as.vector() %>% unique() %>% sort()
metacell.meta<-metacells@meta.data
droplet.meta<-droplets@meta.data

droplet.names<-colnames(droplet.mat)

#AllJTKresult.filtered$p_to_background<-NA
for(i in 1:nrow(Allresult.filtered)){
  feature=Allresult.filtered$CycID[i]
  celltype=Allresult.filtered$celltype[i]
  individual=Allresult.filtered$individual[i]
  #message(paste(feature,celltype,individual,sep = ":"))
  cells_to_compare=rownames(metacell.meta[metacell.meta$type==celltype&metacell.meta$individual==individual,])
  cells_to_compare=grep("TF_",cells_to_compare,value = T)
  observation=metacell.mat[feature,cells_to_compare] %>% as.vector()
  meta.cell.size.observation=metacell.meta[cells_to_compare,"meta.cell.size"]
  observation=observation/meta.cell.size.observation
  if(feature %in% rownames(droplet.mat)){
    background.cells=sample(droplet.names,length(observation),replace = F)
    background=droplet.mat[feature,background.cells] %>% as.vector()
    meta.cell.size.background=droplet.meta[background.cells,"meta.cell.size"]
    background=background/meta.cell.size.background
  }else{
    background=rep(0,length(observation))
  }
  #print(head(observation))
  #print(head(background))
  observation=na.omit(observation)
  background=na.omit(background)
  #if(length(observation<=2)){
  #  next
  #}
  if(i%%10==0){
    print(paste0(i,"/",nrow(Allresult.filtered)))
  }
  #AllJTKresult.filtered$p_to_background[i]=wilcox.test(observation, background, paired = FALSE)$p.value
  Allresult.filtered$fold_to_background[i]=mean(observation)/mean(background)
}

#AllJTKresult.filtered$padj_to_background<-p.adjust(AllJTKresult.filtered$p_to_background, method = "BH")
saveRDS(Allresult.filtered,"/tmpdata/LyuLin/analysis/circadian/R/LS.result.filtered.addp2bkg.bytype.0.08.rds")
