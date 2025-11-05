source('~/script/circadian/circadian_core.R')
integrateSubset<-function(subsrt,method=CCAIntegration,classic.normalize=T,k.weight=100){
  subsrt[["RNA"]]=split(subsrt[["RNA"]], f = subsrt$sample,layers = c("counts"))
  if(classic.normalize){
    subsrt=NormalizeData(subsrt)
    subsrt=FindVariableFeatures(subsrt)
    subsrt=ScaleData(subsrt)
  }else{
    subsrt=SCTransform(subsrt)
  }
  subsrt=RunPCA(subsrt)
  assign("runtime.integrateDubset",subsrt,envir = .GlobalEnv)
  if(classic.normalize){
    subsrt=IntegrateLayers(object=subsrt, method=method, orig.reduction="pca",
                           new.reduction="integrated.cca",verbose = FALSE,k.weight=k.weight)
    subsrt[["RNA"]]=JoinLayers(subsrt[["RNA"]])
  }else{
    subsrt=IntegrateLayers(object=subsrt, method=method, orig.reduction="pca",
                           normalization.method = "SCT",k.weight=k.weight,
                           new.reduction="integrated.cca",verbose = FALSE)
    subsrt[["RNA"]]=JoinLayers(subsrt[["RNA"]])
  }
  subsrt=FindNeighbors(subsrt, reduction = "integrated.cca")
  subsrt=FindClusters(subsrt, resolution = 0.2)
  assign("runtime.integrateDubset",subsrt,envir = .GlobalEnv)
  subsrt=RunUMAP(subsrt, reduction = "integrated.cca",dims=1:50)
  return(subsrt)
}


srts<-readRDS('~/analysis/circadian/R/16individual.srt.annotated.rds')
srts$sample<-paste0(srts$individual,"_",srts$CT)
srts<-subset(srts,cells=sample(Cells(srts),size = 100000,replace = F))
srts<-integrateSubset(srts,classic.normalize = F,k.weight = 50)

saveRDS(srts,'~/analysis/circadian/R/16individual.srt.annotated.integrated.rds')
