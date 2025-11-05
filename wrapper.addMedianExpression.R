source('~/script/circadian/circadian_core.R')

# upstream: wrapper.testRhythmicity.R, metacell2srt
# downstream: NSF
# dependency: circadian_core.R
# caller: NSF
##
# speed: 15min for ~40000 cells x ~40000 features

### read srt object
## metacell data
#metacells<-readRDS('~/analysis/circadian/R/seacell.16individual.rds')

## CYTOF data
#metacells<-readRDS('~/analysis/circadian/R/CYTOF.300k.reannotation.srt.rds')
#metacells$type<-metacells$manual.new

## velocity data
metacells<-readRDS('~/analysis/circadian/R/16individual.unspliced.srt.rds')

### read JTK result
#JTKresult<-readRDS("~/analysis/circadian/R/JTK.result.filtered.16individual.addp2t.addfold2bkg.bytype.rds")
#JTKresult<-readRDS("~/analysis/circadian/R/CYTOF.JTK.newresult.filtered.rds")
JTKresult<-readRDS("~/analysis/circadian/R/JTK.unspliced.filtered.16individual.rds")

#JTKresult$celltype<-JTKresult$celltype %>% gsub(" ","_",.)

metacells<-NormalizeData(metacells,normalization.method = "RC",scale.factor = 1000000)

metacell.mat<-LayerData(metacells,layer="data")
metacell.mat<-rownames_to_column(as.data.frame(metacell.mat),var="CycID")
metacell.mat<-gather(metacell.mat,key="observation",value="value",-CycID)

metadata<-metacells@meta.data[,c("type","individual")]
metadata<-rownames_to_column(as.data.frame(metadata),var="observation")
metacell.mat<-left_join(metacell.mat,metadata,by="observation")
metacell.mat<-metacell.mat %>% group_by(CycID,type,individual) %>% summarise(median_expression=median(value))
colnames(metacell.mat)[2]<-"celltype"
#metacell.mat$celltype<-gsub("_"," ",metacell.mat$celltype)

JTKresult<-left_join(JTKresult,metacell.mat,by=c("CycID","celltype","individual"))

### save result
#saveRDS(JTKresult,'~/analysis/circadian/R/CYTOF.JTK.newresult.filtered.addmedianexp.bytype.rds')
saveRDS(JTKresult,'~/analysis/circadian/R/JTK.unspliced.filtered.16individual.addmedianexp.rds')
