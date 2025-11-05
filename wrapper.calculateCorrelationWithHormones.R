source('~/script/circadian/circadian_core.R')

# This is a wrapper to find correlation of pseudobulk genes and hormones

# upstream: wrapper.testRhythmicity.R
# downstream: NSF
# dependency: NSF
# caller: NSF
##
# speed: 20-30 min

#test<-readRDS('~/analysis/circadian/R/16individual.srt.annotated.rds')
#test$type<-test$predicted.celltype.l2
#test<-NormalizeData(test,normalization.method = "RC",scale.factor = 1000000)

JTK<-readRDS('~/analysis/circadian/R/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
test.mat<-readRDS('~/analysis/circadian/R/16individual.pseudobulk.bytype.CPM.rds')

all_result=NULL
#for (this.individual in c("KD","ZYR","LYH","JJC")) {
  for (this.celltype in unique(JTK$celltype)) {
    #message(paste0(this.individual,this.celltype,collapse = "-"))
    message(paste0(this.celltype,collapse = "-"))
    #this.features=as.data.frame(dplyr::filter(JTK,individual==this.individual,celltype==this.celltype))
    this.features=as.data.frame(dplyr::filter(JTK,celltype==this.celltype))
    if(nrow(this.features)==0){
      next
    }else{
      this.features=unique(this.features$CycID)
      message(length(this.features))
      data=getExpressionZscores(S4.srt = NULL,char.cell.type=this.celltype,vector.features=this.features,vector.individual=c("KD","ZYR","LYH","JJC"),
                                bool.normalize.data=T,bool.renormalize.seurat = F,char.melatonin.file='~/analysis/circadian/R/ELISA_melatonin.rds',
                                char.cortisol.file = '~/analysis/circadian/R/ELISA_cortisol.rds',bool.relative = T,char.data.from.mat = test.mat)
      data.melatonin=data[data$features=="melatonin",c("CT","individual","values")]
      data.cortisol=data[data$features=="cortisol",c("CT","individual","values")]
      data.gene=data[!(data$features%in%c("cortisol","melatonin")),c("features","CT","individual","values")]
      
      colnames(data.melatonin)[3]="melatonin"
      colnames(data.cortisol)[3]="cortisol"
      
      data.final=left_join(data.gene,data.melatonin,by=c(c("CT","individual")))
      data.final=left_join(data.final,data.cortisol,by=c(c("CT","individual")))
      data.final$individual=NULL
      print(head(data.final))
      write.table(x = data.final,file="~/analysis/circadian/tmp/wrapper.calculateCorrelationWithHormones.tmp",quote = F,sep = '\t')
      data.result=data.final %>% group_by(features) %>% summarise(cor2melatonin=cor.test(values,melatonin,method = "spearman")$estimate,p.value.melatonin=cor.test(values,melatonin,method = "spearman")$p.value,
                                                                             cor2cortisol=cor.test(values,cortisol,method = "spearman")$estimate,p.value.cortisol=cor.test(values,cortisol,method = "spearman")$p.value)
      data.result$celltype=this.celltype
      print(head(data.result))
      if(is.null(all_result)){
        all_result=data.result
      }else{
        all_result=rbind(all_result,data.result)
      }
    }
  }
#}

saveRDS(all_result,"~/analysis/circadian/R/Correlation.features2hormone.rds")