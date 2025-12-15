source('~/script/circadian/circadian_core.R')


meta2d.result0<-readRDS('~/analysis/core_data/meta2d.pseudobulkByCelltypeByIndividual.rds')

#bulk_by_type<-readRDS("~/analysis/core_data/16individual.pseudobulk.bytype.CPM.rds")
#top.genes<-CIRCADIAN_GENES_MAIN
#predicted.data<-predictIndividualwithTauFisher(data=bulk_by_type,individual = c("KD","JJC","HZD","ZYX","WLG","XSP"),celltype = "CD14_Mono",
#                                               predictor.genes = top.genes,bool.include.trained = T,bool.return.pred = T)

# calculate dispersion of inter-individual phase shift between genes within individual
for (used.celltype in unique(meta2d.result0$celltype)) {
  plotdata=NULL
  message("debug1")
  gene.list=meta2d.result0[meta2d.result0$celltype==used.celltype&meta2d.result0$meta2d_pvalue<=0.05&!grepl("ENSG|LINC",meta2d.result0$CycID)&meta2d.result0$meta2d_rAMP>0.25,"CycID"]
  message("debug2")
  gene.list=gene.list %>% table()
  gene.list=gene.list[gene.list>=2]
  if(length(gene.list)<2){
    next
  }
  gene.list=names(gene.list)
  message(paste0(used.celltype, "detected ",length(gene.list), " reliable circadian genes."))
  for (this.individual in unique(meta2d.result0$individual)) {
    message(this.individual)
    this.meta2d=meta2d.result0[meta2d.result0$individual==this.individual&meta2d.result0$celltype==used.celltype&meta2d.result0$CycID %in% gene.list&meta2d.result0$meta2d_pvalue<0.3,c("CycID","meta2d_phase")]
    for(this.gene in this.meta2d$CycID){
      this.lag=this.meta2d[this.meta2d$CycID==this.gene,"meta2d_phase"] %>% as.numeric()
      result=this.meta2d
      result$CycID2=this.gene
      result$LAG2=this.lag
      result$lag_diff=result[,"meta2d_phase"]-result$LAG2
      result$individual=this.individual
      if(is.null(plotdata)){
        plotdata=result
      }else{
        plotdata=rbind(plotdata,result)
      }
    }
  }
  plotdata<-plotdata %>% group_by(CycID,CycID2) %>% mutate(n.pair=n())
  plotdata$lag_diff<-ifelse(plotdata$lag_diff< -12,plotdata$lag_diff+24,plotdata$lag_diff)
  plotdata$lag_diff<-ifelse(plotdata$lag_diff> 12,plotdata$lag_diff-24,plotdata$lag_diff)
  plotdata$celltype<-used.celltype
  cell.type.str<-gsub(" ","_",used.celltype)
  saveRDS(plotdata,paste0('~/analysis/R/lag.diffs/lag.diff.',cell.type.str,".rds"))
}