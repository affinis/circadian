source('~/script/circadian/circadian_core.R')

# This is a wrapper to find oscillating genes with a seurat object. 
# Different cell types are made into pseudobulk.
# Only method LS would be used

# upstream: TypeCluster
# downstream: readLS
# dependency: NSF
# caller: NSF
##

exp_table<-readRDS('~/analysis/core_data/CYTOF.5individual.pseudobulk.mat.rds')
out.dir<-"~/analysis/R/FineTypeLevelCYTOFPseudobulk/"

if(!dir.exists(out.dir)){
  dir.create(out.dir)
}

exp_table<-exp_table[,!is.nan(colSums(exp_table))]

exp_table<-rownames_to_column(exp_table,"feature")
exp_table<-gather(exp_table,key="observation",value="value",-feature)
exp_table$type<-getField(exp_table$observation,"-",3)
exp_table$CT<-getField(exp_table$observation,"-",2)
exp_table$individual<-getField(exp_table$observation,"-",1)
celltypes<-sort(unique(exp_table$type))
individuals<-sort(unique(exp_table$individual))

for (this.celltype in celltypes) {
  if(this.celltype=="Eryth"){
    next
  }
  for (this.individual in individuals) {
    prefix=paste0(this.individual,"_",gsub(" ","_",this.celltype))
    if(file.exists(file.path(out.dir,paste0("meta2d_",prefix,".txt")))>=1){
      message(paste0(prefix," data exists, skipping.."))
      next
    }else{
      message(paste0(prefix," analysing"))
    }
    this.exp_table=dplyr::filter(exp_table,type==this.celltype,individual==this.individual)
    if(nrow(this.exp_table)==0){
      message("no expression data, skipping")
      next
    }
    this.exp_table=arrange(this.exp_table,feature,CT)

    min.n=1
    rep.time=nrow(this.exp_table)/min.n
    this.exp_table$Rep=rep(paste0(".Rep",1:min.n),rep.time)
    this.exp_table$Rep.full=paste0(this.exp_table$CT,this.exp_table$Rep)
    print(head(this.exp_table))
    this.timepoints=this.exp_table$CT %>% gsub("CT","",.) %>% gsub("^0","",.) %>% as.numeric() %>% sort() %>% unique()
    if(length(this.timepoints)<5){
      next
    }
    this.exp_table=spread(this.exp_table[c("feature","Rep.full","value")],key=Rep.full,value=value)
    colnames(this.exp_table)[1]="geneSymbol"
    #print(head(this.exp_table))
    
    ##
    prefix=paste0(this.individual,"_",gsub(" ","_",this.celltype))
    #this.timepoints=NULL
    #if(this.individual %in% INDIVIDUALS_BATCH1){
    #  if(this.individual=="ZF"){
    #    this.timepoints=rep(c(13,17,21,25,29), each=min.n)
    #  }else{
    #    this.timepoints=rep(c(9,13,17,21,25,29), each=min.n)
    #  }
    #}else{
    #  this.timepoints=rep(c(11,15,19,23,27,31,35), each=min.n)
    #}
    print(this.timepoints)
    write.table(this.exp_table, file=paste0(out.dir,prefix,".txt"),
                sep="\t", quote=FALSE, row.names=FALSE)
    meta2d(infile=paste0(out.dir,prefix,".txt"), filestyle="txt",cycMethod = c('ARS',"JTK","LS"),
           outdir=out.dir, timepoints=this.timepoints,nCores = 40,
           outIntegration="both")
  }
}