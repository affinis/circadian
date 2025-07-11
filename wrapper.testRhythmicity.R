library("Seurat")
library("MetaCycle")
library("dplyr")
library("tidyverse")
library("groupdata2")

# This is a wrapper to find oscillating genes with a metacell seurat object

# upstream: metacell2srt
# downstream: readJTKFromMetaCells
# dependency: NSF
# caller: NSF

srt.metacell<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/seacell.0.08.rds')
#srt<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/7mixed.integrated.annotated.clean.sct.Azimuth.rds')
#srt.metacell<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/seacell.droplet.0.04.rds')
#srt.metacell$CT<-srt.metacell$CT_ZYR

mat<-as.data.frame(LayerData(srt.metacell,layer = "count"))
cell.size<-srt.metacell[["meta.cell.size"]]$meta.cell.size
mat<-as.data.frame(t(t(mat)/cell.size))

mat<-rownames_to_column(mat,"feature")
exp_table<-gather(mat,key="observation",value="value",-feature)
#out.dir<-"/tmpdata/LyuLin/analysis/circadian/R/seacell.droplet.meta2d.0.04.ZYR/"
out.dir<-"/tmpdata/LyuLin/analysis/circadian/R/seacell.meta2d.0.08.normalizedbyMetacellSize/"

if(!dir.exists(out.dir)){
  dir.create(out.dir)
}

meta.table<-srt.metacell@meta.data[,c("CT","type","individual")]
meta.table<-rownames_to_column(meta.table,"observation")
exp_table<-left_join(exp_table,meta.table)
celltypes<-sort(unique(exp_table$type))
individuals<-sort(unique(exp_table$individual))
for (this.celltype in celltypes) {
  for (this.individual in individuals) {
    this.exp_table=dplyr::filter(exp_table,type==this.celltype,individual==this.individual)
    if(nrow(this.exp_table)==0){
      next
    }
    downsample.table=this.exp_table[c("observation","CT")] %>% unique() %>% downsample(.,"CT")
    this.exp_table=dplyr::filter(this.exp_table,observation %in% downsample.table$observation)
    this.exp_table=arrange(this.exp_table,feature,CT,value)
    if(length(downsample.table$CT %>% table())<7){
      next
    }
    min.n=downsample.table$CT %>% table() %>% min
    if(min.n<3){
      next
    }
    rep.time=nrow(this.exp_table)/min.n
    this.exp_table$Rep=rep(paste0(".Rep",1:min.n),rep.time)
    this.exp_table$Rep.full=paste0(this.exp_table$CT,this.exp_table$Rep)
    print(head(this.exp_table))
    this.exp_table=spread(this.exp_table[c("feature","Rep.full","value")],key=Rep.full,value=value)
    colnames(this.exp_table)[1]="geneSymbol"
    print(head(this.exp_table))
    
    ##
    prefix=paste0(this.individual,"_",gsub(" ","_",this.celltype))
    write.table(this.exp_table, file=paste0(out.dir,prefix,".txt"),
                sep="\t", quote=FALSE, row.names=FALSE)
    meta2d(infile=paste0(out.dir,prefix,".txt"), filestyle="txt",cycMethod = c("JTK"),
           outdir=out.dir, timepoints=rep(c(11,15,19,23,3,7,11), each=min.n),nCores = 64,
           outIntegration="noIntegration")
  }
}
