library("Seurat")
library("MetaCycle")
library("dplyr")
library("tidyverse")
library("groupdata2")

source('/tmpdata/LyuLin/script/circadian/circadian_core.R')

# This is a wrapper to find oscillating genes with a metacell seurat object

# upstream: runSeaCells.sh, runSeaCells.slurm, runSeaCellsForDrop.slurm, runSeaCells.multiplexed.slurm
# downstream: readJTKFromMetaCells
# dependency: NSF
# caller: NSF
##
# speed:

if(!file.exists("/tmpdata/LyuLin/analysis/circadian/R/seacell.0.03.bytype.rds")){
  metacell2srt("/tmpdata/LyuLin/analysis/circadian/R/preparation_seacell_byType.0.03/",
               out.rds.path='/tmpdata/LyuLin/analysis/circadian/R/seacell.0.03.bytype.rds',
               std.cellranger.out = F,empty.drop.mode = F)
}

srt.metacell<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/seacell.0.03.bytype.rds')

srt.metacell<-NormalizeData(srt.metacell,normalization.method = "RC",scale.factor = 1000000)
mat<-as.data.frame(LayerData(srt.metacell,layer = "data"))

mat<-rownames_to_column(mat,"feature")
exp_table<-gather(mat,key="observation",value="value",-feature)
out.dir<-"/tmpdata/LyuLin/analysis/circadian/R/seacell.meta2d.0.03.bytype/"

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
    prefix=paste0(this.individual,"_",gsub(" ","_",this.celltype))
    if(length(list.files(out.dir,pattern = prefix))>=1){
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
    downsample.table=this.exp_table[c("observation","CT")] %>% unique() %>% downsample(.,"CT")
    this.exp_table=dplyr::filter(this.exp_table,observation %in% downsample.table$observation)
    this.exp_table=arrange(this.exp_table,feature,CT,value)
    if(length(downsample.table$CT %>% table())<7){
      message("time points not enough")
      next
    }
    min.n=downsample.table$CT %>% table() %>% min
    if(min.n<3){
      message("n.observation not enough")
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
