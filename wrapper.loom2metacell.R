source('~/script/circadian/circadian_core.R')

# upstream: runSeaCells.slurm, runSeaCellsBySampleByType.slurm, velocyto.array.slurm
# downstream: wrapper.testRhythmicity.R
# dependency: 
# caller:
##
# speed: for 1 individual with ~2400 metacell, 13 min

for(individual in INDIVIDUALS_BATCH1){
  loom2Metacells(individual)
}


spliced_files<-list.files('~/analysis/circadian/R/spliced.mecell',pattern="_spliced.rds",full.names = T)
spliced_srt_list<-list()
i=1
for(spliced_file in spliced_files){
  mat=readRDS(spliced_file)
  rownames(mat)=mat$geneName
  mat$geneName=NULL
  this.srt=CreateSeuratObject(mat)
  spliced_srt_list[[i]]=this.srt
  i=i+1
}

unspliced_files<-list.files('~/analysis/circadian/R/spliced.mecell',pattern="_unspliced.rds",full.names = T)
unspliced_srt_list<-list()
i=1
for(unspliced_file in unspliced_files){
  mat=readRDS(unspliced_file)
  rownames(mat)=mat$geneName
  mat$geneName=NULL
  this.srt=CreateSeuratObject(mat)
  unspliced_srt_list[[i]]=this.srt
  i=i+1
}

spliced_srts<-merge(spliced_srt_list[[1]],spliced_srt_list[2:length(spliced_srt_list)])
unspliced_srts<-merge(unspliced_srt_list[[1]],unspliced_srt_list[2:length(unspliced_srt_list)])

spliced_srts@meta.data$type<-getField(rownames(spliced_srts@meta.data),"_",3)
unspliced_srts@meta.data$type<-getField(rownames(unspliced_srts@meta.data),"_",3)

spliced_srts@meta.data$CT<-getField(rownames(spliced_srts@meta.data),"_",2)
unspliced_srts@meta.data$CT<-getField(rownames(unspliced_srts@meta.data),"_",2)

spliced_srts@meta.data$individual<-getField(rownames(spliced_srts@meta.data),"_",1)
unspliced_srts@meta.data$individual<-getField(rownames(unspliced_srts@meta.data),"_",1)

spliced_srts<-JoinLayers(spliced_srts)
unspliced_srts<-JoinLayers(unspliced_srts)

saveRDS(spliced_srts,"~/analysis/circadian/R/12individual.spliced.metacell.srt.rds")
saveRDS(unspliced_srts,"~/analysis/circadian/R/12individual.unspliced.metacell.srt.rds")
