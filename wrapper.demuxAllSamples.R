# This is a wrapper for demultiplexing sequenced circadian samples
# upstream: <bash>callMTVariantsCellLevelForCellrangerMulti.sh, <bash>wrapper.callMTVariantsBulkBWA.sh
# downstream: <circadian_core.R>mergeMultiplexedSamplesByID
# caller: NSF
# dependency: <seuratWrapper.R>seuratWrap1, <seuratWrapper.R>seuratWrap2, 
#             <seuratWrapper.R>seuratWrap3, <demux_core_functions.R>demux
##
# speed: about 5h for 7 multiplexed sample with each of them containing 50000-60000 cells

source('/tmpdata/LyuLin/script/demuxer/src/demux_core_functions.R')
SCRIPTS<-"/tmpdata/LyuLin/script/wrapper_script"
source(paste(SCRIPTS,"seuratWrapper.R",sep="/"))
source(paste(SCRIPTS,"seuratIntegrateWrapper.R",sep="/"))

samples<-c("X5_250409MIX01","X5_250424MIX02","X5_250424MIX03","X5_250424MIX04","X5_250424MIX05","X5_250424MIX06","X5_250424MIX07")
#samples<-c("X5_250424MIX06","X5_250424MIX07")
sample_individual_CT<-list()
sample_individual_CT[["X5_250409MIX01"]]<-c("KD"="CT15","ZYR"="CT15","JJC"="CT19","LYH"="CT15")
sample_individual_CT[["X5_250424MIX02"]]<-c("KD"="CT35","ZYR"="CT27","JJC"="CT15","LYH"="CT27")
sample_individual_CT[["X5_250424MIX03"]]<-c("KD"="CT23","ZYR"="CT23","JJC"="CT35","LYH"="CT31")
sample_individual_CT[["X5_250424MIX04"]]<-c("KD"="CT31","ZYR"="CT11","JJC"="CT31","LYH"="CT23")
sample_individual_CT[["X5_250424MIX05"]]<-c("KD"="CT27","ZYR"="CT31","JJC"="CT23","LYH"="CT11")
sample_individual_CT[["X5_250424MIX06"]]<-c("KD"="CT11","ZYR"="CT35","JJC"="CT27","LYH"="CT35")
sample_individual_CT[["X5_250424MIX07"]]<-c("KD"="CT19","ZYR"="CT19","JJC"="CT11","LYH"="CT19")

for(sample in samples){
  sample_path=paste0("/tmpdata/LyuLin/analysis/circadian/cellranger/",sample,"/outs/per_sample_outs/",sample,"/count/")
  srt=seuratWrap1(paste0(sample_path,"sample_filtered_feature_bc_matrix/"),min.features = 500)
  srt=seuratWrap2(srt)
  srt=seuratWrap3(srt,res = 0.01)
  
  demux_out=demux(sample_path,'/tmpdata/LyuLin/analysis/circadian/RNASeq')
  
  srt=AddMetaData(srt,demux_out)
  srt@meta.data$CT<-case_when(srt@meta.data$individual == "KD" ~ sample_individual_CT[[sample]][["KD"]],
                              srt@meta.data$individual == "ZYR" ~ sample_individual_CT[[sample]][["ZYR"]],
                              srt@meta.data$individual == "JJC" ~ sample_individual_CT[[sample]][["JJC"]],
                              srt@meta.data$individual == "LYH" ~ sample_individual_CT[[sample]][["LYH"]],
                              TRUE ~ "unknown")
  new_cell_ids <- paste0("TF_",srt@meta.data$individual,"_",srt@meta.data$CT,"_", colnames(srt))
  srt@meta.data$orig.ident=rownames(srt@meta.data)
  srt=RenameCells(srt,new.names = new_cell_ids)
  srt$sample=sample
  saveRDS(srt,paste0('/tmpdata/LyuLin/analysis/circadian/R/',sample,'.demuxed.rds'))
}
