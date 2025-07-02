# Set this option BEFORE installing CRAN packages
# options(repos = c(CRAN = "https://cloud.r-project.org",Bioc = BiocManager::repositories()))

library(Seurat)
library(SeuratObject)
# remotes::install_github("satijalab/seurat-data", quiet = TRUE)
library(SeuratData)
# remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
library(devtools)
library(anndata)
library(org.Hs.eg.db)
library(biomaRt)
library(tidyverse)
library(hash)
library(zeitgebr)
# install_github("velocyto-team/velocyto.R")
# library(velocyto.R)
library(pagoda2)
library(groupdata2)
library(RColorBrewer)
library(MetaCycle)
# remotes::install_github("satijalab/azimuth", quiet = TRUE)
library(Azimuth)
library(vegan)
library(doParallel)
library(foreach)
library(clusterProfiler)
library(destiny)
library(ggcorrplot)
# devtools::install_github("micnngo/tauFisher", build_vignettes = TRUE) 
library(tauFisher)
library(VennDiagram)
library(grid)
library(shazam)
library(dowser)
library(scoper)
library(alakazam)
library(DESeq2)
library(ReactomePA)
library(ggrepel)
library(tidygraph)
library(ggraph)
# devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)
library(ggradar)
library(topGO)
library(dorothea)
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
library(DoubletFinder)
# drc used to analyse ELISA data
library(drc)

SCRIPTS<-"/tmpdata/LyuLin/script/wrapper_script"
source(paste(SCRIPTS,"seuratWrapper.R",sep="/"))
source(paste(SCRIPTS,"seuratIntegrateWrapper.R",sep="/"))
source('/tmpdata/LyuLin/script/circadian/ggplot2_core.R')

#:::::::::::::::::::::#
# GLOBAL VARs/SETTINGs
#:::::::::::::::::::::#

setwd('/tmpdata/LyuLin/analysis/circadian/cellranger')
HOME='/tmpdata/LyuLin/analysis/circadian/cellranger'
#INTERNAL_CONTROL<-c("CDK4")
INTERNAL_CONTROL<-c("ACTB","GAPDH","B2M","RPLP0","TBP","HPRT1","PPIA")
time_point<-c(9,13,17,21,1,5)
time_point_inorder<-c(1,5,9,13,17,21)
CT_TIME<-paste0("CT",time_point)
CT_TIME_ORDER<-paste0("CT",time_point_inorder)
PATIENTS<-c("TFSH190500A_HZD","TFSH190500F_LJQ","TFSH190500I_WLG","TFSH190500I_YXQ","TFSH190500K_QGG","TFSH190501F_XAH",
            "TF_SLY","TF_XSP","TF_ZF","TF_ZXL","TF_ZYJ","TF_ZYX")
AGES<-c(49,50,"WLG",50,"QGG","XAH",30,29,"ZF",29,29,30)
EMI_PATIENTS<-c("TFSH190500I_WLG","TFSH190500K_QGG","TFSH190501F_XAH")
HEALTH<-setdiff(PATIENTS,EMI_PATIENTS)
VALID_INDIVIDUALS=HEALTH[HEALTH!="TF_ZF"]
TCR_LOCUS<-c("TRA:v_call","TRA:j_call","TRB:v_call","TRB:d_call","TRB:j_call")
BCR_LOCUS<-c("IGH:v_call","IGH:d_call","IGH:j_call","IGK:v_call","IGK:j_call","IGL:v_call","IGL:j_call")
CIRCADIAN_GENES_MAIN<-c("CLOCK","BMAL1","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2","DBP","TEF","HLF","CIART")

CIRCADIAN_GENES_PMID_38190520<-c("SERTAD1","ZNF101","PHF21A","STMN3","GNG2","IL1B","IL13RA1",
                         "NR1D1","IRS2","FKBP5","ID3","UBE2B","ZNF438","CLEC4E",
                         "ITGA6","DDIT4","DUSP1","STEAP4","GHRL","NR1D2","NELL2",
                         "AK5","FOSL2","DHRS13","MPPE1","CYB561","CAMKK1","PER1",
                         "AVIL","SMAP2","CDC42EP2","NUDT5","EPHX2","MKNK2","SFXN2",
                         "SESN2","PRRG4")

MAIN_GROUP_MARKERS<-c("CD3D","GNLY","MS4A1","JCHAIN","CD14","MS4A2","PECAM1","CD34","PPBP")

SUBSET_MARKERS<-c("MS4A1","CD27","TCL1A","JCHAIN","EGR1","CD79A","CD19","IGHD","IGHA1",
                  "FCN1","CD14","S100A8","FCGR3A","C1QA","FCER1G","LILRA4","CD1C","CLEC9A",
                  "CD3D","CD4","CD8A","TRBV7-6","FOXP3",
                  "CCR7","LEF1","IL7R","TCF7","PTGER2","FGFBP2","TRDC","TRGC1","CD40LG","GPR183","CXCR3",
                  "SLC4A10","MKI67","GNLY","GZMK","NCAM1","LIMS1","ANXA1","GATA3","MX1")

TNK_markers_NI<-c("CD3G","CCR7","SELL","LEF1","IL7R","TCF7","CD4","GPR183","AQP3","ANXA1","GATA3","PLP2",
                  "FOXP3","CD8A","CD8B","NKG7","GNLY","GZMB","GZMA","GZMK","HAVCR2","TIGIT","CTLA4","CMC1",
                  "DUSP2","CCL5","ZNF683","IKZF2","SLC4A10","TRAV1-2","TRDV2","TRGV9","NCAM1","KLRF1","MKI67",
                  "TRAC","FCGR3A")

TNK_markers_ZhangZeMing<-c("CCR7","LEF1","SELL","TCF7","CD27","CD28","S1PR1", #CD8+TN
                           "CCR7","SELL","IL7R","CD27","CD28","PRF1","GZMA","CCL5","GPR183","S1PR1", #CD8+TCM
                           "KLRG1","CX3CR1","FCGRЗA","FGFBP2","PRF1","GZMH","TBX21","EOMES","S1PR1","S1PR5", #CD8+TEMRA/TEFF
                           "GZMK","CXCR4","CXCR3","CD44", #CD8+TEM
                           "CD6","XCL1","XCL2","MYADM","CAPG","RORA","NR4A1","NR4A2","NR4A3","CD69","ITGAE", #CD8+TRM
                           "CD160","KIR2DL4","TMIGD2","KLRC1","KLRC2","KLRC3","NR4A1","NR4A2","NR4A3","IKZF2","ENTPD1","CD69","ITGAE", #CD8+IEL
                           "HAVCR2","CXCL13","PDCD1","LAYN","TOX","IFNG","GZMB","MIR155HG","TNFRSF9","ITGAE", #CD8+TEX
                           "SLC4A10","KLRB1","ZBTB16","NCRЗ","RORC","RORA" #MAIT
                           )

celltypes<-c("B naive"="B","B memory"="B","Plasmablast"="B","B intermediate"="B",
             "CD14 Mono"="Myeloid","CD16 Mono"="Myeloid","ASDC"="Myeloid","cDC1"="Myeloid","cDC2"="Myeloid","pDC"="Myeloid",
             "CD4 Naive"="T","CD4 CTL"="T","CD4 Proliferating"="T","CD4 TCM"="T","CD4 TEM"="T",
             "CD8 Naive"="T","CD8 TCM"="T","CD8 TEM"="T","dnT"="T","gdT"="T","ILC"="T","MAIT"="T",
             "NK"="NK","NK_CD56ᵇʳⁱᵍʰᵗ"="NK","NK Proliferating"="NK",
             "HSPC"="HSPC")
celltypes_NI<-c("B_naive","CD14⁺Monocytes","CD16⁺Monocytes","CD4_naive_CCR7","CD4_TCM_AQP3",
                "CD4_TEM_ANXA1","CD4_TEM_GNLY","CD4_TEM_GZMK","CD4_Treg_FOXP3","CD8_MAIT_SLC4A10",
                "CD8_naive_LEF1","CD8_TEM_CMC1","CD8_TEM_GNLY","CD8_TEM_ZNF683","dnT_LYST",
                "HSC_CD34","mDC","NK_CD56ᵇʳⁱᵍʰᵗ","NK_CD56ᵈⁱᵐ","pDC","Plasma cell","TNK_proliferatig_MKI67","γδT")

BATCH<-c("JJC-CT11"="X5_250424MIX07","JJC-CT15"="X5_250424MIX02","JJC-CT19"="X5_250409MIX01","JJC-CT23"="X5_250424MIX05","JJC-CT27"="X5_250424MIX06","JJC-CT31"="X5_250424MIX04","JJC-CT35"="X5_250424MIX03",
         "LYH-CT11"="X5_250424MIX05","LYH-CT15"="X5_250409MIX01","LYH-CT19"="X5_250424MIX07","LYH-CT23"="X5_250424MIX04","LYH-CT27"="X5_250424MIX02","LYH-CT31"="X5_250424MIX03","LYH-CT35"="X5_250424MIX06",
         "KD-CT11"="X5_250424MIX06","KD-CT15"="X5_250409MIX01","KD-CT19"="X5_250424MIX07","KD-CT23"="X5_250424MIX03","KD-CT27"="X5_250424MIX05","KD-CT31"="X5_250424MIX04","KD-CT35"="X5_250424MIX02",
         "ZYR-CT11"="X5_250424MIX04","ZYR-CT15"="X5_250409MIX01","ZYR-CT19"="X5_250424MIX07","ZYR-CT23"="X5_250424MIX03","ZYR-CT27"="X5_250424MIX02","ZYR-CT31"="X5_250424MIX05","ZYR-CT35"="X5_250424MIX06")

#::::::::: ::::#
# I/O FUNCTIONS
#:::::::::: :::#

# Function: readAllMetaCell
# create seurat object with metacell
# file.path: a directory containing seacells outputs, each dir relates to a sample
# std.cellranger.out: if seacells was run under cellranger multi path (like 'SAMPLEID/outs/per_sample_outs/SAMPLEID/count/seacells/'), set it TRUE, else false
##
# upstream: <slurm>runSeaCells.slurm
# downstream: <Rscript>wrapper.testRhythmicity.R
# dependency: NSF
# caller: NSF
readAllMetaCell<-function(file.path,std.cellranger.out=T){
  if(std.cellranger.out){
    samples=PATIENTS %>% lapply(.,paste0,"_",1:6) %>% unlist()
  }else{
    samples=system(paste0('find ',file.path,' -name "*_metacells.csv"'),intern = T)
  }
  mat_c=NULL
  meta_c=NULL
  for (sample in samples) {
    mat=NULL
    meta.path=NULL
    if(std.cellranger.out){
      prefix=paste0(getField(sample,"_",c(1,2)),"_",CT_TIME[as.numeric(getField(sample,"_",3))])
      mat.path=paste0(file.path,sample,"/outs/per_sample_outs/",sample,"/count/seacells/",prefix,"_metacells.csv")
      if(!file.exists(mat.path)){
        next
      }
      mat=read.csv(mat.path)
    }else{
      prefix=paste0("TF_",getField(basename(sample),"_",c(1,2)))
      mat=read.csv(sample)
    }
    
    colnames(mat)=gsub(".","-",colnames(mat),fixed = T)
    colnames(mat)[-1]=paste0(prefix,"_",colnames(mat)[-1])
    #print(head(mat[1:5,1:5]))
    
    if(std.cellranger.out){
      meta.path=paste0(file.path,sample,"/outs/per_sample_outs/",sample,"/count/seacells/",prefix,"_metadata.csv")
    }else{
      meta.path=gsub("metacells","metadata",sample)
    }
    
    meta=read.csv(meta.path)
    rownames(meta)=meta$index
    meta$index=NULL
    #print(head(meta[1:5,1:5]))
    
    meta.suma=meta[c("predicted.celltype.l2","SEACell")] %>% table() %>% as.data.frame()
    meta.suma=dplyr::filter(meta.suma,predicted.celltype.l2!="")
    meta.suma=meta.suma %>% group_by(SEACell) %>% mutate(meta.cell.size=sum(Freq))
    meta.suma=meta.suma %>% group_by(SEACell) %>% mutate(predicted.celltype.l2.purity=Freq/meta.cell.size)
    meta.suma=meta.suma %>% group_by(SEACell) %>% mutate(predicted.celltype.l2.main=predicted.celltype.l2[which.max(Freq)])
    meta.suma=dplyr::filter(meta.suma,predicted.celltype.l2==predicted.celltype.l2.main)[c("SEACell","meta.cell.size","predicted.celltype.l2.purity","predicted.celltype.l2.main")]
    
    meta.suma.NI=meta[c("manual_NI","SEACell")] %>% table() %>% as.data.frame()
    meta.suma.NI=dplyr::filter(meta.suma.NI,manual_NI!="")
    meta.suma.NI=meta.suma.NI %>% group_by(SEACell) %>% mutate(meta.cell.size=sum(Freq))
    meta.suma.NI=meta.suma.NI %>% group_by(SEACell) %>% mutate(manual_NI.purity=Freq/meta.cell.size)
    meta.suma.NI=meta.suma.NI %>% group_by(SEACell) %>% mutate(manual_NI.main=manual_NI[which.max(Freq)])
    meta.suma.NI=dplyr::filter(meta.suma.NI,manual_NI==manual_NI.main)[c("SEACell","manual_NI.purity","manual_NI.main")]
    
    meta.suma=left_join(meta.suma,meta.suma.NI,by="SEACell") %>% as.data.frame()
    rownames(meta.suma)=meta.suma$SEACell
    meta.suma$SEACell=NULL
    #print(head(meta.suma[1:3,]))
    
    if(is.null(mat_c)){
      mat_c=mat
    }else{
      mat_c=left_join(mat_c,mat)
    }
    
    if(is.null(meta_c)){
      meta_c=meta.suma
    }else{
      meta_c=rbind(meta_c,meta.suma)
    }
    #srt=CreateSeuratObject(counts=mat)
    #srt=AddMetaData(srt,meta.suma)
    #return(srt)
  }
  mat_c=as.data.frame(mat_c)
  rownames(mat_c)=mat_c$X
  mat_c$X=NULL
  mat_c[is.na(mat_c)]=0
  print(mat_c[1:3,1:3])
  srt=CreateSeuratObject(mat_c)
  if(!std.cellranger.out){
    rownames(meta_c)=gsub("filtered_bc_matrix_","",rownames(meta_c)) %>% gsub("^","TF_",.)
  }
  print(meta_c[1:3,1:3])
  srt=AddMetaData(srt,meta_c)
  srt$type=srt$predicted.celltype.l2.main
  srt$CT=getField(rownames(srt@meta.data),sep = "_",3)
  srt$individual=getField(rownames(srt@meta.data),sep = "_",2)
  srt=AddMetaData(srt,metadata=colSums(LayerData(srt)),col.name = "nCount_RNA")
  return(srt)
}

# Function: readJTKFromMetaCells
# get all meta cell JTK result by individual
# caller: NSF
# dependency: getField
# upstream: wrapper.testRhythmicity.R
# downstream: 
readJTKFromMetaCells<-function(path){
  filesJTK=list.files(path,pattern = "JTKresult_*")
  res=NULL
  for(file in filesJTK){
    this.individual=getField(file,"_",2)
    this.celltype=getField(file,"_",-c(1,2)) %>% gsub(".txt","",.,fixed = T) %>% gsub("_"," ",.)
    message(paste0("reading data from individual: ",this.individual,", celltype: ",this.celltype))
    this.res=read.delim(file.path(path,file),header = T)
    this.res$individual=this.individual
    this.res$celltype=this.celltype
    if(is.null(res)){
      res=this.res
    }else{
      res=rbind(res,this.res)
    }
  }
  print(head(res))
  return(res)
}

# Function: readMergedJTK
# get merged JTK_cycle outs
# outpath: a folder path containing block1, block2 ...
# wanted.celltypes: desired celltype
# caller: NSF
# dependency: getJTK_CYCLEouts
# upstream: NSF
# downstream: 
readMergedJTK<-function(outpath,wanted.celltypes){
  JTK_CYCL=NULL
  for(block_data in list.files(outpath)){
    message(block_data)
    block_path=paste0(outpath,"/",block_data)
    message(block_path)
    block=getJTK_CYCLEouts(block_path,celltypes=wanted.celltypes,filter.data = F)
    if(is.null(JTK_CYCL)){
      JTK_CYCL=block
    }else{
      JTK_CYCL=rbind(JTK_CYCL,block)
    }
  }
  return(JTK_CYCL)
}

# Function: generateAnnotationFile
# merge cell annotations
# input: annotated srt *.rds file (TNK.file,B.file,Myeloid.file)
# output: *.tsv
# note that Platelet is not included so for cells whose toplevel annotation is Platelet, it will not be included in the final tsv file.
# cols: colnames wanted in seurat object@meta.data
##
# caller: NSF
# dependency: NSF
# upstream: TypeCluster
# downstream: <seurat>AddMetaData, <bash>runSeaCells.sh
generateAnnotationFile<-function(TNK.file=NULL,B.file=NULL,Myeloid.file=NULL,Stem.file=NULL,full.file=NULL,cols,new.colnames=NULL,out.path='../analysis/cell.annotations.tsv'){
  all.meta=NULL
  if(!is.null(full.file)){
    if(is.character(full.file)){
      full=readRDS(full.file)
      all.meta=full@meta.data[,cols]
    }else{
      all.meta=full.file@meta.data[,cols]
    }
  }else{
    message("reading TNK file")
    TNK=readRDS(TNK.file)
    TNK.meta=TNK@meta.data[,cols]
    remove(TNK)
    message("reading B file")
    B=readRDS(B.file)
    B.meta=B@meta.data[,cols]
    remove(B)
    message("reading Myeloid file")
    Myeloid=readRDS(Myeloid.file)
    Myeloid.meta=Myeloid@meta.data[,cols]
    remove(Myeloid)
    message("reading Stem file")
    Stem=readRDS(Stem.file)
    Stem.meta=Stem@meta.data[,cols]
    remove(Stem)
    all.meta=base::Reduce(rbind,x=list(TNK.meta,B.meta,Myeloid.meta,Stem.meta))
    if(length(new.colnames)>=1){
      colnames(all.meta)=new.colnames
    }
  }
  all.meta=rownames_to_column(all.meta,var = "cell_id")
  write_delim(all.meta,out.path,delim = '\t')
}

# Function: getWholeDaySplicedData
# get a day's (6 samples) spliced data from an individual 
# patient_id: example: TFSH190500A_HZD
# cells: a vector of cell barcode to include in final data, default is NULL (all cells included)
# output: return a list of loom Matrix
##
# caller: getMultipleWholeDaySplicedData
# dependency: subsetLoomMatrix
# upstream: <slurm>velocyto.array.slurm
# downstream: mergeSplicedData (merge matrix of different time points into a single one)
getWholeDaySplicedData<-function(patient_id,cells=NULL){
  res=list()
  for(i in 1:6){
    sample_id=paste0(patient_id,"_",i)
    sample_loom_path=paste0(sample_id,"/outs/per_sample_outs/",sample_id,"/count/velocyto/")
    sample_loom_file_name=list.files(sample_loom_path)[1]
    sample_loom_full_path=paste0(sample_loom_path,sample_loom_file_name)
    message("loom file: ")
    message(sample_loom_full_path)
    loomMatrix=read.loom.matrices(sample_loom_full_path)
    loomMatrix=subsetLoomMatrix(loomMatrix,cell.preffix=CT_TIME[i],cells=cells)
    colnames(loomMatrix$spliced)=paste(patient_id,colnames(loomMatrix$spliced),sep='_')
    colnames(loomMatrix$unspliced)=paste(patient_id,colnames(loomMatrix$unspliced),sep='_')
    colnames(loomMatrix$ambiguous)=paste(patient_id,colnames(loomMatrix$ambiguous),sep='_')
    res[[i]]=loomMatrix
  }
  return(res)
}

# Function: getMultipleWholeDaySplicedData
# get a day's (6 samples) spliced data from multiple individual 
# individuals: a vector of string indicating different individuals, example: TFSH190500A_HZD
# output: return a list of loom Matrix
##
# caller: NSF
# dependency: getWholeDaySplicedData
# upstream: NSF
# downstream: mergeSplicedData (merge matrix of different time points into a single one)
getMultipleWholeDaySplicedData<-function(individuals){
  res=list()
  for(individual in individuals){
    res.block=getWholeDaySplicedData(individual)
    res=c(res,res.block)
  }
  return(res)
}

# Function: mergeAllSamples
# read RDS of seurat object which annotated at top level
# inputs are .rds which contained in argument data.path
##
# upstream: TypeCluster
# downstream: dimentionalReductionMergedSrt
# dependency: getSamplesOneDayToList
# caller: NSF
mergeAllSamples<-function(data.path="/lustre/home/acct-medll/medll/data/analysis/"){
  suffix=NULL
  srts=list()
  for(patient in PATIENTS){
    file.exist=list.files(data.path) %>% grep(patient,.,value = T) %>% getField(.,"_",3) %>% gsub(".levelTop.rds","",.) %>% as.numeric()
    suffix.patient=paste0(patient,"_",CT_TIME[file.exist])
    suffix=c(suffix,suffix.patient)
    srts=c(srts,getSamplesOneDayToList(patient.id=patient,data.path=data.path))
  }
  message(paste0("num of seurat objects: ",length(srts),", prefix number: ",length(suffix)))
  res=merge(srts[[1]],srts[2:length(srts)],add.cell.ids=suffix)
  return(res)
}

# Function: mergeSamplesOneDay
# merge samples from a single individual, seldom used, this function was created at the very beginning of this project
# patient_id: example: TFSH190500A_HZD
# data.path: the files as input are .rds contained in argument 'data.path'
##
# upstream: TypeCluster
# downstream: plotGeneExpressionByCT
# dependency: getSamplesOneDayToList
# caller: NSF
##
# return merged seurat object
mergeSamplesOneDay<-function(patient.id,data.path="/lustre/home/acct-medll/medll/data/analysis/"){
  samples=paste0(patient.id,"_",1:6,".levelTop.rds")
  suffix=paste0("CT",time_point)
  srts=getSamplesOneDayToList(patient.id=patient.id,data.path=data.path)
  res=merge(srts[[1]],srts[2:6],add.cell.ids=suffix)
  return(res)
}

# Function: mergeRawSamplesOneDay
# merge srt samples from 'patient.id', the difference between this function and 'mergeSamplesOneDay' is
# that this function read from cellranger multi output rather than processed seurat object, it also add module score
# to each cell and annotate the cell with most-like cell type in column 'type', intermediate seurat object of
# each sample will be write to path specified by 'write.path'
# attention doublet will not be removed
##
# upstream: cellranger multi
# downstream: plotGeneExpressionByCT
# dependency: <seuratWrapper.R>seuratWrap1, <seuratWrapper.R>seuratWrap2, <seuratWrapper.R>seuratWrap3, CT_TIME
# caller: NSF
# return merged srt
mergeRawSamplesOneDay<-function(patient.id,data.path="/lustre/home/acct-medll/medll/data/cellranger_out/",save.each.sample=F,over.write=F,
                                write.path="/lustre/home/acct-medll/medll/data/analysis/"){
  srts=list()
  samples_file=paste0(patient.id,"_",1:6,".levelTop.rds")
  samples_file=paste0(write.path,samples_file)
  for(i in 1:6){
    sample_id=paste0(patient.id,"_",i)
    sprt_path=paste0(sample_id,"/outs/per_sample_outs/",sample_id,"/count/sample_filtered_feature_bc_matrix/")
    suffix=CT_TIME[i]
    srt=seuratWrap1(sprt_path)
    srt=seuratWrap2(srt)
    srt=seuratWrap3(srt,res = 0.05)
    srt=AddModuleScore(srt,list("T"=c("CD3D"),"NK"=c("GNLY"),"Plasma"=c("JCHAIN"),"B"=c("MS4A1"),
                                "Myeloid"=c("CD14","FCN1"),"Stem"=c("CD34","THY1","ITGA2"),"Platelet"=c("PPBP"),"Granular"=c("CCR3","IL5RA"),
                                "Plasmatoid"=c("LILRA4")))
    colnames(srt@meta.data)[7:15]=c("T","NK","Plasma","B","Myeloid","Stem","Platelet","Granular","Plasmatoid")
    srt@meta.data$type=NA
    srt@meta.data$type.top.level=NA
    for(j in 0:(levels(srt@meta.data$seurat_clusters) %>% length())-1){
      srt@meta.data$type[srt@meta.data$seurat_clusters==j]=srt@meta.data[srt@meta.data$seurat_clusters==j,
                      c("T","NK","Plasma","B","Myeloid","Stem","Platelet","Granular","Plasmatoid")] %>% colSums() %>% which.max() %>% names()
    }
    srt@meta.data$CT=suffix
    srt@meta.data$patient=patient.id
    srt@meta.data[srt@meta.data$type %in% c("T","NK"),"type.top.level"]="TNK"
    srt@meta.data[srt@meta.data$type %in% c("B","Plasma"),"type.top.level"]="B"
    srt@meta.data[srt@meta.data$type %in% c("Granular","Myeloid","Plasmatoid"),"type.top.level"]="Myeloid"
    srt@meta.data[srt@meta.data$type %in% c("Platelet"),"type.top.level"]="Platelet"
    srt@meta.data[srt@meta.data$type %in% c("Stem"),"type.top.level"]="HSC"
    if(save.each.sample){
      if(file.exists(samples_file[i])&&!over.write){
        stop("File exists, set over.write=T or change path.")
      }
      saveRDS(srt,samples_file[i])
    }
    srts[[i]]=srt
  }
  res=merge(srts[[1]],srts[2:6],add.cell.ids=CT_TIME)
  return(res)
}

#input: srt rds
#return a list of seurat object
getSamplesOneDayToList<-function(patient.id,data.path="/lustre/home/acct-medll/medll/data/analysis/"){
  samples=paste0(patient.id,"_",1:6,".levelTop.rds")
  srts=list()
  i=1
  j=1
  for(sample in samples){
    full.path=paste0(data.path,sample)
    if(!file.exists(full.path)){
      j=j+1
      next
    }
    srt=readRDS(full.path)
    srt@meta.data$CT=CT_TIME[j]
    srt@meta.data$patient=patient.id
    srts[[i]]=srt
    j=j+1
    i=i+1
  }
  return(srts)
}



#::::::::::::::::::#
# DATA MODIFICATION
#::::::::::::::::::#
mergeByRownames<-function(df1,df2){
  df<-merge(df1,df2,by="row.names",all.x=T,all.y=T)
  rownames(df)<-df$Row.names
  df$Row.names<-NULL
  return(df)
}

#add cluster annotation to seurat object's meta.data
TypeCluster<-function(srt,cluster,type,col.name="seurat_clusters",new.meta="type"){
  srt@meta.data[srt@meta.data[,col.name]%in%cluster,new.meta]<-type
  return(srt)
}

# Function: dimentionalReductionMergedSrt
# a preparation step for integration
##
# srt.merged: a merged srt object
##
# upstream: mergeMultiplexedSamplesByID, mergeSamplesOneDay
# downstream: integrateMergedSrt
# dependency: NSF
# caller: NSF
dimentionalReductionMergedSrt<-function(srt.merged){
  srt.merged=JoinLayers(srt.merged)
  srt.merged=NormalizeData(srt.merged)
  srt.merged=ScaleData(srt.merged)
  srt.merged=FindVariableFeatures(srt.merged)
  srt.merged=RunPCA(srt.merged)
  srt.merged=RunUMAP(srt.merged, dims=1:30)
  return(srt.merged)
}

# Function: integrateMergedSrt
# integrate data merged from several srt, should be processed by dimentionalReductionMergedSrt
##
# srt.merged.dim: merged and dimentional reductioned srt
##
# upstream: dimentionalReductionMergedSrt
# downstream: TypeCluster
# dependency: NSF
# caller: NSF
integrateMergedSrt<-function(srt.merged.dim){
  srt.merged.dim[["RNA"]]=split(srt.merged.dim[["RNA"]], f = srt.merged.dim$sample)
  srt.merged.dim=IntegrateLayers(object=srt.merged.dim, method=CCAIntegration, 
                                  orig.reduction="pca", new.reduction="integrated.cca",verbose = FALSE)
  srt.merged.dim[["RNA"]]=JoinLayers(srt.merged.dim[["RNA"]])
  srt.merged.dim=FindNeighbors(srt.merged.dim, reduction = "integrated.cca", dims = 1:30)
  srt.merged.dim=FindClusters(srt.merged.dim, resolution = 0.1)
  srt.merged.dim=RunUMAP(srt.merged.dim, dims = 1:30, reduction = "integrated.cca")
  return(srt.merged.dim)
}

# Function: dimentionalReductionSubsetSrt
# re-perform dimentional reduction to a subseted srt, prepared for re-integration
##
# upstream: <Seurat>subset
# downstream: integrateSubset
# dependency: NSF
# caller: NSF
dimentionalReductionSubsetSrt<-function(srt.sub,res=0.5){
  srt.sub=NormalizeData(srt.sub)
  srt.sub=ScaleData(srt.sub)
  srt.sub=FindVariableFeatures(srt.sub)
  srt.sub=RunPCA(srt.sub)
  srt.sub=RunUMAP(srt.sub, dims=1:30)
  srt.sub=FindNeighbors(srt.sub,dims =1:30,reduction = "pca")
  srt.sub=FindClusters(srt.sub, resolution = res)
  return(srt.sub)
}

# Function: integrateSubset
# re-integarte data after subsetting
##
# upstream: dimentionalReductionSubsetSrt
# downstream: TypeCluster
# dependency: NSF
# caller: NSF
integrateSubset<-function(subsrt,method=CCAIntegration,classic.normalize=T,k.weight=100){
  subsrt[["RNA"]]=split(subsrt[["RNA"]], f = subsrt$sample)
  if(classic.normalize){
    subsrt=NormalizeData(subsrt)
    subsrt=FindVariableFeatures(subsrt)
    subsrt=ScaleData(subsrt)
  }else{
    subsrt=SCTransform(subsrt)
  }
  subsrt=RunPCA(subsrt)
  if(classic.normalize){
    subsrt=IntegrateLayers(object=subsrt, method=method, orig.reduction="pca", new.reduction="integrated.cca",verbose = FALSE,k.weight=k.weight)
    subsrt[["RNA"]]=JoinLayers(subsrt[["RNA"]])
  }else{
    subsrt=IntegrateLayers(object=subsrt, method=method, orig.reduction="pca", normalization.method = "SCT",k.weight=k.weight,
                           new.reduction="integrated.cca",verbose = FALSE)
    subsrt[["RNA"]]=JoinLayers(subsrt[["RNA"]])
  }
  subsrt=FindNeighbors(subsrt, reduction = "integrated.cca", dims = 1:30)
  subsrt=FindClusters(subsrt, resolution = 0.2)
  subsrt=RunUMAP(subsrt, dims = 1:30, reduction = "integrated.cca")
  return(subsrt)
}

#subset a loom matrix by cells and/or genes, we can add suffix to barcodes here
#input are loom matrix from "read.loom.matrices"
#return is loom matrix
subsetLoomMatrix<-function(loom.matrix,cell.preffix=NULL,cells=NULL,genes=NULL){
  colnames(loom.matrix$spliced)=gsub("sample_alignments_.*:","",colnames(loom.matrix$spliced)) %>% gsub("x$","-1",.)
  colnames(loom.matrix$unspliced)=gsub("sample_alignments_.*:","",colnames(loom.matrix$unspliced)) %>% gsub("x$","-1",.)
  colnames(loom.matrix$ambiguous)=gsub("sample_alignments_.*:","",colnames(loom.matrix$ambiguous)) %>% gsub("x$","-1",.)
  if(!is.null(cells)&&!is.null(genes)){
    loom.matrix$spliced=loom.matrix$spliced[rownames(loom.matrix$spliced) %in% genes,colnames(loom.matrix$spliced) %in% cells]
    loom.matrix$unspliced=loom.matrix$unspliced[rownames(loom.matrix$unspliced) %in% genes,colnames(loom.matrix$unspliced) %in% cells]
    loom.matrix$ambiguous=loom.matrix$ambiguous[rownames(loom.matrix$ambiguous) %in% genes,colnames(loom.matrix$ambiguous) %in% cells]
  }else if(!is.null(cells)&&is.null(genes)){
    loom.matrix$spliced=loom.matrix$spliced[,colnames(loom.matrix$spliced) %in% cells]
    loom.matrix$unspliced=loom.matrix$unspliced[,colnames(loom.matrix$unspliced) %in% cells]
    loom.matrix$ambiguous=loom.matrix$ambiguous[,colnames(loom.matrix$ambiguous) %in% cells]
  }else if(is.null(cells)&&!is.null(genes)){
    loom.matrix$spliced=loom.matrix$spliced[rownames(loom.matrix$spliced) %in% genes,]
    loom.matrix$unspliced=loom.matrix$unspliced[rownames(loom.matrix$unspliced) %in% genes,]
    loom.matrix$ambiguous=loom.matrix$ambiguous[rownames(loom.matrix$ambiguous) %in% genes,]
  }else{
    message("No subsetting performed, return whole matrix directly.")
  }
  if(!is.null(cell.preffix)){
    colnames(loom.matrix$spliced)=gsub("^",paste0(cell.preffix,"_"),colnames(loom.matrix$spliced))
    colnames(loom.matrix$unspliced)=gsub("^",paste0(cell.preffix,"_"),colnames(loom.matrix$unspliced))
    colnames(loom.matrix$ambiguous)=gsub("^",paste0(cell.preffix,"_"),colnames(loom.matrix$ambiguous))
  }
  return(loom.matrix)
}


# Function: mergeMultiplexedSamplesByID
# this is a function for merge multiplexed samples in this study
# file.path: a folder containing multiple rds file to be merged
# sample.ids: sample IDs like 'X5_250409MIX01'
##
# upstream: <Rscript>wrapper.demuxAllSamples.R
# downstream: dimentionalReductionMergedSrt
# dependency: NSF
# caller: NSF
##
# speed: about 1 min for 1G rds
mergeMultiplexedSamplesByID<-function(dir.path,sample.ids){
  message(date())
  files=list.files(dir.path) %>% grep(paste(sample.ids,collapse="|"),.,value = T) %>% grep('rds$',.,value=T)
  all.data=list()
  i=1
  for(file in files){
    this.srt=readRDS(file.path(dir.path,file))
    this.srt=subset(this.srt,individual!="unknown")
    all.data[[i]]=this.srt
    i=i+1
  }
  all.data=merge(all.data[[1]],all.data[2:length(all.data)])
  message(date())
  return(all.data)
}

#::::::::::::::::::::::#
#VISUALIZATION FUNCTIONS
#::::::::::::::::::::::#
#plot Segment pairing by CT
plotSegmentPairingByCT<-function(res.immacantation,seg,cell.type=NULL){
  if(!is.null(cell.type)){
    res.immacantation=dplyr::filter(res.immacantation,cell_type==cell.type)
  }
  res.immacantation.pair.info=table(res.immacantation$cell_id)
  res.immacantation.pair.info=res.immacantation.pair.info[res.immacantation.pair.info==2]
  res.immacantation=dplyr::filter(res.immacantation,cell_id %in% names(res.immacantation.pair.info))
  res.immacantation=res.immacantation[c("v_call","d_call","j_call","cell_id")]
  if(grepl("TR",res.immacantation$v_call[1])){
    s_chain=res.immacantation[grepl("TRA",res.immacantation$v_call),]
    l_chain=res.immacantation[grepl("TRB",res.immacantation$v_call),]
    res.immacantation=left_join(s_chain,l_chain,by="cell_id")
    res.immacantation$IR="TCR"
  }else{
    s_chain=res.immacantation[grepl("IGH",res.immacantation$v_call),]
    l_chain=res.immacantation[grepl("IGK|IGL",res.immacantation$v_call),]
    res.immacantation=left_join(s_chain,l_chain,by="cell_id")
    res.immacantation$IR="BCR"
  }
   print(head(res.immacantation))
   res.immacantation$CT=getField(res.immacantation$cell_id,"_",3)
   res.immacantation=res.immacantation[apply(res.immacantation, 1, function(row) any(grepl(seg, row, ignore.case = TRUE))), ]
   res.immacantation=balanceObsevations(res.immacantation,main.group = "IR",balance.group = "CT")
   #ggplot(res.immacantation)+geom_bar(aes(x=CT,))
}

#plot segment usage
plotSegmentUsageByCT<-function(res.immacantation,chain,seg,cell_restriction=NULL,scale=T,remove.na=T,angle=60){
  if(!is.null(cell_restriction)){
    res.immacantation=dplyr::filter(res.immacantation,cell_type==cell_restriction)
  }
  data=dplyr::filter(res.immacantation,locus==chain) %>% balanceObsevations(.,"locus","CT")
  data=data %>% group_by(CT,across(seg)) %>% summarise(count=n()) %>% arrange(by=count)
  data=data[!is.na(data[,seg] %>% unlist() %>% as.vector()),]
  seg_order=(data %>% group_by(across(seg)) %>% summarise(sum=sum(count)) %>% arrange(desc(sum)))[,seg] %>% unlist() %>% as.vector()
  data=data %>% group_by(across(seg)) %>% arrange(across(seg))
  print(head(data,12))
  if(remove.na){
    na.segments.data=data[seg] %>% table()
    na.segments=na.segments.data[na.segments.data<6] %>% names()
    if(length(na.segments)>0){
      message(paste0("Segment: ",na.segments, " have NA data points."))
      #return(data)
      data=data[!(data[seg] %>% unlist() %in% na.segments),]
      seg_order=seg_order[!(seg_order %in% na.segments)]
    }
  }
  data.scale=reframe(data,z_score=scale(count))
  print(head(data.scale,12))
  data$z_score=data.scale$z_score[,1]
  print(head(seg_order))
  if(scale){
    ggplot(data)+geom_tile(aes(x=get(seg),y=CT,fill=z_score))+scale_y_discrete(limits=rev(CT_TIME_ORDER))+
      scale_x_discrete(limits=seg_order)+xlab("")+ylab("")+
      scale_fill_gradient2(low = "#507AAF",high = "#BE5C37")+
      theme(axis.text.x = element_text(angle = angle,hjust=1,color="black"),
            axis.text.y = element_text(color="black"))+
      ggtitle(paste0(cell_restriction," ",chain,":",seg))
  }else{
    ggplot(data)+geom_tile(aes(x=get(seg),y=CT,fill=count))+scale_y_discrete(limits=rev(CT_TIME_ORDER))+
      scale_x_discrete(limits=seg_order)+xlab("")+ylab("")+
      #scale_fill_gradient2(low = "#507AAF",high = "#BE5C37")+
      theme(axis.text.x = element_text(angle = angle,hjust=1,color="black"),
            axis.text.y = element_text(color="black"))+
      ggtitle(paste0(cell_restriction," ",chain,":",seg))
  }
}



#input: merged spliced data
plotSplicingCounts<-function(loom.matrix,features=CIRCADIAN_GENES_MAIN,cells=NULL,summarise="sum"){
  if(is.null(cells)){
    cells=colnames(loom.matrix$spliced)
  }
  spliced=loom.matrix$spliced[features,cells] %>% as.data.frame() %>% 
    rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  spliced$status="spliced"
  unspliced=loom.matrix$unspliced[features,cells] %>% as.data.frame() %>% 
    rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  unspliced$status="unspliced"
  ambiguous=loom.matrix$ambiguous[features,cells] %>% as.data.frame() %>% 
    rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  ambiguous$status="ambiguous"
  data=rbind(spliced,unspliced)
  print(head(data))
  data$time=getField(data$cell,sep="_",3)
  if(summarise=="mean"){
    plotdata=data %>% group_by(.,feature,status,time) %>% summarise(.,mean_count=mean(count))
    ylab="mean_count"
  }else{
    plotdata=data %>% group_by(.,feature,status,time) %>% summarise(.,total_count=sum(count))
    ylab="total_count"
  }
  print(head(plotdata))
  ggplot(plotdata)+geom_point(aes(x=time,y=get(ylab),color=status))+
    geom_line(aes(x=time,y=get(ylab),group=status,color=status))+
    facet_wrap(~feature,scales = "free_y")+ylab(ylab)+
    scale_x_discrete(limits=CT_TIME_ORDER)+theme(axis.text.x = element_text(angle=45,hjust=1))
}

#input: merged spliced data
#do not use normalize_at 'sample' level, code logic error, as input is merged spliced data, normalize at sample level should have been done
plotSingleSplicingCounts<-function(loom.matrix,feature=CIRCADIAN_GENES_MAIN[1],cells=NULL,summarise="mean",normalize=F,normalize_at="cell",return.data=FALSE){
  if(is.null(cells)){
    cells=colnames(loom.matrix$spliced)
  }
  spliced=loom.matrix$spliced[feature,cells] %>% t %>% as.data.frame() %>% 
    rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  if(normalize){
    if(normalize_at=="cell"){
      nCount_cell_spliced=loom.matrix$spliced[,cells] %>% t %>% as.data.frame() %>% rowSums()
      spliced$count=spliced$count*10000/nCount_cell_spliced
    }else if(normalize_at=="sample"){
      nCount_cell_spliced=loom.matrix$spliced[,cells] %>% sum()
      spliced$count=spliced$count*10000000/nCount_cell_spliced
    }else{
      stop("normalize_at should be 'cell' or 'sample'")
    }
  }
  spliced$individual=getField(spliced$cell,sep="_",c(1,2))
  spliced$status="spliced"
  unspliced=loom.matrix$unspliced[feature,cells] %>% t %>% as.data.frame() %>% 
    rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  if(normalize){
    if(normalize_at=="cell"){
      nCount_cell_unspliced=loom.matrix$unspliced[,cells] %>% t %>% as.data.frame() %>% rowSums()
      unspliced$count=unspliced$count*10000/nCount_cell_unspliced
    }else if(normalize_at=="sample"){
      nCount_cell_unspliced=loom.matrix$unspliced[,cells] %>% sum()
      unspliced$count=unspliced$count*10000000/nCount_cell_unspliced
    }else{
      stop("normalize_at should be 'cell' or 'sample'")
    }
  }
  unspliced$individual=getField(unspliced$cell,sep="_",c(1,2))
  unspliced$status="unspliced"
  
  #total=(loom.matrix$spliced[feature,cells]+loom.matrix$unspliced[feature,cells]+
  #       loom.matrix$ambiguous[feature,cells]) %>% t %>% as.data.frame() %>% 
  #       rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  #total$individual=getField(total$cell,sep="_",c(1,2))
  #total$status="total"
  #ambiguous=loom.matrix$ambiguous[features,cells] %>% as.data.frame() %>% 
  #  rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  #ambiguous$status="ambiguous"
  data=rbind(spliced,unspliced)
  #data=rbind(data,total)
  print(head(data))
  data$time=getField(data$cell,sep="_",3)
  if(summarise=="mean"){
    plotdata=data %>% group_by(.,individual,status,time) %>% summarise(.,mean_count=mean(count))
    ylab="mean_count"
  }else{
    plotdata=data %>% group_by(.,individual,status,time) %>% summarise(.,total_count=sum(count))
    ylab="total_count"
  }
  print(head(plotdata))
  if(return.data){
    return(plotdata)
  }
  ggplot(plotdata)+geom_point(aes(x=time,y=get(ylab),color=status))+
    geom_line(aes(x=time,y=get(ylab),group=status,color=status))+
    facet_wrap(~individual,scales = "free_y")+ylab(ylab)+
    scale_x_discrete(limits=CT_TIME_ORDER)+theme(axis.text.x = element_text(angle=45,hjust=1))
}

#gene_table columns: type, individual, feature
plotSplicingCountsByGene<-function(loom.matrix,meta,gene_table,summarise="mean",plot.all=FALSE,normalize=F,return.data=FALSE){
  feature=gene_table$feature[1]
  plot.data=NULL
  if(plot.all){
   for(i in unique(gene_table$individual)){
    for(j in unique(gene_table$type)){
      type.this=j
      print(type.this)
      individual.this=i
      print(individual.this)
      #cells.this=Cells(subset(srt,manual.level2==type.this&patient==individual.this))
      cells.this=rownames(dplyr::filter(meta,manual.level2==type.this,patient==individual.this))
      message(paste0(length(cells.this)," cells"))
      print(head(cells.this))
      expression.this=plotSingleSplicingCounts(loom.matrix,feature=feature,cells=cells.this,return.data=T,normalize = normalize)
      expression.this$cell.type=type.this
      if(is.null(plot.data)){
        plot.data=expression.this
      }else{
        plot.data=rbind(plot.data,expression.this)
      }
    }
   }
  }else{
    for(i in 1:nrow(gene_table)){
        type.this=gene_table[i,"type"]
        print(type.this)
        individual.this=gene_table[i,"individual"]
        print(individual.this)
        #cells.this=Cells(subset(srt,manual.level2==type.this&patient==individual.this))
        cells.this=rownames(dplyr::filter(meta,manual.level2==type.this,patient==individual.this))
        message(paste0(length(cells.this)," cells"))
        print(head(cells.this))
        expression.this=plotSingleSplicingCounts(loom.matrix,feature=feature,cells=cells.this,return.data=T,normalize=normalize)
        expression.this$cell.type=type.this
        if(is.null(plot.data)){
          plot.data=expression.this
        }else{
          plot.data=rbind(plot.data,expression.this)
        }
      }
  }
  plot.data.pseudo.point=plot.data[plot.data$time=="CT1",]
  plot.data.pseudo.point$time="CT25"
  plot.data=rbind(plot.data,plot.data.pseudo.point)
  print(head(plot.data))
  ggplot(plot.data)+geom_point(aes(x=time,y=mean_count,color=status))+
    geom_line(aes(x=time,y=mean_count,group=status,color=status))+
    facet_grid(cell.type~individual,scales = "free_y")+ylab("mean_count")+
    scale_x_discrete(limits=c(CT_TIME_ORDER,"CT25"))+theme(axis.text.x = element_text(angle=45,hjust=1))
}

#input: merged spliced data
plotSplicedVersusUnspliced<-function(loom.matrix,features=CIRCADIAN_GENES_MAIN,cells=NULL,normalize=F,summarise="mean",return.data=FALSE){
  if(is.null(cells)){
    cells=colnames(loom.matrix$spliced)
  }
  if(is.null(features)){
    features=rownames(loom.matrix$spliced)
  }
  if(features[!features %in% rownames(loom.matrix$spliced)] %>% length()>0){
    features=intersect(features,rownames(loom.matrix$spliced))
    emit_features=features[!features %in% rownames(loom.matrix$spliced)]
    message(paste0(length(emit_features), " gene(s) not in spliced data"))
  }
  spliced=loom.matrix$spliced[features,cells] %>% as.data.frame() %>% 
    rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  if(normalize){
    nCount_cell_spliced=loom.matrix$spliced[,cells] %>% t %>% as.data.frame() %>% rowSums()
    spliced$count=spliced$count*10000/nCount_cell_spliced
  }
  spliced$status="spliced"
  unspliced=loom.matrix$unspliced[features,cells] %>% as.data.frame() %>% 
    rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  if(normalize){
    nCount_cell_unspliced=loom.matrix$unspliced[,cells] %>% t %>% as.data.frame() %>% rowSums()
    unspliced$count=unspliced$count*10000/nCount_cell_unspliced
  }
  unspliced$status="unspliced"
  #ambiguous=loom.matrix$ambiguous[features,cells] %>% as.data.frame() %>% 
  #  rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  #ambiguous$status="ambiguous"
  spliced$time=getField(spliced$cell,sep="_",3)
  spliced$individual=getField(spliced$cell,sep="_",c(1,2))
  unspliced$time=getField(spliced$cell,sep="_",3)
  unspliced$individual=getField(unspliced$cell,sep="_",c(1,2))
  plotdata=NULL
  if(summarise=="mean"){
    spliced=spliced %>% group_by(.,feature,time) %>% summarise(.,spliced=mean(count))
    unspliced=unspliced %>% group_by(.,feature,time) %>% summarise(.,unspliced=mean(count))
    plotdata=left_join(spliced,unspliced,by=c("feature","time"))
    ylab="mean_count"
  }else{
    spliced=spliced %>% group_by(.,feature,time) %>% summarise(.,spliced=sum(count))
    unspliced=unspliced %>% group_by(.,feature,time) %>% summarise(.,unspliced=sum(count))
    plotdata=left_join(spliced,unspliced,by=c("feature","time"))
    ylab="total_count"
  }
  plotdata$time_numeric=gsub("CT","",plotdata$time) %>% as.numeric()
  print(head(plotdata))
  if(return.data){
    return(plotdata)
  }else{
    ggplot(plotdata)+geom_point(aes(x=spliced,y=unspliced,color=time),shape=1,size=5)+
      geom_text(aes(x=spliced,y=unspliced,label=time_numeric),size=3)+
      facet_wrap(~feature,scales = "free")+NoLegend()
  }
}

#input: merged spliced data
#the arg 'feature' should only provide 1 gene 
plotSingleSplicedVersusUnspliced<-function(loom.matrix,feature=CIRCADIAN_GENES_MAIN[1],cells=NULL,summarise="mean",return.data=FALSE){
  if(is.null(cells)){
    cells=colnames(loom.matrix$spliced)
  }
  spliced=loom.matrix$spliced[feature,cells] %>% t %>% as.data.frame() %>% 
    rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  spliced$status="spliced"
  unspliced=loom.matrix$unspliced[feature,cells] %>% t %>% as.data.frame() %>% 
    rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  unspliced$status="unspliced"
  ambiguous=loom.matrix$ambiguous[feature,cells] %>% t %>% as.data.frame() %>% 
    rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  ambiguous$status="ambiguous"
  spliced$time=getField(spliced$cell,sep="_",3)
  spliced$individual=getField(spliced$cell,sep="_",c(1,2))
  unspliced$time=getField(unspliced$cell,sep="_",3)
  unspliced$individual=getField(unspliced$cell,sep="_",c(1,2))
  plotdata=NULL
  if(summarise=="mean"){
    spliced=spliced %>% group_by(.,feature,time,individual) %>% summarise(.,spliced=mean(count))
    unspliced=unspliced %>% group_by(.,feature,time,individual) %>% summarise(.,unspliced=mean(count))
    plotdata=left_join(spliced,unspliced,by=c("feature","time","individual"))
    ylab="mean_count"
  }else{
    spliced=spliced %>% group_by(.,feature,time,individual) %>% summarise(.,spliced=sum(count))
    unspliced=unspliced %>% group_by(.,feature,time,individual) %>% summarise(.,unspliced=sum(count))
    plotdata=left_join(spliced,unspliced,by=c("feature","time","individual"))
    ylab="total_count"
  }
  plotdata$time_numeric=gsub("CT","",plotdata$time) %>% as.numeric()
  print(head(plotdata))
  if(return.data){
    return(plotdata)
  }else{
    ggplot(plotdata)+geom_point(aes(x=spliced,y=unspliced,color=time),shape=1,size=5)+
      geom_text(aes(x=spliced,y=unspliced,label=time_numeric),size=3)+
      facet_wrap(~individual,scales = "free")+NoLegend()
  }
}

#a clock like figure
plotPhaseofCellType<-function(df=NULL,p.cutoff=0.01,celltypes=NULL,colors,out.dir="../analysis/OSgenePredict/",method="JTK",embeding="box"){
  if(is.null(df)){
    plotdata=data.frame(CycID=character(0),PhaseShift=numeric(0),type=character(0))
  }else{
    plotdata=df
    if(method=="JTK"){
      colnames(plotdata)[c(5,7)]=c("PhaseShift","type")
    }else{
      colnames(plotdata)[c(2,12)]=c("PhaseShift","type")
    }
  }
  if(!is.null(celltypes)){
    for(type in celltypes){
      if(grepl("Temra/Teff",type)){
        type_char=gsub("Temra/Teff","Temra#Teff",type)
      }else{
        type_char=type
      }
      if(method=="LS"){
        file.path=paste0(out.dir,'LSresult_',type_char,".txt")
        p.str="p"
      }else if(method=="JTK"){
        file.path=paste0(out.dir,'JTKresult_',type_char,".txt")
        p.str="ADJ.P"
      }else{
        stop("Unknown method.")
      }
      if(file.exists(file.path)){
        outdata=read.delim(file.path)
        if(method=="JTK"){
          phase=dplyr::filter(outdata,get(p.str)<=p.cutoff)[c(1,5)]
          colnames(phase)[2]="PhaseShift"
        }else if(method=="LS"){
          phase=dplyr::filter(outdata,get(p.str)<=p.cutoff)[1:2]
        }else{
          stop("Unknown method.")
        }
        if(nrow(phase)==0){
          plotdata=rbind(plotdata,data.frame(CycID=NA,PhaseShift=0,type=type))
        }else{
          phase$type=type
          plotdata=rbind(plotdata,phase)
        }
      }else{
        plotdata=rbind(plotdata,data.frame(CycID=NA,PhaseShift=0,type=type))
      }
    }
  }
  print(head(plotdata))
  if(embeding=="box"){
    ggplot(plotdata)+geom_boxplot(aes(x=PhaseShift,y=type,fill=type),outlier.shape = NA,color=NA)+
      #geom_segment(aes(x=0,xend=0,y=0,yend=20),arrow=arrow(type="closed",angle=15))+
      scale_fill_manual(values=colors)+coord_polar()+scale_y_discrete(expand=c(0.1,0))+
      scale_x_continuous(n.breaks = 12,limits=c(0,24))+ylab("")+xlab("")+
      theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
            axis.text.x=element_text(size=11,vjust =-4,color="black"),
            panel.background=element_rect(fill = "white"),
            panel.grid.minor.x = element_line(color="grey",linewidth=1),
            plot.margin=margin(0,0,0,0,"cm"))
  }else if(embeding=="violin"){
    ggplot(plotdata)+geom_violin(aes(x=PhaseShift,y=type,fill=type))+
      geom_segment(aes(x=0,xend=0,y=0,yend=20),arrow=arrow(type="closed",angle=15))+
      scale_fill_manual(values=colors)+coord_polar()+scale_y_discrete(expand=c(0.1,0))+
      scale_x_continuous(n.breaks = 12,limits=c(0,24))+ylab("")+xlab("")+
      theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
            axis.text.x = element_text(size=12,vjust = -4),
            panel.background=element_rect(fill = "white"),
            panel.grid.minor.x = element_line(color="grey",linewidth=1))
  }else{
    stop("Unknown embeding.")
  }
}

#barplot showing number of oscillating genes of each cell type
plotCircadianGeneCount<-function(p.cutoff,celltypes,colors,out.dir="../analysis/OSgenePredict/",method="JTK"){
  plotdata=data.frame(type=character(0),count=numeric(0))
  for(type in celltypes){
    if(grepl("Temra/Teff",type)){
      type_char=gsub("Temra/Teff","Temra#Teff",type)
    }else{
      type_char=type
    }
    if(method=="JTK"){
      file.path=paste0(out.dir,'JTKresult_',type_char,".txt")
      p.col="ADJ.P"
    }else if(method=="LS"){
      file.path=paste0(out.dir,'LSresult_',type_char,".txt")
      p.col="p"
    }
    if(file.exists(file.path)){
      outdata=read.delim(file.path)
      count=dplyr::filter(outdata,get(p.col)<=p.cutoff) %>% nrow()
      plotdata.type=data.frame(type=type,count=count)
      plotdata=rbind(plotdata,plotdata.type)
    }else{
      plotdata.type=data.frame(type=type,count=0)
      plotdata=rbind(plotdata,plotdata.type)
    }
  }
  ggplot(plotdata)+geom_bar(aes(x=type,y=count,fill=type),stat="identity",color="black")+
    scale_fill_manual(values=colors)+scale_y_continuous(expand=c(0,0))+scale_x_discrete(limits=rev(celltypes))+
    xlab("number of circadian genes")+ylab("")+NoLegend()+coord_flip()
}


plotGeneExpressionHeatmapByCT<-function(srt,cell.type,gene.list=CIRCADIAN_GENES_MAIN,return.data=F){
  cell.list=rownames(srt@meta.data[srt@meta.data$type==cell.type,])
  data=t(LayerData(srt,layer = "count",cells=cell.list,features=gene.list)) %>% as.data.frame()
  data=rownames_to_column(data,var="cell_id")
  data=gather(data,key="genes",value="expression",-cell_id)
  meta=srt@meta.data[data$cell_id %>% unique(),c("cell_id","CT")]
  data=left_join(data,meta,by="cell_id")
  print(head(data))
  data=downsample(data,"CT")
  data=data %>% group_by(genes,CT) %>% summarise(.,expression_all=log(sum(expression))+1)
  plotdata=data %>% group_by(genes) %>% reframe(relative_expression=scale(expression_all))
  plotdata$CT=data$CT
  #plotdata=data
  if(return.data){
    return(plotdata)
  }else{
    print(head(plotdata))
  }
  
  ggplot(plotdata)+geom_bin_2d(aes(x=genes,y=CT,fill=relative_expression))+
   scale_fill_gradientn(colors = rev(brewer.pal(11,"RdBu")))+
    scale_x_discrete(limits=gene.list)+scale_y_discrete(limits=rev(CT_TIME_ORDER))+
    ylab("")+xlab("")+theme(axis.text.x=element_blank(),axis.ticks=element_blank(),legend.position = "top")
}

#violin plot showing expression of genes by celltype
#arg1 is integrated and annotated srt, celltype data is in manual.level2
plotExpressionByCelltype<-function(int.srt,gene.list=CIRCADIAN_GENES_MAIN,colors=generateColor(25)){
  data=FetchData(int.srt,vars=gene.list)
  data=rownames_to_column(data,var="observation")
  data=gather(data,key=feature,value=count,-observation)
  meta=int.srt@meta.data[data$observation,c("patient","manual.level2")]
  meta=rownames_to_column(meta,var="observation")
  data=left_join(data,meta,by="observation")
  data=dplyr::filter(data,count>0)
  print(head(data))
  x.labels=data$manual.level2 %>% getField(.,"_",1) %>% unique() %>% sort()
  print(x.labels)
  ggplot(data)+geom_violin(aes(x=manual.level2,y=count,fill=manual.level2))+scale_fill_manual(values=colors)+
    scale_x_discrete(labels=x.labels)+facet_wrap(~feature)+
    theme(axis.ticks.x=element_blank(),axis.text.x=element_text(angle=60,hjust=1),
    panel.background=element_rect(fill="white",color="black"),panel.grid = element_line(color="grey"),
    legend.key.size=unit(10,"pt"))+labs(fill="cell_type",x="")+
    guides(fill=guide_legend(ncol = 1))
}

# plotCellProliferationByCTByIndividual
plotCellProliferationByCTByIndividual<-function(srt,meta.data.col="Proliferation_Score1",return.data=F){
  meta.data=srt@meta.data[,c("patient","CT","manual_NI",meta.data.col)]
  meta.data=meta.data[meta.data$manual_NI=="TNK_proliferatig_MKI67",]
  meta.data=meta.data %>% group_by(patient,CT) %>% summarise(median_proliferation_score=median(.data[[meta.data.col]]))
  print(head(meta.data))
  if(return.data){
    return(meta.data)
  }
  ggplot(meta.data)+geom_point(aes(x=CT,y=median_proliferation_score,color=patient))+
    geom_line(aes(x=CT,y=median_proliferation_score,group=patient))+
    scale_x_discrete(limits=CT_TIME_ORDER)
}

#arg1 is the return of "mergeSamplesOneDay"
#the ggplot's function geom_smooth change the model when sample size increase automatically, if we specify
#a fixed model arbitrarily, the parameter k will change, as a result, the fitting line is not comparable and
#the result is not reliable, this function is not so useful.
plotGeneExpressionByCT<-function(merged.srt,gene.list=CIRCADIAN_GENES_MAIN,cell.type=NULL,k=10,alpha=0.25,
                                 CT_field=3,patient_filed=1,return_data=F,use.mean=T,layer="data"){
  data=fetchMergedDataOneCellType(merged.srt,gene.list=gene.list,cell.type=cell.type,CT_field=CT_field,patient_filed=patient_filed,layer=layer)
  if(return_data){
    return(data)
  }
  data_summar=data %>% group_by(.,features,time_numeric) %>% summarise(.,mean=mean(values))
  plot.tittle=paste0(cell.type," cells")
  if(use.mean){
    ggplot(data)+geom_point(aes(x=time_numeric,y=values,group=patient,color=patient),alpha=alpha)+
      geom_path(data=data_summar,aes(x=time_numeric,y=mean),color="red")+
      facet_wrap(vars(features),nrow=ceiling(length(gene.list)/6),scales="free_y")+xlab("")+ylab("")+
      ggtitle(plot.tittle)
  }else{
    if(max(data[,c("features","time")] %>% table())>100){
      ggplot(data,aes(x=time_numeric,y=values,group=patient,color=patient))+
        geom_smooth(method="gam",formula=y~s(x,bs="cs",k=k))+
        geom_point(alpha=alpha)+facet_wrap(vars(features),nrow=ceiling(length(gene.list)/6),scales = "free_y")+xlab("")+ylab("")+
        ggtitle(plot.tittle)
    }else{
      ggplot(data,aes(x=time_numeric,y=values,group=patient,color=patient))+
        geom_smooth()+
        geom_point(alpha=alpha)+facet_wrap(vars(features),nrow=ceiling(length(gene.list)/6),scales = "free_y")+xlab("")+ylab("")+
        ggtitle(plot.tittle)
    }
  }
}

viewRNAVelosity<-function(loom.matrix,min.transcripts.per.cell,min.cells.per.gene,deltaT=1,return=F){
  test=loom.matrix
  emat=test$spliced
  emat=emat[,colSums(emat)>=1e3]
  emat=emat[!duplicated(rownames(emat)),]
  r=Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T,min.transcripts.per.cell=min.transcripts.per.cell,min.cells.per.gene=min.cells.per.gene)
  r$adjustVariance(plot=T,do.par=T,gam.k=10)
  r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
  r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine')
  r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
  r$getEmbedding(type='PCA',embeddingType='UMAP',perplexity=50,verbose=T)
  
  emat=test$spliced
  nmat=test$unspliced
  emat=emat[,rownames(r$counts)]
  nmat=nmat[,rownames(r$counts)] # restrict to cells that passed p2 filter
  # take cluster labels
  cluster.label=r$clusters$PCA[[1]]
  print(head(cluster.label))
  values_label=names(cluster.label) %>% strsplit(.,"_") %>% lapply(.,`[`,3) %>% unlist()%>% as.vector()
  print(head(values_label))
  names(values_label)=names(cluster.label)
  values_label=factor(values_label,levels = CT_TIME)
  cell.colors=sccore::fac2col(values_label)
  # take embedding
  emb=r$embeddings$PCA$UMAP
  cell.dist=as.dist(1-armaCor(t(r$reductions$PCA)))
  emat=filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
  nmat=filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
  shared=intersect(rownames(emat),rownames(nmat))
  shared.cell=intersect(colnames(emat),colnames(nmat))
  fit.quantile=0.02
  rvel.cd=gene.relative.velocity.estimates(emat[shared,shared.cell],nmat[shared,shared.cell],deltaT=deltaT,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)
  show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),
                                 cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,
                                 arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
  if(return){
    return(rvel.cd)
  }
}

#input are a list of seurat object
plotCellRatioVersusCT<-function(srts,colors=pal_d3("category20")(20),type.col.name="type"){
  plot.data=data.frame("type"=character(0),"CT"=character(0))
  for(srt in srts){
    new.meta=getMetaData(srt,c(type.col.name,"CT"))
    colnames(new.meta)=c("type","CT")
    plot.data=rbind(plot.data,new.meta)
  }
  ddplot(plot.data,reverse.col=T,x.order=CT_TIME_ORDER)+scale_fill_manual(values=colors)
}

#input is merged or integrated srt
plotClockByGenes<-function(srt,gene.list=CIRCADIAN_GENES_MAIN,cell.type=NULL,CT_field=3,patient_filed=c(1,2),summa="median"){
  data=fetchMergedDataOneCellType(srt,gene.list=gene.list,cell.type=cell.type,CT_field=CT_field,patient_filed=patient_filed)
  print(head(data))
  if(summa=="median"){
    data_summar=data %>% group_by(.,features,time_numeric,patient) %>% summarise(.,summa=median(values))
  }else if(summa=="sum"){
    data_summar=data %>% group_by(.,features,time_numeric,patient) %>% summarise(.,summa=sum(values))
  }else{
    stop("Invalid summa")
  }
  print(head(data_summar))
  ggplot(data)+geom_point(aes(x=time_numeric,y=values),alpha=0.1)+
      geom_path(data=data_summar,aes(x=time_numeric,y=summa,group=patient,color=patient))+
      facet_wrap(vars(features),nrow=ceiling(length(gene.list)/6),scales = "free_y")+xlab("")+ylab("")+NoLegend()+
      ggtitle(cell.type)
}

plotRhythmicityPvalue<-function(celltypes,p.cutoff=1,features=CIRCADIAN_GENES_MAIN,colors,out.dir="../analysis/OSgenePredict/",method="JTK"){
  plotdata=NULL
  color_bool=rep(T,length(colors))
  index=0
  for(type in celltypes){
    #print(type)
    index=index+1
    if(grepl("Temra/Teff",type)){
      type_char=gsub("Temra/Teff","Temra#Teff",type)
    }else{
      type_char=type
    }
    if(method=="JTK"){
      file.path=paste0(out.dir,'JTKresult_',type_char,".txt")
      p.col="ADJ.P"
    }else if(method=="LS"){
      file.path=paste0(out.dir,'LSresult_',type_char,".txt")
      p.col="p"
    }
    if(file.exists(file.path)){
      outdata=read.delim(file.path)
      if(nrow(outdata)==0){
        color_bool[index]=F
        next
      }
      plotdata.type=dplyr::filter(outdata,get(p.col)<=p.cutoff,CycID %in% features)
      #print(head(plotdata.type))
      plotdata.type$cell_type=type
      if(is.null(plotdata)){
        plotdata=plotdata.type
      }else{
        plotdata=rbind(plotdata,plotdata.type)
      }
    }else{
      color_bool[index]=F
      next
    }
  }
  print(head(plotdata))
  print(color_bool)
  ggplot(plotdata)+geom_point(aes(x=CycID,y=cell_type,size=-log2(get(p.col)),fill=cell_type),shape=21)+
    xlab("")+ylab("")+scale_y_discrete(limits=rev(celltypes))+
    scale_fill_manual(values=colors[color_bool])+
    theme(axis.text.x=element_text(angle=45,hjust=1),
          panel.background=element_rect(fill="white",colour = "black"),
          panel.grid=element_line(color="grey"))+
   guides(size=guide_legend(title="-log2(p)",override.aes=list(fill="black")),fill=FALSE)
}

plotMetaCellByIndividual<-function(srt,cell.type,feature,layer="counts",norm.dist=T){
  data=fetchMergedDataOneCellType(srt,cell.type=cell.type,gene.list=feature,CT_field=3,patient_filed = c(1,2),layer=layer,filter_data=F)
  control=srt[["nCount_RNA"]] %>% as.data.frame()
  control=rownames_to_column(control,var = "observations")
  colnames(control)[2]="control_expression"
  print(head(data))
  print(head(control))
  data=left_join(data,control)
  #data$normalized=data$values/data$control_expression
  data$normalized=1000000*data$values/data$control_expression
  data=data %>% group_by(time,patient) %>% mutate(mean_expression=mean(log(normalized+1))) %>% group_by(patient) %>%
    mutate(relative_expression=(mean_expression-mean(mean_expression))/sd(mean_expression))
  data=data %>% group_by(time) %>% mutate(mean_relative=mean(relative_expression),sd_relative=sd(relative_expression))
  data=data %>% mutate(relative_q25=quantile(relative_expression,na.rm=T)[2],
                       relative_q75=quantile(relative_expression,na.rm=T)[4],
                       relative_median=median(relative_expression))
  print(head(data))
  if(layer=="data"){
    ggplot(data)+geom_boxplot(aes(x=time,y=values))+
      geom_jitter(aes(x=time,y=values))+facet_wrap(~patient,scale="free_y")+
      ggtitle(paste0(cell.type,": ",feature))
  }else{
    ggplot(data)+geom_boxplot(aes(x=time,y=normalized))+
      geom_jitter(aes(x=time,y=normalized))+facet_wrap(~patient,scale="free_y")+
      ggtitle(paste0(cell.type,": ",feature))
  }

  #if(norm.dist){
  #  ggplot(data)+geom_point(aes(x=time,y=mean_relative))+
  #    geom_point(aes(x=time,y=relative_expression))+
  #    geom_line(aes(x=time,y=mean_relative,group=features))+
  #    geom_errorbar(aes(x=time,ymin=mean_relative-sd_relative,ymax=mean_relative+sd_relative),width=0.2)+
  #    ggtitle(paste0(cell.type,": ",feature))
  #}else{
  #  ggplot(data)+geom_point(aes(x=time,y=relative_median))+
  #    geom_line(aes(x=time,y=relative_median,group=features))+
  #    geom_errorbar(aes(x=time,ymin=relative_q25,ymax=relative_q75),width=0.2)+
  #    ggtitle(paste0(cell.type,": ",feature))
  #}
}

plotPseudobulkByCTByIndividual<-function(srt,cell.type,features,individuals=NULL,normalize.data=F,time.points=CT_TIME_ORDER,
                                         return.data=F,proportion=1,split.by="no",layer="count",relative=T,log=F,plot="raw"){
  data=fetchMergedDataOneCellType(srt,cell.type=cell.type,gene.list=features,CT_field=3,patient_filed = c(1,2),layer=layer,filter_data=F)
  print(head(data))
  if(is.null(individuals)){
    individuals=unique(data$patient)
  }
  plotdata=NULL
  for(individual in individuals){
    message("looping ..")
    message(individual)
    #pmeta=srt@meta.data[srt@meta.data$patient==individual,c("nCount_RNA","CT")] %>% downsample(.,"CT") %>% 
    #  group_by(CT) %>% summarise(sum=sum(nCount_RNA))
    #pmeta$min=min(pmeta$sum)
    #pmeta$normalize_cof=pmeta$sum/pmeta$min
    #colnames(pmeta)[1]="time"
    for(feature in features){
      block_data=dplyr::filter(data,features==feature,patient==individual)
      block_data_min_cell_count=min(table(block_data$time))
      block_data=rbind(block_data,data.frame(observations=rep("pseudo",block_data_min_cell_count*proportion),
                                             features=rep(feature,block_data_min_cell_count*proportion),
                                             values=rep(0,block_data_min_cell_count*proportion),
                                             time=rep("CTpseudo",block_data_min_cell_count*proportion),
                                             patient=rep("pseudo",block_data_min_cell_count*proportion),
                                             time_numeric=rep("pseudo",block_data_min_cell_count*proportion),
                                             cell_type=rep(cell.type,block_data_min_cell_count*proportion)))
      block_data=downsample(block_data,"time")
      internal_control_counts=NULL
      if(normalize.data){
        cells_block=block_data$observations
        cells_block=cells_block[cells_block!="pseudo"]
        internal_control_counts=LayerData(srt,layer=layer,cells=cells_block,features=INTERNAL_CONTROL)
        internal_control_counts=internal_control_counts %>% as.data.frame() %>% t %>% as.data.frame()
        internal_control_counts=rowSums(internal_control_counts) %>% as.data.frame()
        colnames(internal_control_counts)[1]<-"count"
        internal_control_counts$time<-getField(rownames(internal_control_counts),sep = "_",3)
        internal_control_counts=internal_control_counts %>% group_by(time) %>% summarise(.,sum_control=sum(count))
      }
      block_data=block_data %>% group_by(.,time) %>% summarise(.,count=sum(values))
      block_data=block_data[block_data$time!="CTpseudo",]
      if(normalize.data){
        block_data=left_join(block_data,internal_control_counts,by="time")
        block_data$coef=block_data$sum_control/min(block_data$sum_control)
        block_data$normalized_count=block_data$count/block_data$coef
      }
      #print(block_data)
      block_data$feature=feature
      block_data$individual=individual
      #block_data=left_join(block_data,pmeta,by="time")
      #block_data$raw_count=block_data$count
      #block_data$count=block_data$count/block_data$normalize_cof
      block_data$cell_count=block_data_min_cell_count
      if(is.null(plotdata)){
        plotdata=block_data
      }else{
        plotdata=rbind(plotdata,block_data)
      }
    }
  }
  if(relative){
    if(normalize.data){
      #plotdata$relative_expression=(plotdata %>% group_by(feature,individual) %>% reframe(relative_expression=scale(normalized_count)))$relative_expression[,1]
      if(log){
        plotdata=plotdata %>% group_by(feature,individual) %>% mutate(relative_expression=(log2(normalized_count+1)-mean(log2(normalized_count+1)))/sd(log2(normalized_count+1)))
      }else{
        plotdata=plotdata %>% group_by(feature,individual) %>% mutate(relative_expression=(normalized_count-mean(normalized_count))/sd(normalized_count))
      }
    }else{
      #plotdata$relative_expression=(plotdata %>% group_by(feature,individual) %>% reframe(relative_expression=scale(count)))$relative_expression[,1]
      if(log){
        plotdata=plotdata %>% group_by(feature,individual) %>% mutate(relative_expression=(log2(count+1)-mean(log2(count+1)))/sd(log2(count+1)))
      }else{
        plotdata=plotdata %>% group_by(feature,individual) %>% mutate(relative_expression=(count-mean(count))/sd(count))
      }
    }
  }else{
    plotdata$relative_expression=plotdata$normalized_count
  }
  plotdata$time_numeric=plotdata$time %>% gsub("CT","",.) %>% as.numeric()
  plotdata=plotdata %>% group_by(time) %>% mutate(CT_mean=mean(relative_expression),CT_sd=sd(relative_expression),
                                                  CT_median=median(relative_expression),CT_q25=quantile(relative_expression,na.rm=T)[2],
                                                  CT_q75=quantile(relative_expression,na.rm=T)[4])
  if(return.data){
    return(plotdata)
  }
  print(head(plotdata))
  if(split.by=="individual"){
      ggplot(plotdata)+geom_line(aes(x=time,y=relative_expression,group=feature))+
        geom_point(aes(x=time,y=relative_expression,size=count,color=feature))+
        scale_x_discrete(limits=time.points)+scale_color_manual(values=generateColor(length(features)))+
        facet_wrap(~individual,scales="free_y")+guides(size=guide_legend(title="UMI count"))+ggtitle(paste0(cell.type))
  }else{
    if(plot=="raw"){
      ggplot(plotdata)+geom_line(aes(x=time,y=relative_expression,group=individual))+
        geom_point(aes(x=time,y=relative_expression,size=count,color=feature),alpha=0.5)+
        scale_x_discrete(limits=time.points)+scale_color_manual(values=generateColor(length(features)))+
        facet_wrap(~feature,scales="free_y")+guides(size=guide_legend(title="UMI count"))+ggtitle(paste0(cell.type))+ylab("z-score")
    }else if(plot=="errorbar"){
      if(log){
        ggplot(plotdata)+geom_line(aes(x=time,y=CT_mean,group=feature))+
          geom_point(aes(x=time,y=CT_mean))+
          geom_errorbar(aes(x=time,ymin=CT_mean-CT_sd,ymax=CT_mean+CT_sd),width=0.2)+
          scale_x_discrete(limits=time.points)+facet_wrap(~feature,scales="free_y")+
          ggtitle(paste0(cell.type))+ylab("z-score")
      }else{
        ggplot(plotdata)+geom_line(aes(x=time,y=CT_median,group=feature))+
          geom_point(aes(x=time,y=CT_median))+
          geom_errorbar(aes(x=time,ymin=CT_q25,ymax=CT_q75),width=0.2)+
          scale_x_discrete(limits=time.points)+facet_wrap(~feature,scales="free_y")+
          ggtitle(paste0(cell.type))+ylab("z-score")
      }

    }else if(plot=="box"){
      ggplot(plotdata)+
        geom_boxplot(aes(x=time,y=relative_expression))+
        scale_x_discrete(limits=time.points)+facet_wrap(~feature,scales="free_y")+
        ggtitle(paste0(cell.type))+ylab("z-score")
    }else{
      message("unsupported plot type, supporting: raw/errorbar")
    }
  }
}

plotBatchByCT<-function(srt,cell.type,this.patient=PATIENTS[1],slot="counts"){
  data=fetchMergedDataOneCellType(srt,cell.type=cell.type,gene.list=c("PTEN","ARNTL","PER1","CRY1","DNAJA1"),
        CT_field=3,patient_filed = c(1,2),layer = slot,filter_data=F)
  data.meta=srt@meta.data[,c("nCount_RNA","nFeature_RNA")]
  data.meta=rownames_to_column(data.meta,var="observations")
  data=left_join(data,data.meta,by="observations")
  #print(head(data))
  plotdata.typecellcount=data %>% group_by(.,time,patient) %>% summarise(.,type.cell.count=n())
  plotdata.type.UMIcount=data %>% group_by(.,time,patient) %>% summarise(.,type.umi.count=sum(nCount_RNA))
  plotdata.gene.UMIcount=data %>% group_by(.,time,patient,features) %>% summarise(.,count=sum(values))
  plotdata.gene.UMIcount=spread(plotdata.gene.UMIcount,key = features,value=count)
  print(head(plotdata.gene.UMIcount))
  
  
  plotdata.totalcellcount=srt@meta.data %>% group_by(.,CT,patient) %>% summarise(.,total.cell.count=n())
  plotdata.total.UMIcount=srt@meta.data %>% group_by(.,CT,patient) %>% summarise(.,total.umi.count=sum(nCount_RNA))
  colnames(plotdata.total.UMIcount)[which(colnames(plotdata.total.UMIcount)=="CT")]="time"
  colnames(plotdata.totalcellcount)[which(colnames(plotdata.totalcellcount)=="CT")]="time"
  
  plotdata=left_join(plotdata.typecellcount,plotdata.type.UMIcount,by=c("time","patient"))
  #plotdata=left_join(plotdata,plotdata.totalcellcount,by=c("time","patient"))
  #plotdata=left_join(plotdata,plotdata.total.UMIcount,by=c("time","patient"))
  plotdata=left_join(plotdata,plotdata.gene.UMIcount,by=c("time","patient"))
  #print(head(plotdata))
  #cell.count.plot=ggplot(plotdata)+geom_point(aes(x=time,y=type.cell.count))+
  #  geom_line(aes(x=time,y=type.cell.count,group=patient),stat = "identity")+
  #  scale_x_discrete(limits=CT_TIME_ORDER)+
  #  facet_wrap(~patient,scales = "free_y")+
  #  theme(axis.text.x=element_text(angle=60,hjust=1))
  
  #type.umi.plot=ggplot(plotdata)+geom_point(aes(x=time,y=type.umi.count/type.cell.count))+
  #  geom_line(aes(x=time,y=type.umi.count/type.cell.count,group=patient),stat = "identity")+
  #  scale_x_discrete(limits=CT_TIME_ORDER)+
  #  facet_wrap(~patient,scales = "free_y")+
  #  theme(axis.text.x=element_text(angle=60,hjust=1))
  
  #total.cell.plot=ggplot(plotdata)+geom_point(aes(x=time,y=total.cell.count))+
  #  geom_line(aes(x=time,y=total.cell.count,group=patient),stat = "identity")+
  #  scale_x_discrete(limits=CT_TIME_ORDER)+
  #  facet_wrap(~patient,scales="free_y")+
  #  theme(axis.text.x=element_text(angle=60,hjust=1))
  
  #total.umi.plot=ggplot(plotdata)+geom_point(aes(x=time,y=total.umi.count/total.cell.count))+
  #  geom_line(aes(x=time,y=total.umi.count/total.cell.count,group=patient),stat = "identity")+
  #  scale_x_discrete(limits=CT_TIME_ORDER)+
  #  facet_wrap(~patient,scales = "free_y")+
  #  theme(axis.text.x=element_text(angle=60,hjust=1))
  
  #ggarrange(cell.count.plot,type.umi.plot,total.cell.plot,total.umi.plot)
  
  #plots=list()
  #i=1
  #for(this.patient in unique(plotdata$patient)){
  #  sub.plotdata=dplyr::filter(plotdata,patient==this.patient)[3:ncol(plotdata)]
  #  corr=round(cor(sub.plotdata), 2)
  #  p.mat=cor_pmat(sub.plotdata)
  #  plots[[i]]=ggcorrplot(corr, type = "lower", lab = TRUE)
  #  i=i+1
  #}
  #ggarrange(plotlist = plots,ncol=4,nrow=3)
  plots=list()
  sub.plotdata=dplyr::filter(plotdata,patient==this.patient)[3:ncol(plotdata)]
  sub.plotdata2=sub.plotdata
  corr=round(cor(sub.plotdata),2)
  plots[[1]]=ggcorrplot(corr, type = "lower", lab = TRUE,title = "No normalization")
  
  sub.plotdata$cell.norm.coefficient=sub.plotdata$type.cell.count/min(sub.plotdata$type.cell.count)
  sub.plotdata$type.umi.count=sub.plotdata$type.umi.count/sub.plotdata$cell.norm.coefficient
  sub.plotdata$PTEN=sub.plotdata$PTEN/sub.plotdata$cell.norm.coefficient
  sub.plotdata$ARNTL=sub.plotdata$ARNTL/sub.plotdata$cell.norm.coefficient
  sub.plotdata$PER1=sub.plotdata$PER1/sub.plotdata$cell.norm.coefficient
  sub.plotdata$CRY1=sub.plotdata$CRY1/sub.plotdata$cell.norm.coefficient
  sub.plotdata$DNAJA1=sub.plotdata$DNAJA1/sub.plotdata$cell.norm.coefficient
  sub.plotdata$cell.norm.coefficient=NULL
  corr=round(cor(sub.plotdata),2)
  
  plots[[2]]=ggcorrplot(corr, type = "lower", lab = TRUE,title="Normalized by num of cell")
  
  sub.plotdata2$umi.norm.coefficient=sub.plotdata2$type.umi.count/min(sub.plotdata2$type.umi.count)
  sub.plotdata2$PTEN=sub.plotdata2$PTEN/sub.plotdata2$umi.norm.coefficient
  sub.plotdata2$ARNTL=sub.plotdata2$ARNTL/sub.plotdata2$umi.norm.coefficient
  sub.plotdata2$PER1=sub.plotdata2$PER1/sub.plotdata2$umi.norm.coefficient
  sub.plotdata2$CRY1=sub.plotdata2$CRY1/sub.plotdata2$umi.norm.coefficient
  sub.plotdata2$DNAJA1=sub.plotdata2$DNAJA1/sub.plotdata2$umi.norm.coefficient
  sub.plotdata2$umi.norm.coefficient=NULL
  corr=round(cor(sub.plotdata2),2)
  plots[[3]]=ggcorrplot(corr,type = "lower",lab = TRUE,title="Normalized by umi")
  
  sub.plotdata$umi.norm.coefficient=sub.plotdata$type.umi.count/min(sub.plotdata$type.umi.count)
  sub.plotdata$PTEN=sub.plotdata$PTEN/sub.plotdata$umi.norm.coefficient
  sub.plotdata$ARNTL=sub.plotdata$ARNTL/sub.plotdata$umi.norm.coefficient
  sub.plotdata$PER1=sub.plotdata$PER1/sub.plotdata$umi.norm.coefficient
  sub.plotdata$CRY1=sub.plotdata$CRY1/sub.plotdata$umi.norm.coefficient
  sub.plotdata$DNAJA1=sub.plotdata$DNAJA1/sub.plotdata$umi.norm.coefficient
  sub.plotdata$umi.norm.coefficient=NULL
  corr=round(cor(sub.plotdata),2)
  
  #print(head(sub.plotdata))
  
  plots[[4]]=ggcorrplot(corr,type = "lower",lab = TRUE,title="Normalized by num of cell and umi")
  ggarrange(plotlist = plots)
}

plotPseudobulkByCT<-function(srt,cell.type,features=CIRCADIAN_GENES_MAIN,return.data=F,rep=3,proportion=0.3,layer="count",normalize.data=F){
  data=fetchMergedDataOneCellType(srt,cell.type=cell.type,gene.list=features,CT_field=3,patient_filed = c(1,2),layer = layer,filter_data=F)
  print(head(data))
  #patients=data$patient %>% table() %>% names()
  plotdata=NULL
  downsample_idx=1
  for(this.rep in 1:rep){
    message("looping ..")
    message(paste0(downsample_idx,"/",rep))
    for(feature in features){
      #message(patient)
      block_data=dplyr::filter(data,features==feature)
      block_data_min_cell_count=min(table(block_data$time))
      #message(block_data_min_cell_count)
      block_data=rbind(block_data,data.frame(observations=rep("pseudo",block_data_min_cell_count*proportion),
                                             features=rep(feature,block_data_min_cell_count*proportion),
                                             values=rep(0,block_data_min_cell_count*proportion),
                                             time=rep("CTpseudo",block_data_min_cell_count*proportion),
                                             patient=rep("pseudo",block_data_min_cell_count*proportion),
                                             time_numeric=rep("pseudo",block_data_min_cell_count*proportion),
                                             cell_type=rep(cell.type,block_data_min_cell_count*proportion)))
      block_data=downsample(block_data,"time")
      #print(head(block_data))
      internal_control_counts=NULL
      if(normalize.data){
        cells_block=block_data$observations
        cells_block=cells_block[cells_block!="pseudo"]
        internal_control_counts=LayerData(srt,layer=layer,cells=cells_block,features=INTERNAL_CONTROL)
        internal_control_counts=internal_control_counts %>% as.data.frame() %>% t %>% as.data.frame()
        internal_control_counts=rowSums(internal_control_counts) %>% as.data.frame()
        colnames(internal_control_counts)[1]<-"count"
        internal_control_counts$time<-getField(rownames(internal_control_counts),sep = "_",3)
        internal_control_counts=internal_control_counts %>% group_by(time) %>% summarise(.,sum_control=sum(count))
      }
      block_data=block_data %>% group_by(.,time) %>% summarise(.,count=sum(values))
      block_data=block_data[block_data$time!="CTpseudo",]
      if(normalize.data){
        block_data=left_join(block_data,internal_control_counts,by="time")
        block_data$coef=block_data$sum_control/min(block_data$sum_control)
        block_data$normalized_count=block_data$count/block_data$coef
      }
      print(block_data)
      block_data$downsample=downsample_idx
      block_data$feature=feature
      if(is.null(plotdata)){
        plotdata=block_data
      }else{
        plotdata=rbind(plotdata,block_data)
      }
    }
    downsample_idx=downsample_idx+1
  }
  
  if(return.data){
    return(plotdata)
  }
  plotdata$downsample<-as.character(plotdata$downsample)
  plotdata$time_numeric=plotdata$time %>% gsub("CT","",.) %>% as.numeric()
  if(normalize.data){
    #normalize_data=srt@meta.data[c("nCount_RNA","CT")] %>% group_by(.,CT) %>% summarise(sum=sum(nCount_RNA))
    #normalize_data$sum=normalize_data$sum/min(normalize_data$sum)
    #colnames(normalize_data)=c("time","sum")
    #plotdata=left_join(plotdata,normalize_data,by="time")
    #plotdata$normalized_count=plotdata$count/plotdata$sum
    #print(head(plotdata))
    ggplot(plotdata)+geom_smooth(aes(x=time_numeric,y=normalized_count),span=1)+
      geom_point(aes(x=time_numeric,y=normalized_count))+
      facet_wrap(~feature,scales="free_y")+
      #scale_x_discrete(limits=CT_TIME_ORDER)+ggtitle(cell.type)+
      theme(axis.text.x=element_text(angle=45,hjust=1),
            panel.background=element_rect(fill="white",colour = "black"),
            panel.grid = element_line(color="grey"))+
      NoLegend()+xlab("")+ggtitle(cell.type)
  }else{
    print(head(plotdata))
    ggplot(plotdata)+geom_smooth(aes(x=time_numeric,y=count),span=1)+
      geom_point(aes(x=time_numeric,y=count))+
      facet_wrap(~feature,scales="free_y")+
      #scale_x_discrete(limits=CT_TIME_ORDER)+ggtitle(cell.type)+
      theme(axis.text.x=element_text(angle=45,hjust=1),
            panel.background=element_rect(fill="white",colour = "black"),
            panel.grid = element_line(color="grey"))+
      NoLegend()+xlab("")+ggtitle(cell.type)
  }
  #plotdata$normalized_count=plotdata$normalized_count/(1+min(plotdata$normalized_count))
}

plotVenn<-function(set1,set2,set.names=c("set1","set2"),colors=c("red","blue")){
  VennDiagram::venn.diagram(list(set1,set2),category.names = set.names,
  filename = NULL,disable.logging = T,cat.pos=c(-30,30),cat.dist=0.03,fill=colors) %>% grid.draw()
}


#::::::::::::::::::::#
# AUXILLIARY FUNCTION
#::::::::::::::::::::#
#check if each individual have 6 loom files
#return a vector of individual ids that have 6 loom files
checkLoomFilesForIndividuals<-function(){
  for(individual in PATIENTS){
    for(i in 1:6){
      sample_id=paste0(individual,"_",i)
      sample_loom_path=paste0(sample_id,"/outs/per_sample_outs/",sample_id,"/count/velocyto/")
      sample_loom_file_name=list.files(sample_loom_path)[1]
      sample_loom_full_path=paste0(sample_loom_path,sample_loom_file_name)
      message("loom file: ")
      message(sample_loom_full_path)
    }
  }
}



predictCTofCells<-function(srt,features,cell.type){
  lm.test.data=FetchData(subset(srt,type==cell.type),vars=features)
  for(feature in features){
    feature.data=lm.test.data[,feature]
    feature.data=dplyr::filter(feature.data,get(feature)>0)
    feature.data$CT=rownames(feature.data) %>% strsplit(.,"_") %>% lapply(.,`[`,3) %>% unlist() %>% gsub("CT","",.) %>% as.numeric()
    colnames(feature.data)=c("expression","CT")
    lm.model=geneIsOscillated(feature.data)
  }
}

#get field
#vec: vector of string
#sep: separator
#field: wanted field
getField<-function(vec,sep,field){
  if(length(field)>1){
    strsplit(vec,sep,fixed=T) %>% lapply(.,`[`,field) %>% lapply(.,paste0,collapse=sep) %>% unlist()
  }else{
    strsplit(vec,sep,fixed=T) %>% lapply(.,`[`,field) %>% unlist()
  }
}


#test oscillating
#input is data.frame, colnames are CT(numeric), gene_names, expression
#return a list
geneIsOscillated<-function(data){
  res=list()
  newdata=data.frame("CT"=1:24)
  c_model=lm(expression~sin(CT*pi/12)+cos(CT*pi/12),data)
  newdata$predict=predict(c_model,newdata)
  rawdata=ggplot(data)+geom_point(aes(x=CT,y=expression),alpha=0.1)+geom_line(data=newdata,aes(x=CT,y=predict),color="blue")
  residual=ggplot(c_model$residuals %>% as.data.frame())+geom_density(aes(x=.))
  res$model=c_model
  res$residual=residual
  res$rawdata=rawdata
  res$normality_of_residual=shapiro.test(c_model$residuals)
  return(res)
}

#this function not work
convertMouseGeneList<-function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x ,
                   mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

#input are from "mergeSamplesOneDay"
#return a data frame, colnames are observations, features, values, time (CT), cell_type
fetchMergedDataOneCellType<-function(merged.srt,gene.list=CIRCADIAN_GENES_MAIN,cell.type=NULL,CT_field=1,patient_filed=NULL,layer="data",filter_data=T){
  data=LayerData(merged.srt,layer=layer,cells=rownames(merged.srt@meta.data[merged.srt@meta.data$type==cell.type,]),features=gene.list) %>% t %>% as.data.frame()
  #print(colnames(data))
  #print(gene.list)
  if(length(colnames(data))!=length(gene.list)){
    message(paste0("Gene ",setdiff(gene.list,colnames(data)))," not found in data.")
    gene.list=gene.list[gene.list %in% setdiff(gene.list,colnames(data))]
  }
  data=data[gene.list]
  #data=FetchData(merged.srt,vars=gene.list,cells=rownames(merged.srt@meta.data[merged.srt@meta.data$type==cell.type,]),layer=layer)
  data=rownames_to_column(data,var = "observations") %>% gather(.,key="features",value="values",-observations)
  data$time=strsplit(data$observations,"_") %>% lapply(.,`[`,CT_field) %>% unlist()
  if(!is.null(patient_filed)){
    data$patient=strsplit(data$observations,"_") %>% lapply(.,`[`,patient_filed) %>% lapply(.,paste0,collapse="_") %>% unlist()
  }
  if(filter_data){
    data=data[data$values!=0,]
  }
  data$time_numeric=data$time %>% gsub("CT","",.) %>% as.numeric()
  data$cell_type=cell.type
  return(data)
}

#input are the return of "getWholeDaySplicedData" or "getMultipleWholeDaySplicedData"
#return a single loom.matrix from merging of the list by cbind
#downstream: plotSplicingCounts
mergeSplicedData<-function(loom.matrix.list,normalize=F,num.scale=300000000,normalize_by_splicing_stat=F){
  genes=list()
  res=list()
  if(normalize){
    for(i in 1:length(loom.matrix.list)){
      if(normalize_by_splicing_stat){
        total_spliced=sum(loom.matrix.list[[i]]$spliced)
        total_unspliced=sum(loom.matrix.list[[i]]$unspliced)
        total_ambiguous=sum(loom.matrix.list[[i]]$ambiguous)
        loom.matrix.list[[i]]$spliced=num.scale*loom.matrix.list[[i]]$spliced/total_spliced
        loom.matrix.list[[i]]$unspliced=num.scale*loom.matrix.list[[i]]$unspliced/total_unspliced
        loom.matrix.list[[i]]$ambiguous=num.scale*loom.matrix.list[[i]]$ambiguous/total_ambiguous
      }else{
        total_count=sum(sum(loom.matrix.list[[i]]$spliced),sum(loom.matrix.list[[i]]$unspliced),sum(loom.matrix.list[[i]]$ambiguous))
        loom.matrix.list[[i]]$spliced=num.scale*loom.matrix.list[[i]]$spliced/total_count
        loom.matrix.list[[i]]$unspliced=num.scale*loom.matrix.list[[i]]$unspliced/total_count
        loom.matrix.list[[i]]$ambiguous=num.scale*loom.matrix.list[[i]]$ambiguous/total_count
      }
    }
  }
  for(i in 1:length(loom.matrix.list)){
    genes[[i]]=rownames(loom.matrix.list[[i]]$spliced)
  }
  shared_genes=base::Reduce(intersect,genes)
  for(i in 1:length(loom.matrix.list)){
    loom.matrix.list[[i]]=subsetLoomMatrix(loom.matrix.list[[i]],genes=shared_genes)
  }
  listsplice=list()
  listunsplice=list()
  listambiguous=list()
  for(i in 1:length(loom.matrix.list)){
    listsplice=c(listsplice,loom.matrix.list[[i]]$spliced)
    listunsplice=c(listunsplice,loom.matrix.list[[i]]$unspliced)
    listambiguous=c(listambiguous,loom.matrix.list[[i]]$ambiguous)
  }
  res$spliced=base::Reduce(cbind,listsplice)
  res$unspliced=base::Reduce(cbind,listunsplice)
  res$ambiguous=base::Reduce(cbind,listambiguous)
  #res$spliced=base::Reduce(cbind,list(loom.matrix.list[[1]]$spliced,loom.matrix.list[[2]]$spliced,
  #                                    loom.matrix.list[[3]]$spliced,loom.matrix.list[[4]]$spliced,
  #                                    loom.matrix.list[[5]]$spliced,loom.matrix.list[[6]]$spliced))
  #res$unspliced=base::Reduce(cbind,list(loom.matrix.list[[1]]$unspliced,loom.matrix.list[[2]]$unspliced,
  #                                      loom.matrix.list[[3]]$unspliced,loom.matrix.list[[4]]$unspliced,
  #                                      loom.matrix.list[[5]]$unspliced,loom.matrix.list[[6]]$unspliced))
  #res$ambiguous=base::Reduce(cbind,list(loom.matrix.list[[1]]$ambiguous,loom.matrix.list[[2]]$ambiguous,
  #                                      loom.matrix.list[[3]]$ambiguous,loom.matrix.list[[4]]$ambiguous,
  #                                      loom.matrix.list[[5]]$ambiguous,loom.matrix.list[[6]]$ambiguous))
  return(res)
}

#input are seurat object
#return data.frame, colnames are specified colnames in arguments with "meta.cols"
getMetaData<-function(srt,meta.cols){
  return(srt@meta.data[,meta.cols])
}

detectCircadianVDJ<-function(result.immcantation,out.dir,immune.receptor=NULL){
  locus.testing=NULL
  if(is.null(immune.receptor)){
    locus.testing=c(TCR_LOCUS,BCR_LOCUS)
  }else if(immune.receptor=="bcr"){
    locus.testing=BCR_LOCUS
  }else if(immune.receptor=="tcr"){
    locus.testing=TCR_LOCUS
  }else{
    stop("'immune.receptor' should be NULL/bcr/tcr")
  }
  for(this.locus in locus.testing){
    locus.to.use=getField(this.locus,":",1)
    segment.to.use=getField(this.locus,":",2)
    data.to.use=dplyr::filter(result.immcantation,locus==locus.to.use)
    data.to.use=balanceObsevations(data.to.use,main.group = c("individual"),balance.group=c("CT"))
    data.to.use=data.to.use %>% group_by(individual,CT,across(segment.to.use)) %>% summarise(count=n())
    
    #check if every individual have 6 points
    individual.point.num=(data.to.use[c("individual","CT")] %>% unique())$individual %>% table()
    individual.point.num.valid=individual.point.num[individual.point.num==6] %>% names()
    data.to.use=dplyr::filter(data.to.use,individual %in% individual.point.num.valid)
    
    #change individual ids to RepX
    data.to.use$individual=catSub(data.to.use$individual,unique(data.to.use$individual),paste0("Rep",1:length(unique(data.to.use$individual))))
    print(head(data.to.use))
    data.to.use$reps=paste0(data.to.use$CT,".",data.to.use$individual)
    data.to.use=spread(data.to.use[c(segment.to.use,"reps","count")],key = "reps",value="count")
    
    colnames(data.to.use)[1]="geneSymbol"
    rep.num=(ncol(data.to.use)-1)/6
    data.to.use=data.to.use[c("geneSymbol",lapply(CT_TIME_ORDER,paste0,".",paste0("Rep",1:rep.num)) %>% unlist())]
    data.to.use=data.to.use[!is.na(data.to.use$geneSymbol),]

    print(head(data.to.use))
    data.to.use[is.na(data.to.use)]=0
    write.table(data.to.use, file=paste0(out.dir,locus.to.use,"_",segment.to.use,".txt"),
                sep="\t", quote=FALSE, row.names=FALSE)
    meta2d(infile=paste0(out.dir,locus.to.use,"_",segment.to.use,".txt"), filestyle="txt",
           outdir=out.dir, timepoints=rep(seq(1, 21, by=4), each=rep.num),
           outIntegration="noIntegration")
  }
}

#
detectCircadianGenesByPseudoBulk<-function(int.srt,features=CIRCADIAN_GENES_MAIN,usecelltype,out.dir=NULL,pseudobulk="sampling",rep=1,proportion=1){
  if(grepl("Temra/Teff",usecelltype)){
    usecelltype.char=gsub("Temra/Teff","Temra#Teff",usecelltype)
  }else{
    usecelltype.char=usecelltype
  }
  if(pseudobulk=="sampling"){
    test.data=plotPseudobulkByCT(int.srt,features =features, cell.type = usecelltype,return.data = T,rep = rep,proportion = proportion)
    test.data=test.data[,c("feature","time","downsample","count")]
    test.data$downsample=paste0("Rep",test.data$downsample)
    test.data$colname.new=paste0(test.data$time,".",test.data$downsample)
  }else if(pseudobulk=="individual"){
    test.data=plotPseudobulkByCTByIndividual(int.srt,features =features, cell.type = usecelltype,return.data = T)
    test.data=test.data[,c("feature","time","individual","relative_expression")]
    test.data[is.na(test.data$relative_expression),"relative_expression"]=0
    colnames(test.data)[4]="count"
    test.data.summ=test.data[c("time","individual")] %>% unique() %>% table() %>% as.data.frame()
    reps.absent=test.data.summ[test.data.summ$Freq==0,"individual"] %>% unique() %>% as.vector()
    test.data=test.data[!(test.data$individual %in% reps.absent),]
    test.data$individual=cateSub(test.data$individual,unique(test.data$individual),paste0("Rep",1:length(unique(test.data$individual))))
    test.data$colname.new=paste0(test.data$time,".",test.data$individual)
    rep=length(unique(test.data$individual))
  }else{
    stop("argument 'pseudobulk' should be one of c('sampling','individual')")
  }

  test.data=test.data[c("feature","colname.new","count")]
  test.data=test.data %>% spread(.,key=colname.new,value=count)
  colnames(test.data)[1]="geneSymbol"
  col.order=NULL
  for(i in 1:length(CT_TIME_ORDER)){
    col.order=c(col.order,CT_TIME_ORDER[i] %>% paste0(.,".",paste0("Rep",1:rep)))
  }
  test.data=test.data[c("geneSymbol",col.order)]
  write.table(test.data, file=paste0(out.dir,usecelltype.char,".txt"),
              sep="\t", quote=FALSE, row.names=FALSE)
  meta2d(infile=paste0(out.dir,usecelltype.char,".txt"), filestyle="txt",
         outdir=out.dir, timepoints=rep(seq(1, 21, by=4), each=rep),
         outIntegration="noIntegration")
  
}

#period patterns detection, note that we should assign manual.level2 to type in srt
detectCircadianGenes<-function(int.srt,usecelltype,individuals,features=NULL,HVarGene=T,out.dir="../analysis/OSgenePredict_HEALTH/"){
  rep=length(individuals)
  if(grepl("Temra/Teff",usecelltype)){
    usecelltype.char=gsub("Temra/Teff","Temra#Teff",usecelltype)
  }else{
    usecelltype.char=usecelltype
  }
  if(!file.exists(paste0(out.dir,usecelltype.char,".txt"))){
    metacycle.test=subset(int.srt,type==usecelltype) %>% subset(.,patient %in% individuals)
    metacycle.test=NormalizeData(metacycle.test)
    metacycle.test=ScaleData(metacycle.test)
    if(!is.null(features)){
      wanted_genes=features
    }else{
      if(HVarGene){
        metacycle.test=FindVariableFeatures(metacycle.test)
        wanted_genes=VariableFeatures(metacycle.test)
      }else{
        wanted_genes=rowSums(metacycle.test) %>% sort() %>% tail(.,3000) %>% names()
      }
    }
    metacycle.data=fetchMergedDataOneCellType(metacycle.test,cell.type=usecelltype,gene.list=wanted_genes,CT_field=3,patient_filed=c(1,2))
    metacycle.data.mean=metacycle.data %>% group_by(.,features,time,patient) %>% summarise(.,mean=median(values))
    metacycle.data.mean$rep=case_when(metacycle.data.mean$patient==individuals[1] ~ "Rep1",
                                      metacycle.data.mean$patient==individuals[2] ~ "Rep2",
                                      metacycle.data.mean$patient==individuals[3] ~ "Rep3",
                                      metacycle.data.mean$patient==individuals[4] ~ "Rep4",
                                      metacycle.data.mean$patient==individuals[5] ~ "Rep5",
                                      metacycle.data.mean$patient==individuals[6] ~ "Rep6",)
    metacycle.data.mean$new.colnames=paste(metacycle.data.mean$time,metacycle.data.mean$rep,sep = '.')
    metacycle.data.mean.mtx=metacycle.data.mean[,c("features","mean","new.colnames")] %>% spread(.,key=new.colnames,value=mean)
    colnames(metacycle.data.mean.mtx)[1]="geneSymbol"
    rep.str=paste0("Rep",1:rep)
    order.used=lapply(CT_TIME_ORDER,paste,rep.str,sep='.') %>% unlist()
    if(ncol(metacycle.data.mean.mtx)!=length(order.used)+1){
      return(NULL)
    }
    metacycle.data.mean.mtx=metacycle.data.mean.mtx[c("geneSymbol",order.used)]
    metacycle.data.mean.mtx.rmna=metacycle.data.mean.mtx[!is.na(metacycle.data.mean.mtx[,order.used] %>% rowSums()),]
    if(nrow(metacycle.data.mean.mtx.rmna)==0){
      return(NULL)
    }
    write.table(metacycle.data.mean.mtx.rmna, file=paste0(out.dir,usecelltype.char,".txt"),
                sep="\t", quote=FALSE, row.names=FALSE)
    remove(metacycle.test)
    gc()
  }
  meta2d(infile=paste0(out.dir,usecelltype.char,".txt"), filestyle="txt",
         outdir=out.dir, timepoints=rep(seq(1, 21, by=4), each=rep),
         outIntegration="noIntegration")
}

#filter for elements that only occured once
getOnceOccured<-function(vec){
  duplicated.elements=vec[duplicated(vec)]
  non.duplicated=vec[!(vec %in% duplicated.elements)]
  return(non.duplicated)
}

#inputs are points, points (x1,y1) and (x2,y2) form a segment, and (x3,y3) and (x4,y4)
#form another, this function will determine if the two segment have a point of interection 
havePointIntersect<-function(seg1,seg2){
  x1=seg1[1]
  y1=seg1[2]
  x2=seg1[3]
  y2=seg1[4]
  x3=seg2[1]
  y3=seg2[2]
  x4=seg2[3]
  y4=seg2[4]
  if(((y4-y3)*(x2-x1)-(y2-y1)*(x4-x3))==0){
    return(TRUE)
  }
  y0=((y2-y1)*(y4*x3-y3*x4)-(y4-y3)*(y2*x1-y1*x2))/((y4-y3)*(x2-x1)-(y2-y1)*(x4-x3))
  if(y1==y2&&y3!=y4){
    x0=(y4*x3+(x4-x3)*y0-y3*x4)/(y4-y3)
  }else if(y1!=y2&&y3==y4){
    x0=(y2*x1+(x2-x1)*y0-y1*x2)/(y2-y1)
  }else if(y1!=y2&&y3!=y4){
    x0=(y2*x1+(x2-x1)*y0-y1*x2)/(y2-y1)
  }else{
    return(FALSE)
  }
  if((x0>x1&x0>x2)|(x0>x3&x0>x4)|(y0>y1&y0>y2)|(y0>y3&y0>y4)){
    return(FALSE)
  }else if((x0<x1&x0<x2)|(x0<x3&x0<x4)|(y0<y1&y0<y2)|(y0<y3&y0<y4)){
    return(FALSE)
  }else{
    return(TRUE)
  }
}


isConvex <- function(gene_data) {
  points=gene_data[,c(3,4)]
  n=nrow(points)
  
  # Check for coincident points
  if (nrow(unique(points)) != n) {
    return(FALSE)
  }
  orientation=NULL
  # Check for collinearity and orientation
  for (i in 1:n) {
    v1=points[i,] - points[ifelse(i == 1, n, i - 1),]
    v2=points[ifelse(i == n, 1, i + 1),] - points[i,]
    cross_product=v1[1] * v2[2] - v1[2] * v2[1]
    
    if (cross_product == 0) {
      return(FALSE)  # Points are collinear
    }
    #message(cross_product)
    if (is.null(orientation)) {
      orientation=cross_product > 0
    } else {
      if ((cross_product > 0) != orientation) {
        return(FALSE) # Not all points are in the same orientation
      }
    }
  }
  return(TRUE)
}

isConvexParallel<-function(gene_data_list) {
  cores=detectCores()
  cl=makeCluster(cores)
  registerDoParallel(cl)
  res=foreach(i = 1:length(gene_data_list), .export = "isConvex",.inorder = TRUE) %dopar% {
    isConvex(gene_data_list[[i]])
  }
  stopCluster(cl)
  return(res)
}

isPolygon<-function(gene.data){
  segs=list()
  for(i in 1:(nrow(gene.data)-1)){
    thisseg=c(gene.data$spliced[i],gene.data$unspliced[i],gene.data$spliced[i+1],gene.data$unspliced[i+1])
    segs=c(segs,list(thisseg))
  }
  lastseg=c(gene.data$spliced[1],gene.data$unspliced[1],gene.data$spliced[nrow(gene.data)],gene.data$unspliced[nrow(gene.data)])
  segs=c(segs,list(lastseg))
  polygon=TRUE
  for(j in 3:length(segs)){
    for(k in 1:(j-2)){
      if(havePointIntersect(segs[[k]],segs[[j]])){
        polygon=FALSE
        return(polygon)
      }
    }
    if(!polygon){
      return(polygon)
    }
  }
  return(polygon)
}

isPolygonParallel<-function(gene_data_list) {
  cores=detectCores()
  cl=makeCluster(cores)
  registerDoParallel(cl)
  res=foreach(i = 1:length(gene_data_list), .export = c("isPolygon","havePointIntersect"),.inorder = TRUE) %dopar% {
    isPolygon(gene_data_list[[i]])
  }
  stopCluster(cl)
  return(res)
}

#input: merged spliced data
fetchSingleSplicedCounts<-function(loom.matrix,feature=CIRCADIAN_GENES_MAIN[1],cells=NULL,summarise="mean"){
  if(is.null(cells)){
    cells=colnames(loom.matrix$spliced)
  }
  #spliced=loom.matrix$spliced[feature,cells] %>% t %>% as.data.frame() %>% 
  #  rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  #spliced$individual=getField(spliced$cell,sep="_",c(1,2))
  #spliced$status="spliced"
  #data=spliced
  total=(loom.matrix$spliced[feature,cells]+loom.matrix$unspliced[feature,cells]+
           loom.matrix$ambiguous[feature,cells]) %>% t %>% as.data.frame() %>% 
    rownames_to_column(.,var="feature") %>% gather(.,cell,count,-feature)
  total$individual=getField(total$cell,sep="_",c(1,2))
  total$status="total"
  data=total
  data$time=getField(data$cell,sep="_",3)
  if(summarise=="mean"){
    data=data %>% group_by(.,individual,status,time) %>% summarise(.,mean_count=mean(count))
    ylab="mean_count"
  }else{
    data=data %>% group_by(.,individual,status,time) %>% summarise(.,total_count=sum(count))
    ylab="total_count"
  }
  return(data)
}

getJTK_CYCLEouts<-function(data.dir,celltypes,filter.data=T,period=24,padj=0.01){
  JTK.result=NULL
  for(celltype in celltypes){
    if(grepl("Temra/Teff",celltype)){
      celltype=gsub("Temra/Teff","Temra#Teff",celltype)
    }else{
      celltype=celltype
    }
    if(file.exists(paste0(data.dir,"/",'JTKresult_',celltype,".txt"))){
      JTK.result.part=read.delim(paste0(data.dir,"/",'JTKresult_',celltype,".txt"))
      JTK.result.part$celltype=celltype
    }else{
      message(paste0(paste0(data.dir,"/",'JTKresult_',celltype,".txt")," doesn't exist"))
      next
    }
    if(is.null(JTK.result)){
      JTK.result=JTK.result.part
    }else{
      JTK.result=rbind(JTK.result,JTK.result.part)
    }
  }
  if(filter.data){
    JTK.result=JTK.result%>% dplyr::filter(.,ADJ.P<=padj,PER==period)
  }
  
  return(JTK.result)
}

getAllJTK_CYCLEouts<-function(path,celltypes){
  JTK_CYCL=NULL
  for(block_data in list.files(path)){
    block_path=paste0(path,"/",block_data)
    block=getJTK_CYCLEouts(block_path,celltypes=celltypes,filter.data = F)
    if(is.null(JTK_CYCL)){
      JTK_CYCL=block
    }else{
      JTK_CYCL=rbind(JTK_CYCL,block)
    }
  }
  return(JTK_CYCL)
}

getLSouts<-function(data.dir,celltypes,filter.data=T,period=24,p.cut=0.01){
  LS.result=NULL
  for(celltype in celltypes){
    if(grepl("Temra/Teff",celltype)){
      celltype=gsub("Temra/Teff","Temra#Teff",celltype)
    }else{
      celltype=celltype
    }
    if(file.exists(paste0(data.dir,'LSresult_',celltype,".txt"))){
      LS.result.part=read.delim(paste0(data.dir,'LSresult_',celltype,".txt"))
      LS.result.part$celltype=celltype
    }else{
      next
    }
    if(is.null(LS.result)){
      LS.result=LS.result.part
    }else{
      LS.result=rbind(LS.result,LS.result.part)
    }
  }
  if(filter.data){
    LS.result=LS.result%>% dplyr::filter(.,p<=p.cut,Period>=period-2,Period<=period+2)
  }
  return(LS.result)
}

#add labels to srt by existing labels
#label_mapping is like:
#type (existing col)     new_label
#T cell              1
#B cell              2
#Monocyte            3
addLabels<-function(srt,label_mapping){
  srt@meta.data[,colnames(label_mapping)[2]]<-NA
  for(i in 1:nrow(label_mapping)){
    srt@meta.data[srt@meta.data[colnames(label_mapping)[1]]==label_mapping[i,1],colnames(label_mapping)[2]]=label_mapping[i,2]
  }
  return(srt)
}

#get cells by meta
CellsByMeta<-function(srt,meta.col,meta.content,sample.size){
  cells=rownames(srt@meta.data[srt@meta.data[,meta.col]==meta.content,])
  cells=sample(cells,sample.size)
  return(cells)
}

#transform velosity list to merged srt
velosityList2msrt<-function(list){
  res=list()
  i=1
  message("transforming loom files...")
  for(loom in list){
    message(paste0("loom file: ",i))
    res[[i]]=CreateSeuratObject(loom)
    i=i+1
  }
  message("all loom files transformed.")
  message("merging...")
  res_srt=merge(res[[1]],res[2:length(res)])
  return(res_srt)
}

#balance number of observations(num of rows) in a data.frame
balanceObsevations<-function(df,main.group,balance.group){
  res=NULL
  if(length(main.group)>1){
    main.categories=df[main.group[1]] %>% unique() %>% unlist() %>% as.vector()
    for(subgroup in main.categories){
      sub.df=dplyr::filter(df,get(main.group[1])==subgroup)
      sub.res=balanceObsevations(sub.df,main.group[-1],balance.group)
      if(is.null(res)){
        res=sub.res
      }else{
        res=rbind(res,sub.res)
      }
    }
  }else{
    main.categories=df[main.group] %>% unique() %>% unlist() %>% as.vector()
    for(subgroup in main.categories){
      #message(subgroup)
      sub.df=dplyr::filter(df,get(main.group)==subgroup)
      min.observation=sub.df[balance.group] %>% table() %>% min()
      if(min.observation==0){
        message(paste0("subgroup: ",subgroup," have 0 observation as min count"))
        next
      }
      sub.df=sub.df %>% group_by(.,across(all_of(balance.group))) %>% sample_n(size = min.observation,replace = F)
      #print(head(sub.df))
      if(is.null(res)){
        res=sub.df
      }else{
        res=rbind(res,sub.df)
      }
    }
  }
  return(res)
}

#substitute whole string/int in columns by designated strings/ints with a mapping logic
# col1          col2                 col1        col2
# individual_1   3                   Rep1         3
# individual_1   5                   Rep1         5
# individual2    8          =>       Rep2         8
# individual:3   2                   Rep3         2
# individual:3   9                   Rep3         9
#usage: df$col1<-catSub(df$col1,unique(df$col1),paste0("Rep",1:length(unique(df$col1))))
cateSub<-function(vector,categories.ori,categories.nov){
  if(length(categories.ori)!=length(categories.nov)){
    stop("length of categories is not equal")
  }
  i=1
  for(i in 1:length(categories.ori)){
    vector=gsub(categories.ori[i],categories.nov[i],vector)
  }
  return(vector)
}

srt2bulkMatrix<-function(srt,split.by,layer="count",normalize=F){
  res=NULL
  if(length(split.by)==1){
    reps=srt[[split.by]] %>% unique() %>% unlist() %>% as.vector()
    for(rep in reps){
      this.cells=rownames(srt@meta.data[srt@meta.data[split.by]==rep,])
      df=LayerData(srt,cells=this.cells,layer=layer)%>% rowSums() %>% as.data.frame()
      colnames(df)=rep
      if(normalize){
        df[,rep]=df[,rep]*1000000/sum(df[,rep])
      }
      message(colnames(df))
      df=rownames_to_column(df,var = "genes")
      if(is.null(res)){
        res=df
      }else{
        res=left_join(res,df,by="genes")
      }
      gc()
    }
  }else if(length(split.by)==2){
    reps1=srt[[split.by[1]]] %>% unique() %>% unlist() %>% as.vector()
    reps2=srt[[split.by[2]]] %>% unique() %>% unlist() %>% as.vector()
    for(rep1 in reps1){
      for(rep2 in reps2){
        this.cells=rownames(srt@meta.data[srt@meta.data[split.by[1]]==rep1&srt@meta.data[split.by[2]]==rep2,])
        df=LayerData(srt,cells=this.cells,layer=layer)%>% rowSums() %>% as.data.frame()
        new_colname=paste0(rep1,"-",rep2) %>% gsub(" ","_",.)
        colnames(df)=new_colname
        if(normalize){
          df[,new_colname]=df[,new_colname]*1000000/sum(df[,new_colname])
        }
        message(colnames(df))
        df=rownames_to_column(df,var = "genes")
        if(is.null(res)){
          res=df
        }else{
          res=left_join(res,df,by="genes")
        }
        gc()
      }
    }
  }else if(length(split.by)==3){
    reps1=srt[[split.by[1]]] %>% unique() %>% unlist() %>% as.vector()
    reps2=srt[[split.by[2]]] %>% unique() %>% unlist() %>% as.vector()
    reps3=srt[[split.by[3]]] %>% unique() %>% unlist() %>% as.vector()
    for(rep1 in reps1){
      for(rep2 in reps2){
        for(rep3 in reps3){
          this.cells=rownames(srt@meta.data[srt@meta.data[split.by[1]]==rep1&srt@meta.data[split.by[2]]==rep2&srt@meta.data[split.by[3]]==rep3,])
          df=LayerData(srt,cells=this.cells,layer=layer)%>% rowSums() %>% as.data.frame()
          new_colname=paste0(rep1,"-",rep2,"-",rep3) %>% gsub(" ","_",.)
          colnames(df)=new_colname
          if(normalize){
            df[,new_colname]=df[,new_colname]*1000000/sum(df[,new_colname])
          }
          message(colnames(df))
          df=rownames_to_column(df,var = "genes")
          if(is.null(res)){
            res=df
          }else{
            res=left_join(res,df,by="genes")
          }
          gc()
        }
      }#for(rep2 in reps2)
    }
  }else{
    stop("split.by exceeds its max length.")
  }
  rownames(res)=res$genes
  res$genes=NULL
  return(res)
}

#gnerate Adjacency Matrix for net work buiding
generateAdjacencyMatrix<-function(JTK_CYCL_out,spread=F){
  celltypes=JTK_CYCL_out$celltype %>% unique()
  share.data=NULL
  for(type.this in celltypes){
    for(type.that in celltypes){
      n.gene.type.this=nrow(JTK_CYCL_out[JTK_CYCL_out$celltype==type.this,])
      n.share.gene=length(intersect(JTK_CYCL_out[JTK_CYCL_out$celltype==type.this,"CycID"],
                                    JTK_CYCL_out[JTK_CYCL_out$celltype==type.that,"CycID"]))
      if(n.gene.type.this==0){
        this.share.data=data.frame(from.type=type.this,to.type=type.that,shared.ratio=NA)
      }else{
        this.share.data=data.frame(from.type=type.this,to.type=type.that,shared.ratio=n.share.gene/n.gene.type.this)
      }
      if(is.null(share.data)){
        share.data=this.share.data
      }else{
        share.data=rbind(share.data,this.share.data)
      }
    }
  }
  if(spread){
    share.data=spread(share.data,key=to.type,value=shared.ratio)
    rownames(share.data)=share.data$from.type
    share.data$from.type=NULL
  }
  print(head(share.data[1:3,1:3]))
  return(share.data)
}

#plotNetWorkCircadian
plotNetWorkCircadian<-function(srt,JTK.result.filtered){
  nFeatureType=srt@meta.data[,c("nFeature_RNA","type")] %>% group_by(type) %>% summarise(median_genes=max(nFeature_RNA))
  nCircadianType=JTK.result.filtered$celltype %>% table %>% as.data.frame() %>% `colnames<-`(c("type","circadian_genes"))
  meta.data=left_join(nCircadianType,nFeatureType) %>% mutate(rhythmic_gene_ratio=circadian_genes/median_genes)
  colnames(meta.data)[1]="from.type"
  
  AM=generateAdjacencyMatrix(JTK.result.filtered)
  AM=AM %>% dplyr::filter(shared.ratio>0.05)
  AM=left_join(AM,meta.data) %>% as.data.frame()
  
  graph=graph_from_data_frame(AM, directed = TRUE)
  V(graph)$rhythmic_gene_ratio=AM$rhythmic_gene_ratio[match(V(graph)$name, AM$from.type)]
  # Plot with ggraph
  ggraph(graph, layout = "circle") +  # Using "circle" layout for clarity
    # Use geom_edge_fan to make edges dodge each other
    geom_edge_fan(aes(width = shared.ratio,alpha=shared.ratio),
                  arrow = arrow(length = unit(2, "mm")),  # Arrows for direction
                  end_cap = circle(3, "mm"),
                  start_cap =circle(2, "mm")) +           # Avoid overlap with nodes
    # Scale edge width
    scale_edge_width(range = c(0.5, 2)) +
    # Add nodes and labels
    geom_node_point(aes(size=rhythmic_gene_ratio))+
    geom_node_text(aes(label = name, x = x*1.25, y = y*1.15), 
                   size = 4,  # Adjust label size if needed
                   repel = FALSE)+
    # Customize theme
    theme_graph()
}

#plotRadarCellType
plotRadarCellType<-function(srt,celltypes,feature,colors=generateColor(25)){
  plotdata=NULL
  for(celltype in celltypes){
    this.data=plotPseudobulkByCTByIndividual(srt=srt,cell.type=celltype,features=feature,return.data = T)
    this.data=this.data[c("time","CT_median")] %>% unique()
    this.data$cell.type=celltype
    if(is.null(plotdata)){
      plotdata=this.data
    }else{
      plotdata=rbind(plotdata,this.data)
    }
  }
  plotdata=spread(plotdata,key=cell.type,value=CT_median) %>% as.data.frame()
  rownames(plotdata)=plotdata$time
  plotdata$time=NULL
  plotdata=10^plotdata
  plotdata=t(plotdata) %>% as.data.frame()
  plotdata$cell_type=rownames(plotdata)
  plotdata=plotdata[,c("cell_type",CT_TIME_ORDER)]
  print(head(plotdata))
  plotlist=list()
  for(i in 1:length(celltypes)){
    max.value=max(plotdata[i,CT_TIME_ORDER])
    median.value=1
    min.value=min(plotdata[i,CT_TIME_ORDER])
    plotlist[[i]]=ggradar(plotdata[i,],grid.min = min.value,grid.mid = 1,grid.max = max.value,fill = T,fill.alpha = 0.5,label.gridline.min = F,label.gridline.mid = F,
            label.gridline.max = F,group.colours = colors[i],group.point.size = 2,group.line.width = 1,legend.text.size=12,axis.label.size=3,
            grid.label.size = 5)+ggtitle(celltypes[i])+theme(plot.title = element_text(size=11))
  }
  ggarrange(plotlist=plotlist)
  #plotdata$time<-factor(plotdata$time,levels = CT_TIME_ORDER)
  #ggplot(plotdata,aes(x=time,y=CT_median,group=cell.type,color=cell.type))+coord_equal()+coord_polar()+
  #  geom_point()+geom_line()+
  #  theme_classic()
}

# Function: plotCountOscilattingByMetaCelltype
# plot bar plot showing count of circadian genes for each cell type
##
# upstream: readJTKFromMetaCells
# downstream: NSF
# dependency: NSF
# caller: NSF
plotCountOscilattingByMetaCelltype<-function(JTK.result,threshold.ADJ.P=0.05,PER.floor=20,PER.ceiling=28){
  AllJTKresult.filtered=filter(JTK.result,ADJ.P<threshold.ADJ.P,PER>=PER.floor,PER<=PER.ceiling)
  #filter for genes oscillate in at least 2 people
  filtered.genes=AllJTKresult.filtered[c("CycID","celltype")] %>% table() %>% as.data.frame() %>% dplyr::filter(.,Freq>=2)
  filtered.genes=filtered.genes$CycID %>% unique()
  AllJTKresult.filtered=AllJTKresult.filtered[AllJTKresult.filtered$CycID %in% filtered.genes,]
  plotdata=(AllJTKresult.filtered[c("CycID","celltype")] %>% unique())$celltype %>% table() %>% as.data.frame()
  colnames(plotdata)[1]="cell_type"
  plotdata$main_type=celltypes[as.vector(plotdata$cell_type)] %>% as.character()
  print(head(plotdata))
  ggplot(plotdata)+geom_bar(aes(x=cell_type,y=Freq,fill=main_type),stat="identity",color="black")+
    scale_y_continuous(expand = c(0,0),limits = c(0,max(plotdata$Freq)+100))+scale_x_discrete(limits=names(celltypes))+
    theme(axis.text.x = element_text(angle=60,hjust=1),panel.background = element_rect(fill="white",color="black"),
          panel.grid = element_line(color="grey"))+
    scale_fill_manual(values = generateColor(n = length(unique(plotdata$main_type))))+
    ylab("# of oscillating genes")+ xlab("")+NoLegend()
}

# Function: plotCountOscilattingByMetaCelltypeByCT
# plot bar plot showing count of circadian genes at each CT that reaching peaks for each cell type
##
# upstream: readJTKFromMetaCells
# downstream: NSF
# dependency: NSF
# caller: NSF
plotCountOscilattingByMetaCelltypeByCT<-function(JTK.result,threshold.ADJ.P=0.05,PER.floor=20,PER.ceiling=28){
  AllJTKresult.filtered=filter(JTK.result,ADJ.P<threshold.ADJ.P,PER>=PER.floor,PER<=PER.ceiling)
  #filter for genes oscillate in at least 2 people
  filtered.genes=AllJTKresult.filtered[c("CycID","celltype")] %>% table() %>% as.data.frame() %>% dplyr::filter(.,Freq>=2)
  filtered.genes=filtered.genes$CycID %>% unique()
  AllJTKresult.filtered=AllJTKresult.filtered[AllJTKresult.filtered$CycID %in% filtered.genes,]
  plotdata=AllJTKresult.filtered[c("celltype","LAG","individual")] %>% table() %>% as.data.frame()
  print(head(plotdata))
  ggplot(plotdata)+geom_bar(aes(x=LAG,y=Freq,fill=celltype),stat="identity")+
    facet_wrap(~individual+celltype,scale="free_y")
}

# Function: plotEnrichGObyCT
#
##
# upstream: readJTKFromMetaCells
# downstream: NSF
# dependency: enrichGObyHGNC
# caller: NSF
plotEnrichGObyCT<-function(JTK.result,type,threshold.ADJ.P=0.05,PER.floor=20,PER.ceiling=28,return.data=F){
  JTKresult.filtered=filter(JTK.result,ADJ.P<threshold.ADJ.P,PER>=PER.floor,PER<=PER.ceiling)
  #filter for genes oscillate in at least 2 people
  filtered.genes=JTKresult.filtered[c("CycID","celltype")] %>% table() %>% as.data.frame() %>% dplyr::filter(.,Freq>=2)
  filtered.genes=filtered.genes$CycID %>% unique()
  JTKresult.filtered=JTKresult.filtered[JTKresult.filtered$CycID %in% filtered.genes,]
  print(head(JTKresult.filtered))
  plotdata=NULL
  for(lag in seq(0,24,2)){
    genes=dplyr::filter(JTKresult.filtered,LAG==lag,celltype==type)$CycID %>% unique
    print(head(genes))
    if(length(genes)==0){
      next
    }
    out=enrichGObyHGNC(genes)
    out=dplyr::filter(out@result,p.adjust<0.01)
    print(head(out))
    if(nrow(out)==0){
      next
    }
    out$Lag=lag
    if(is.null(plotdata)){
      plotdata=out
    }else{
      plotdata=rbind(plotdata,out)
    }
  }
  plotdata$category=Ontology(plotdata$ID)
  #plotdata=plotdata %>% group_by(Lag) %>% arrange(p.adjust) %>% top_n(-10,p.adjust)
  
  #plotdata=plotdata[,c("Lag","ID")]
  #print(head(plotdata))
  #n_total_ID=length(plotdata$ID)
  #n_blank=floor(0.3*n_total_ID)
  #plotdata=rbind(plotdata,data.frame("Lag"=sort(rep(seq(0,22,2),n_blank)),"ID"=rep(paste0("AO:",1:n_blank),12)))
  #plotdata$Lag=as.factor(plotdata$Lag)
  #plotdata$plotpadding=ifelse(grepl("AO",plotdata$ID),"blank","data")
  if(return.data){
    return(plotdata)
  }
  ggplot(plotdata)+geom_tile(aes(x=Lag,y=fct_infreq(ID),fill=plotpadding,color=plotpadding))+NoLegend()+
      geom_text(data=plotdata[plotdata$plotpadding!="blank",],aes(x=Lag,y=fct_infreq(ID),label=fct_infreq(ID)))+
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.background = element_rect(fill="white"))+
      scale_fill_manual(values = c("white","lightblue"))+scale_color_manual(values = c("white","black"))+ylab("")+xlab("")
}

# Function: enrichGObyHGNC
# change HGNC gene symbols to ENTREZ ID and then query for GO terms
##
# dependency: NSF
# caller: plotEnrichGObyCT
# upstream: NSF
# downstream: NSF
enrichGObyHGNC<-function(genes,backgroud.genes=NULL,ont="BP"){
  gene_ids=bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID
  if(!is.null(backgroud.genes)){
    gene_ids_bg=bitr(backgroud.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID
    out=enrichGO(gene_ids,OrgDb = org.Hs.eg.db,universe = gene_ids_bg,ont = ont,pool = T)
  }else{
    out=enrichGO(gene_ids,OrgDb = org.Hs.eg.db,ont = ont,pool = T)
  }
  return(out)
}

#FetchExpressedGenesCellType
#fetch a list of cell type expressed genes
FetchExpressedGenesCellType<-function(srt,cell.type,pct.express=0.05,idents="type"){
  cells=rownames(srt@meta.data)
  labels=srt@meta.data[,idents]
  names(labels)=cells
  Idents(srt)=labels
  subsrt=subset(srt,idents = cell.type)
  mat=LayerData(subsrt,layer = "counts")
  gene_expression_percentage=Matrix::rowSums(mat > 1) / ncol(mat)
  expressed_genes=names(gene_expression_percentage[gene_expression_percentage >= pct.express])
  remove(subsrt)
  gc()
  return(expressed_genes)
}



# Function: modAzimuthAnnotation
# Create predicted.celltype.l1.5 on the basis of predicted.celltype.l2, only work for data annotated with pbmcref
##
# upstream: <Azimuth>RunAzimuth
# downstream: NSF
# dependency: NSF
# caller: NSF
modAzimuthAnnotation<-function(srt){
  srt$predicted.celltype.l1.5=case_when(srt$predicted.celltype.l2=="ASDC" ~ "ASDC",
                                        srt$predicted.celltype.l2=="B intermediate" ~ "B memory",
                                        srt$predicted.celltype.l2=="B memory" ~ "B memory",
                                        srt$predicted.celltype.l2=="B naive" ~ "B naive",
                                        srt$predicted.celltype.l2=="CD14 Mono" ~ "CD14 Mono",
                                        srt$predicted.celltype.l2=="CD16 Mono" ~ "CD16 Mono",
                                        srt$predicted.celltype.l2 %in% c("CD4 CTL","CD4 Naive","CD4 Proliferating","CD4 TCM","CD4 TEM","Treg") ~ "CD4 T",
                                        srt$predicted.celltype.l2 %in% c("CD8 Naive","CD8 Proliferating","CD8 TCM","CD8 TEM") ~ "CD8 T",
                                        srt$predicted.celltype.l2=="cDC1" ~ "cDC1",
                                        srt$predicted.celltype.l2=="cDC2" ~ "cDC2",
                                        srt$predicted.celltype.l2=="dnT" ~ "dnT",
                                        srt$predicted.celltype.l2=="Eryth" ~ "Eryth",
                                        srt$predicted.celltype.l2=="gdT" ~ "gdT",
                                        srt$predicted.celltype.l2=="HSPC" ~ "HSPC",
                                        srt$predicted.celltype.l2=="ILC" ~ "ILC",
                                        srt$predicted.celltype.l2=="MAIT" ~ "MAIT",
                                        srt$predicted.celltype.l2 %in% c("NK","NK Proliferating","NK_CD56bright") ~ "NK",
                                        srt$predicted.celltype.l2=="pDC" ~ "pDC",
                                        srt$predicted.celltype.l2=="Plasmablast" ~ "Plasma",
                                        srt$predicted.celltype.l2=="Platelet" ~ "Platelet"
  )
  return(srt)
}

