library("Seurat")
library("MetaCycle")
library("dplyr")
library("tidyverse")
library("groupdata2")

time_point<-c(9,13,17,21,1,5)
time_point_inorder<-c(1,5,9,13,17,21)
CT_TIME<-paste0("CT",time_point)
CT_TIME_ORDER<-paste0("CT",time_point_inorder)

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
    test.data=plotPseudobulkByCTByIndividual(int.srt,features =features, cell.type = usecelltype,return.data = T,proportion=1)
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

plotPseudobulkByCTByIndividual<-function(srt,cell.type,features,individuals=NULL,normalize.data=F,
                                         return.data=F,proportion=0.8,split.by="no",layer="count",relative=T){
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
      plotdata=plotdata %>% group_by(feature,individual) %>% mutate(relative_expression=(normalized_count-mean(normalized_count))/sd(normalized_count))
      #plotdata$relative_expression=(plotdata %>% group_by(feature,individual) %>% reframe(relative_expression=scale(normalized_count)))$relative_expression[,1]
    }else{
      plotdata=plotdata %>% group_by(feature,individual) %>% mutate(relative_expression=(count-mean(count))/sd(count))
      #plotdata$relative_expression=(plotdata %>% group_by(feature,individual) %>% reframe(relative_expression=scale(count)))$relative_expression[,1]
    }
  }else{
    plotdata$relative_expression=plotdata$normalized_count
  }

  plotdata$time_numeric=plotdata$time %>% gsub("CT","",.) %>% as.numeric()
  if(return.data){
    return(plotdata)
  }
  print(head(plotdata))
  if(split.by=="individual"){
    ggplot(plotdata)+geom_line(aes(x=time,y=relative_expression,group=feature))+
      geom_point(aes(x=time,y=relative_expression,size=count,color=feature))+
      scale_x_discrete(limits=CT_TIME_ORDER)+scale_color_manual(values=generateColor(length(features)))+
      facet_wrap(~individual,scales="free_y")+guides(size=guide_legend(title="UMI count"))+ggtitle(paste0(cell.type))
  }else{
    ggplot(plotdata)+geom_line(aes(x=time,y=relative_expression,group=individual))+
      geom_point(aes(x=time,y=relative_expression,size=count,color=feature),alpha=0.5)+
      scale_x_discrete(limits=CT_TIME_ORDER)+scale_color_manual(values=generateColor(length(features)))+
      facet_wrap(~feature,scales="free_y")+guides(size=guide_legend(title="UMI count"))+ggtitle(paste0(cell.type))
  }
}

fetchMergedDataOneCellType<-function(merged.srt,gene.list=CIRCADIAN_GENES_MAIN,cell.type=NULL,CT_field=1,patient_filed=NULL,layer="data",filter_data=T){
  data=LayerData(merged.srt,layer=layer,cells=rownames(merged.srt@meta.data[merged.srt@meta.data$type==cell.type,]),features=gene.list) %>% as.data.frame() %>% t %>% as.data.frame()
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

plotPseudobulkByCT<-function(srt,cell.type,features=CIRCADIAN_GENES_MAIN,return.data=F,rep=3,proportion=0.25){
  data=fetchMergedDataOneCellType(srt,cell.type=cell.type,gene.list=features,CT_field=3,patient_filed = c(1,2),layer = "count",filter_data=F)
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
      block_data=downsample(block_data,"time") %>% group_by(.,time) %>% summarise(.,count=sum(values))
      #block_data=downsample(block_data,"time") %>% group_by(.,time)
      #block_data=block_data[block_data$values>0,] %>% summarise(.,count=n())
      block_data=block_data[block_data$time!="CTpseudo",]
      #print(block_data)
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
  normalize_data=srt@meta.data[c("nCount_RNA","CT")] %>% group_by(.,CT) %>% summarise(sum=sum(nCount_RNA))
  normalize_data$sum=normalize_data$sum/min(normalize_data$sum)
  colnames(normalize_data)=c("time","sum")
  plotdata=left_join(plotdata,normalize_data,by="time")
  plotdata$normalized_count=plotdata$count/plotdata$sum
  plotdata$normalized_count=plotdata$normalized_count/(1+min(plotdata$normalized_count))
  print(head(plotdata))
  ggplot(plotdata)+geom_boxplot(aes(x=time,y=normalized_count))+
    geom_point(aes(x=time,y=normalized_count,color=downsample))+
    facet_wrap(~feature,scales="free_y")+
    scale_x_discrete(limits=CT_TIME_ORDER)+ggtitle(cell.type)+scale_color_aaas()+
    theme(axis.text.x=element_text(angle=45,hjust=1),
          panel.background=element_rect(fill="white",colour = "black"),
          panel.grid = element_line(color="grey"))+
    NoLegend()+xlab("")
}


args<-commandArgs(trailingOnly = TRUE)
gene_block_index<-args[1] %>% as.numeric()
#gene_block_index<-4

setwd("/lustre/home/acct-medll/medll/data/analysis/OSgenePredict_PseudoBulk_20k_all_genes_NI_by_individual/")
command_exe<-paste0("mkdir block",gene_block_index)
system(command_exe)
outdir<-paste0("/lustre/home/acct-medll/medll/data/analysis/OSgenePredict_PseudoBulk_20k_all_genes_NI_by_individual/block",gene_block_index,"/")

test<-readRDS("/lustre/home/acct-medll/medll/data/analysis/all.samples.oddsremoved.tsneran.Azimuth.healthonly.rds")
#DefaultAssay(test)<-"unspliced"

celltypes<-test$manual_NI %>% unique() %>% sort()
test$type<-test$manual_NI

gene_ids<-((gene_block_index-1)*1000+1):(gene_block_index*1000)
total_gene_list<-Features(test)

genes_wanted<-NULL
if(gene_ids[length(gene_ids)]>length(total_gene_list)){
  genes_wanted=total_gene_list[((gene_block_index-1)*1000+1):length(total_gene_list)]
}else{
  genes_wanted=total_gene_list[gene_ids]
}


for(i in 1:length(celltypes)){
  usecelltype=celltypes[i]
  if(file.exists(paste0(outdir,usecelltype,".txt"))|usecelltype=="Eryth"){
    next
  }
  message(usecelltype)
  detectCircadianGenesByPseudoBulk(test,features=genes_wanted,
      usecelltype = usecelltype,out.dir = outdir,pseudobulk="individual")
}
