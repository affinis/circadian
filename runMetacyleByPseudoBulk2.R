library("Seurat")
library("MetaCycle")
library("dplyr")
library("tidyverse")
library("groupdata2")

time_point<-c(9,13,17,21,1,5)
time_point_inorder<-c(1,5,9,13,17,21)
CT_TIME<-paste0("CT",time_point)
CT_TIME_ORDER<-paste0("CT",time_point_inorder)

detectCircadianGenesByPseudoBulk<-function(int.srt,features=CIRCADIAN_GENES_MAIN,usecelltype,out.dir=NULL,rep=1,proportion=1){
  if(grepl("Temra/Teff",usecelltype)){
    usecelltype.char=gsub("Temra/Teff","Temra#Teff",usecelltype)
  }else{
    usecelltype.char=usecelltype
  }
  test.data=plotPseudobulkByCT(int.srt,features =features, cell.type = usecelltype,return.data = T,rep = rep,proportion = proportion)
  test.data=test.data[,c("feature","time","downsample","count")]
  #return(test.data)
  test.data$downsample=paste0("Rep",test.data$downsample)
  #for(i in 1:length(PATIENTS)){
  #  test.data[test.data$downsample==as.character(i),"downsample"]=paste0("Rep",i)
  #}
  test.data$colname.new=paste0(test.data$time,".",test.data$downsample)
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
#gene_block_index<-1

setwd("/lustre/home/acct-medll/medll/data/analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5_spliced/")
command_exe<-paste0("mkdir block",gene_block_index)
system(command_exe)
outdir<-paste0("/lustre/home/acct-medll/medll/data/analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5_spliced/block",gene_block_index,"/")

test<-readRDS("/lustre/home/acct-medll/medll/data/analysis/all.samples.oddsremoved.tsneran.Azimuth.velosity.rds")
DefaultAssay(test)<-"spliced"

celltypes<-test$predicted.celltype.l1.5 %>% unique() %>% sort()
test$type<-test$predicted.celltype.l1.5

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
      usecelltype = usecelltype,out.dir = outdir,
      rep=3,proportion=0.3)
}
