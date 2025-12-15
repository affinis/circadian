source('~/script/circadian/circadian_core.R')
#vdj analysis
#read BCR data
BCR_data<-NULL
for(individual in INDIVIDUALS_BATCH1){
  for(i in CT_TIME_ORDER){
    file=file.path(paste0(individual,"_",i),"outs/per_sample_outs",paste0(individual,"_",i),"vdj_b/vdj_results/filtered_contig_igblast_db-pass.tsv")
    if(file.exists(file)){
      BCR_sample_data=readChangeoDb(file)
      BCR_sample_data$individual=individual
      BCR_sample_data$CT=i
      BCR_sample_data$cell_id=paste0(individual,"_",i,"_",BCR_sample_data$cell_id)
      if(is.null(BCR_data)){
        BCR_data=BCR_sample_data
      }else{
        BCR_data=rbind(BCR_data,BCR_sample_data)
      }
    }else{
      message(paste0(file," not found"))
      next
    }
  }
}

test2<-readRDS('~/analysis/circadian/R/4multiplexedindividual.srt.annotated.rds')
BCR_data2=NULL
for(sample in names(BATCH_KD)){
    file=file.path(sample,"outs/per_sample_outs",sample,"vdj_b/vdj_results/filtered_contig_igblast_db-pass.tsv")
    if(file.exists(file)){
      BCR_sample_data=readChangeoDb(file)
      BCR_sample_data$sample=sample
      BCR_sample_meta=test2@meta.data[test2@meta.data$sample==sample,c("individual","CT")]
      BCR_sample_meta$cell_id=getField(rownames(BCR_sample_meta),"_",4)
      BCR_sample_data=left_join(BCR_sample_data,BCR_sample_meta,by="cell_id")
      BCR_sample_data$cell_id=paste0(BCR_sample_data$individual,"_",BCR_sample_data$CT,"_",BCR_sample_data$cell_id)
      BCR_sample_data=dplyr::filter(BCR_sample_data,!is.na(individual))
      if(is.null(BCR_data2)){
        BCR_data2=BCR_sample_data
      }else{
        BCR_data2=rbind(BCR_data2,BCR_sample_data)
      }
    }else{
      message(paste0(file," not found"))
      next
    }
}
BCR_data2$sample<-NULL
BCR_data<-rbind(BCR_data,BCR_data2)

BCR_prim<-BCR_data %>% dplyr::filter(productive)
# remove cells with multiple heavy chain
multi_heavy<-table(dplyr::filter(BCR_data, locus == "IGH")$cell_id)
multi_heavy_cells<-names(multi_heavy)[multi_heavy > 1]
BCR_prim <- dplyr::filter(BCR_data,!cell_id %in% multi_heavy_cells)

# split cells by heavy and light chains
heavy_cells<-dplyr::filter(BCR_prim, locus == "IGH")$cell_id
light_cells<-dplyr::filter(BCR_prim, locus == "IGK" | locus == "IGL")$cell_id
no_heavy_cells<-light_cells[which(!light_cells %in% heavy_cells)]

BCR_prim<-dplyr::filter(BCR_prim, !cell_id %in% no_heavy_cells)

## about 5 min
dist_nearest<-distToNearest(dplyr::filter(BCR_prim,locus=="IGH"),nproc=1,cellIdColumn = "cell_id")

# find threshold for cloning automatically
## about 5 min
threshold_output <- findThreshold(dist_nearest$dist_nearest,
                                  method = "gmm", model = "gamma-norm",
                                  cutoff = "user", spc = 0.995,progress = T)
#threshold=0.1218473
threshold<-threshold_output@threshold

# call clones using hierarchicalClones, need package 'scoper'
## about 5min
results_raw <- hierarchicalClones(BCR_prim, cell_id = 'cell_id',
                              threshold = threshold, only_heavy = FALSE,
                              split_light = TRUE, summarize_clones = FALSE)
results<-results_raw
results$sample<-getField(results$cell_id,"_",1:2)
results$sequence_id<-paste(results$sequence_id,results$sample,sep="_")
celltype.meta<-test$type %>% as.data.frame()
colnames(celltype.meta)[1]<-"cell_type"
celltype.meta<-rownames_to_column(celltype.meta,var="cell_id")
celltype.meta$cell_id<-gsub("TF_","",celltype.meta$cell_id)
results<-left_join(results,celltype.meta)
results<-dplyr::filter(results,!is.na(cell_type))
results<-dplyr::filter(results,cell_type %in% c("B memory","B naive","B intermediate",'Plasmablast'))
saveRDS(results,"~/analysis/circadian/R/BCR_filtered.clone_type.16individual.rds")

#read TCR data
TCR_data<-NULL
for(individual in INDIVIDUALS_BATCH1){
  for(i in CT_TIME_ORDER){
    file=file.path(paste0(individual,"_",i),"outs/per_sample_outs",paste0(individual,"_",i),"vdj_t/vdj_results/filtered_contig_igblast_db-pass.tsv")
    if(file.exists(file)){
      TCR_sample_data=readChangeoDb(file)
      TCR_sample_data$individual=individual
      TCR_sample_data$CT=i
      TCR_sample_data$cell_id=paste0(individual,"_",i,"_",TCR_sample_data$cell_id)
      if(is.null(TCR_data)){
        TCR_data=TCR_sample_data
      }else{
        TCR_data=rbind(TCR_data,TCR_sample_data)
      }
    }else{
      message(paste0(file," not found"))
      next
    }
  }
}

test2<-readRDS('~/analysis/circadian/R/4multiplexedindividual.srt.annotated.rds')
TCR_data2=NULL
for(sample in names(BATCH_KD)){
  file=file.path(sample,"outs/per_sample_outs",sample,"vdj_t/vdj_results/filtered_contig_igblast_db-pass.tsv")
  if(file.exists(file)){
    TCR_sample_data=readChangeoDb(file)
    TCR_sample_data$sample=sample
    TCR_sample_meta=test2@meta.data[test2@meta.data$sample==sample,c("individual","CT")]
    TCR_sample_meta$cell_id=getField(rownames(TCR_sample_meta),"_",4)
    TCR_sample_data=left_join(TCR_sample_data,TCR_sample_meta,by="cell_id")
    TCR_sample_data$cell_id=paste0(TCR_sample_data$individual,"_",TCR_sample_data$CT,"_",TCR_sample_data$cell_id)
    TCR_sample_data=dplyr::filter(TCR_sample_data,!is.na(individual))
    if(is.null(TCR_data2)){
      TCR_data2=TCR_sample_data
    }else{
      TCR_data2=rbind(TCR_data2,TCR_sample_data)
    }
  }else{
    message(paste0(file," not found"))
    next
  }
}
TCR_data2$sample<-NULL
TCR_data<-rbind(TCR_data,TCR_data2)

TCR_data<-TCR_data %>% dplyr::filter(productive)
TCR_prim<-TCR_data
rm(TCR_data)
rm(TCR_data2)
gc()
# remove cells with multiple alpha/beta chain
multi_alpha<-table(dplyr::filter(TCR_prim, locus == "TRA")$cell_id)
multi_beta<-table(dplyr::filter(TCR_prim, locus == "TRB")$cell_id)
multi_chain_cells<-c(names(multi_alpha)[multi_alpha > 1],names(multi_beta)[multi_beta > 1]) %>% unique()
TCR_prim <- dplyr::filter(TCR_prim,!cell_id %in% multi_chain_cells)

# split cells by heavy and light chains
alpha_cells<-dplyr::filter(TCR_prim, locus == "TRA")$cell_id
beta_cells<-dplyr::filter(TCR_prim, locus == "TRB")$cell_id
no_pair_cells<-alpha_cells[which(!alpha_cells %in% beta_cells)]

TCR_prim<-dplyr::filter(TCR_prim, !(cell_id %in% no_pair_cells))
# used 6 hours
dist_nearest<-distToNearest(TCR_prim,nproc=1,locusValues = c("TRB"),cellIdColumn = "cell_id")
saveRDS(dist_nearest,'~/analysis/circadian/R/VDJ.T.dist_nearest.tmp.rds')
dist_nearest<-dist_nearest[!is.na(dist_nearest$dist_nearest),]
# find threshold for cloning automatically
threshold_output <- findThreshold(dist_nearest$dist_nearest,
                                  method = "gmm", model = "gamma-norm",
                                  cutoff = "user", spc = 0.995,progress = T)
#threshold=0.035
threshold<-0.035
#threshold<-threshold_output@threshold

# call clones using hierarchicalClones
# about 8 hours
results <- hierarchicalClones(TCR_prim, cell_id = 'cell_id',
                              threshold = threshold, only_heavy = FALSE,
                              split_light = TRUE, summarize_clones = FALSE)
results$sample<-getField(results$cell_id,"_",1:2)
results$sequence_id<-paste(results$sequence_id,results$sample,sep="_")
test<-readRDS('~/analysis/circadian/R/16individual.srt.annotated.rds')
celltype.meta<-test$predicted.celltype.l2 %>% as.data.frame()
colnames(celltype.meta)[1]<-"cell_type"
celltype.meta<-rownames_to_column(celltype.meta,var="cell_id")
celltype.meta$cell_id<-gsub("TF_","",celltype.meta$cell_id)
results.bkp<-results
results<-results.bkp
results<-left_join(results,celltype.meta)
results.filtered<-dplyr::filter(results,cell_type %in% c("CD4 CTL","CD4 Naive","CD4 Proliferating",
                                                         "CD4 TCM","CD4 TEM","CD8 Naive","CD8 Proliferating",
                                                         "CD8 TCM","CD8 TEM","dnT","gdT","ILC","MAIT","Treg"))
saveRDS(results.filtered,"~/analysis/circadian/R/TCR_filtered.clone_type.16individual.rds")
saveRDS(results,"~/analysis/circadian/R/TCR_RAW.clone_type.16individual.rds")

# merge initial 12 individual data
test.srt<-mergeRawSamplesOneDay("HZD",data.path = '~/analysis/cellranger/',save.each.sample = F,over.write = F,write.path = '~/analysis/circadian/R/',addmodulescore = F,no.processing = T)
rm(test.srt)
gc()
srts<-list()
individual_i=1
for (individual in names(AGES)) {
  srts[[individual_i]]=mergeRawSamplesOneDay(individual,data.path = '~/analysis/cellranger/',save.each.sample = F,over.write = F,write.path = '~/analysis/circadian/R/',addmodulescore = F,no.processing = T)
  individual_i=individual_i+1
}
srts=merge(srts[[1]],srts[2:length(srts)])
saveRDS(srts,'~/analysis/circadian/R/12individual.srt.raw.rds')

# add annotations to raw srt
srt<-readRDS('~/analysis/circadian/R/12individual.raw.rds')
anno<-read.delim('~/analysis/circadian/R/cell.annotations.manual1_2_NI.predicted2.modforSeacells.tsv')
rownames(anno)<-anno$cell_id
srt<-AddMetaData(srt,anno)
srt
srt<-subset(srt,!is.na(manual.level2))
saveRDS(srt,'~/analysis/circadian/R/12individual.srt.annotated.rds')


# merge 4 individual in batch3
i=1
srts<-list()
for (file in list.files('~/analysis/circadian/R/',pattern = "X5_")){
  this.srt=readRDS(file.path('~/analysis/circadian/R',file))
  srts[[i]]=subset(this.srt,individual!="unknown")
  i=i+1
}
srts=merge(srts[[1]],srts[2:length(srts)])
saveRDS(srts,'~/analysis/circadian/R/4multiplexedindividual.srt.raw.rds')
anno<-read.delim('~/analysis/circadian/R/cell.annotation.clean.sct.Azimuth.tsv')
rownames(anno)<-anno$cell_id
srts<-AddMetaData(srts,anno)
srts<-subset(srts,!is.na(predicted.celltype.l2))
srts$sample<-BATCH[paste0(srts$individual,"-",srts$CT)] %>% as.vector()
saveRDS(srts,'~/analysis/circadian/R/4multiplexedindividual.srt.annotated.rds')

srts<-readRDS('~/analysis/circadian/R/4multiplexedindividual.srt.annotated.rds')
srts2<-readRDS('~/analysis/circadian/R/12individual.srt.annotated.rds')
srts<-merge(srts,srts2)
srts[["RNA"]]<-JoinLayers(srts[["RNA"]])
saveRDS(srts,'~/analysis/circadian/R/16individual.srt.annotated.rds')

srts<-readRDS('~/analysis/circadian/R/16individual.srt.annotated.rds')
srts<-integrateSubset(srts)
# redo gene enrichment analysis
# regardless of celltype
JTK_result_all<-readRDS('~/analysis/circadian/R/JTK.16individual.cor2batch0.2.fold2bkg1.5.fold2trough1.5.rds')
JTK.result.filtered<-JTK_result_all[JTK_result_all$ADJ.P<=0.05,]
results<-enrichGObyHGNC(genes=unique(JTK.result.filtered$CycID))
clusterProfiler::dotplot(results)

celltypes<-names(CELL_TYPES)
#specific.genes<-getOnceOccured(JTK.result.filtered$CycID)
#JTK.result.filtered<-JTK.result.filtered[JTK.result.filtered$CycID %in% specific.genes,]
result_list<-list()
for(celltype in celltypes){
  celltype.specific.genes=JTK.result.filtered[JTK.result.filtered$celltype==celltype,"CycID"] %>% unique()
  #celltype.genes=FetchExpressedGenesCellType(test,celltype)
  results=enrichGObyHGNC(celltype.specific.genes)
  result_list[[celltype]]=results
}
saveRDS(result_list,"~/analysis/circadian/R/enrichGO.16individual.rds")

genes <- mapIds(org.Hs.eg.db,
                keys = "GO:0033151",
                column = "SYMBOL",
                keytype = "GO",
                multiVals = "list")[[1]]
genes

# analyse CYTOF data
library(CATALYST)
library(flowCore)
library(arrow)
srt<-feather2SeuratObject('~/data/CYTOF/exp_mat.feather','~/data/CYTOF/metadata_.feather')
srt[["RNA"]]$data <- srt[["RNA"]]$counts
srt<-RenameCells(srt,add.cell.id=srt@meta.data$sample)
srt<-RenameCells(srt,add.cell.id="TF")
srt$type<-srt$merge %>% gsub("\u00A0","_",.)
srt$CT<-getField(srt$sample,"_",2)
srt$individual<-getField(srt$sample,"_",1)
saveRDS(srt,'~/analysis/circadian/R/CYTOF.5individual.srt.rds')

srt<-readRDS('~/analysis/circadian/R/CYTOF.5individual.srt.rds')
srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 40)
srt <- ScaleData(srt, do.scale = TRUE, do.center = TRUE)
srt <- RunPCA(srt, features = VariableFeatures(object = srt))
srt <- RunUMAP(srt, dims = 1:20)
DimPlot(srt, reduction = "umap",group.by = "merge",cols = generateColor(25))




# find effective oscillators like LYZ
srt.metacell<-readRDS('~/analysis/circadian/R/seacell.16individual.rds')
JTK<-readRDS("~/analysis/circadian/R/JTK.result.filtered.16individual.addp2t.addfold2bkg.bytype.rds")
JTK.raw<-readRDS("~/analysis/circadian/R/JTK.result.filtered.16individual.addp2t.addfold2bkg.bytype.rds")
JTK<-JTK[JTK$fold_to_background>=2&!is.na(JTK$CycID)&JTK$PER<=28&JTK$PER>=20&JTK$fold_peak2trough>=2,]

JTK<-JTK %>% group_by(CycID,celltype) %>% mutate(n.individual=n(),mean_fold_peak2trough=mean(fold_peak2trough))
JTK[c("CycID","celltype","n.individual","mean_fold_peak2trough")] %>% unique() %>% dplyr::filter(n.individual>2) %>% arrange(desc(n.individual)) %>% View()

plotMetaCellByIndividual(srt.metacell,cell.type = "CD14 Mono",feature = "LYZ",merge = F) 

plotMetaCellByIndividual(srt.metacell,cell.type = "CD14 Mono",feature = "CCL3",merge = F) 


# using tauFisher to predict internal time
individual1.Mono<-read.delim('~/analysis/circadian/R/seacell.meta2d.0.04.12Healthy.bytype/JTKresult_HZD_CD14_Mono.txt')
individual2.Mono<-read.delim('~/analysis/circadian/R/seacell.meta2d.0.04.12Healthy.bytype/JTKresult_WLG_CD14_Mono.txt')
genes.used<-c("BMAL1","CLOCK","PER1","PER2","PER3","CRY1","CRY2","DBP")

predictIndividualwithTauFisher(matrix.train = '~/analysis/circadian/R/seacell.meta2d.0.04.12Healthy.bytype/HZD_CD14_Mono.txt', 
                               matrix.predict ='~/analysis/circadian/R/seacell.meta2d.0.04.12Healthy.bytype/LJQ_CD14_Mono.txt',
                               predictor.genes = genes.used)

JTK<-readRDS('~/analysis/circadian/R/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')

srt<-readRDS('~/analysis/circadian/R/16individual.srt.annotated.rds')
bulk_by_type<-readRDS('~/analysis/circadian/R/16individual.pseudobulk.bytype.rds')


JTK.filtered<-JTK[JTK$PER>=20&JTK$fold_to_background>2&JTK$fold_peak2trough>1.5&JTK$median_expression>10,] %>% arrange(ADJ.P)
JTK.filtered<-JTK.filtered[!is.na(JTK.filtered$CycID),]
JTK.filtered<-JTK.filtered %>% group_by(CycID,celltype) %>% mutate(n.individual=n()) %>% arrange(desc(n.individual))

JTK.used<-JTK.filtered[JTK.filtered$celltype=="CD4 TCM"&JTK.filtered$individual=="KD",]
genes.used<-(JTK.used %>% group_by(LAG) %>% top_n(.,wt=n.individual,1))$CycID
predictIndividualwithTauFisher(data=bulk_by_type,individual = "KD",celltype = "CD4_TCM",predictor.genes = genes.used)

JTK.used<-JTK.filtered[JTK.filtered$celltype=="CD14 Mono"&JTK.filtered$individual=="JJC"&JTK.filtered$PER==24,]
genes.used<-(JTK.used %>% group_by(LAG) %>% top_n(.,wt=n.individual,1))$CycID
predictIndividualwithTauFisher(data=bulk_by_type,individual = "KD",celltype = "CD14_Mono",predictor.genes = genes.used)
predictIndividualwithTauFisher(data=bulk_by_type,individual = "KD",celltype = "CD14_Mono",
                               predictor.genes = c("BMAL1","CLOCK","PER1","PER2","CRY1","CRY2","NR1D1","NR1D2"))

JTK.used<-JTK.filtered[JTK.filtered$celltype=="CD14 Mono"&JTK.filtered$individual%in%c("HZD","WLG"),] %>% group_by(CycID,celltype) %>% mutate(n.individual=n())
JTK.used<-JTK.used[JTK.used$n.individual==2,]
genes.used<-unique((JTK.used %>% group_by(LAG) %>% top_n(.,wt=desc(ADJ.P),1))$CycID)
predictIndividualwithTauFisher(matrix.train = '~/analysis/circadian/R/seacell.meta2d.0.04.12Healthy.bytype/HZD_CD14_Mono.txt', 
                               matrix.predict ='~/analysis/circadian/R/seacell.meta2d.0.04.12Healthy.bytype/HZD_CD16_Mono.txt',
                               predictor.genes = genes.used)

# TOP10 by number of individual
JTK_result_all<-readRDS('~/analysis/circadian/R/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
JTK_result_all<-JTK_result_all %>% group_by(CycID,celltype) %>% mutate(n.occurance=n())
JTK_result_filtered<-JTK_result_all[JTK_result_all$fold_peak2trough>=1.5&JTK_result_all$median_expression>10&JTK_result_all$fold_to_background>2,]
used.celltype<-"CD14 Mono"
JTK_result_used<-JTK_result_filtered[!is.na(JTK_result_filtered$CycID)&JTK_result_filtered$celltype==used.celltype&JTK_result_filtered$PER>=20,]
JTK_result_used<-JTK_result_used %>% group_by(CycID,celltype) %>% mutate(n.individual=n())
JTK_result_used<-JTK_result_used %>% group_by(individual) %>% mutate(percentile_rank=percent_rank(ADJ.P))
JTK_result_used<-JTK_result_used %>% group_by(CycID) %>% mutate(median.percentile=median(percentile_rank),median.Fpeak2trough=median(fold_peak2trough))
top.genes<-(JTK_result_used[c("CycID","n.individual")] %>% unique() %>% arrange(desc(n.individual)) %>% head(.,100))$"CycID"
top.genes<-(JTK_result_used[c("CycID","median.percentile")] %>% unique() %>% arrange(median.percentile) %>% head(.,50))$"CycID"
top.genes<-c("PER1","FMN1","NR1D2",'RNF144B','FKBP5','TSC22D3','MIR181A1HG','CEBPD')
predictIndividualwithTauFisher(data=bulk_by_type,individual = "JJC",celltype = "CD14_Mono",
                               predictor.genes = top.genes)

## interestingly, many individual's points go out of the FDA embeded PCA, it indicates that gene to gene difference belongs to multiple
## timing system, which lead to different phase shift between individuals

# top n by prevalence/median percentile of adjusted p
plotdata=NULL
for (this.individual in JTK_result_used$individual %>% unique()) {
  this.plotdata=data.frame("n.top.genes"=5:30,"accuracy"=NA,"RMSE"=NA)
  for (n.top.genes in 5:30) {
    message(n.top.genes)
    #top.genes=(JTK_result_used[c("CycID","n.individual")] %>% unique() %>% arrange(desc(n.individual)) %>% head(.,n.top.genes))$"CycID"
    top.genes=(JTK_result_used[c("CycID","median.percentile")] %>% unique() %>% arrange(median.percentile) %>% head(.,n.top.genes))$"CycID"
    res=predictIndividualwithTauFisher(data=bulk_by_type,individual = this.individual,celltype = "CD14_Mono",
                                       predictor.genes = top.genes,return.res = T)
    this.plotdata[this.plotdata$n.top.genes==n.top.genes,"accuracy"]=res$Accuracy
    this.plotdata[this.plotdata$n.top.genes==n.top.genes,"RMSE"]=res$RMSE
  }
  this.plotdata$individual=this.individual
  if(is.null(plotdata)){
    plotdata=this.plotdata
  }else{
    plotdata=rbind(plotdata,this.plotdata)
  }
}
plotdata.filtered<-plotdata[plotdata$n.top.genes %in% c(5,10,15,20,30),]
# Assuming your data is in a dataframe called 'df'
# First, reshape the data to long format for ggplot
df_long <- plotdata %>%
  pivot_longer(cols = c(accuracy, RMSE), 
               names_to = "metric", 
               values_to = "value")

# For separate panels
ggplot(df_long, aes(x = factor(n.top.genes), y = value, fill = metric)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    width = 0.6
  ) +
  geom_point(
    position = position_jitter(width = 0.2),
    size = 1.5,
    alpha = 0.6,
    aes(color = metric)
  ) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("accuracy" = "#2E86AB", "RMSE" = "#A23B72")) +
  scale_color_manual(values = c("accuracy" = "#2E86AB", "RMSE" = "#A23B72")) +
  labs(
    x = "Number of Top Genes",
    y = "Value"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(),
    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2)
  )

# number of individual as variable
plotdata=NULL
for (n.individual in 1:8) {
  for (rep.id in 1:3) {
    #message(n.top.genes)
    test.individuals=sample(unique(JTK_result_used$individual),n.individual,replace=F)
    top.genes=(JTK_result_used[c("CycID","n.individual")] %>% unique() %>% arrange(desc(n.individual)) %>% head(.,10))$"CycID"
    #top.genes=(JTK_result_used[c("CycID","median.percentile")] %>% unique() %>% arrange(median.percentile) %>% head(.,n.top.genes))$"CycID"
    #top.genes=c("BMAL1","CLOCK","PER1","PER2","CRY1","CRY2","NR1D1","NR1D2","DBP")
    top.genes=c(top.genes,c("BMAL1","CLOCK","PER1","PER2","CRY1","CRY2","NR1D1","NR1D2","DBP")) %>% unique()
    res=predictIndividualwithTauFisher(data=bulk_by_type,individual = test.individuals,celltype = "CD14_Mono",
                                       predictor.genes = top.genes,return.res = T)
    this.plotdata=data.frame("n.individual"=n.individual,"accuracy"=res$Accuracy,"RMSE"=res$RMSE,"replicate"=rep.id)
    if(is.null(plotdata)){
      plotdata=this.plotdata
    }else{
      plotdata=rbind(plotdata,this.plotdata)
    }
  }
}
predictIndividualwithTauFisher(data=bulk_by_type,individual = c("KD","JJC","HZD","ZYX","WLG","XSP"),celltype = "CD14_Mono",
                               predictor.genes = top.genes)

df_long <- plotdata %>%
  pivot_longer(cols = c(accuracy, RMSE), 
               names_to = "metric", 
               values_to = "value")

# For separate panels
ggplot(df_long, aes(x = factor(n.individual), y = value, fill = metric)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    width = 0.6
  ) +
  geom_point(
    position = position_jitter(width = 0.2),
    size = 1.5,
    alpha = 0.6,
    aes(color = metric)
  ) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("accuracy" = "#2E86AB", "RMSE" = "#A23B72")) +
  scale_color_manual(values = c("accuracy" = "#2E86AB", "RMSE" = "#A23B72")) +
  labs(
    x = "Number of Individual used",
    y = "Value"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(),
    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2)
  )

# calculate dispersion of inter-individual phase shift between genes within individual
plotdata=NULL
for (this.individual in unique(JTK_result_used$individual)) {
  message(this.individual)
  this.JTK=JTK_result_used[JTK_result_used$individual==this.individual,c("CycID","LAG")]
  for(this.gene in this.JTK$CycID){
    this.lag=this.JTK[this.JTK$CycID==this.gene,"LAG"] %>% as.numeric()
    result=this.JTK
    result$CycID2=this.gene
    result$LAG2=this.lag
    result$lag_diff=result$LAG-result$LAG2
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
saveRDS(plotdata,paste0('~/analysis/circadian/R/lag.diff.',cell.type.str,".rds"))

plotdata<-readRDS('~/analysis/circadian/R/lag.diff.CD8_TEM.rds')
plotdata.filtered<-plotdata[plotdata$n.pair>5,]
plotdata.filtered<-plotdata.filtered %>% group_by(CycID,CycID2) %>% mutate(sd.lag.diff=sd(lag_diff))
plotdata.filtered<-plotdata.filtered[c("CycID","CycID2","sd.lag.diff")] %>% unique()
plotdata.mat<-spread(plotdata.filtered,key="CycID2",value="sd.lag.diff")
plotdata.mat<-as.data.frame(plotdata.mat)
plotdata.mat<-column_to_rownames(plotdata.mat,"CycID")

df<-plotdata.filtered
# First, create the distance matrix as before
all_cycids <- unique(c(df$CycID, df$CycID2))
dist_matrix <- matrix(NA, nrow = length(all_cycids), ncol = length(all_cycids),
                      dimnames = list(all_cycids, all_cycids))

for(i in 1:nrow(df)) {
  dist_matrix[df$CycID[i], df$CycID2[i]] <- df$sd.lag.diff[i]
  dist_matrix[df$CycID2[i], df$CycID[i]] <- df$sd.lag.diff[i]
}
diag(dist_matrix) <- 0

# Impute missing values with the mean distance
mean_distance <- mean(dist_matrix, na.rm = TRUE)
dist_matrix_imputed <- dist_matrix
dist_matrix_imputed[is.na(dist_matrix_imputed)] <- mean_distance

# Or use a more sophisticated imputation: maximum distance + small buffer
max_distance <- max(dist_matrix, na.rm = TRUE)
dist_matrix_imputed <- dist_matrix
dist_matrix_imputed[is.na(dist_matrix_imputed)] <- max_distance * 1.1

# Now create the heatmap with imputed matrix
col_fun <- colorRamp2(c(0, max(df$sd.lag.diff, na.rm = TRUE)), 
                      c("#003366", "#FFEE99"))

ht <- Heatmap(
  dist_matrix_imputed,
  name = "Distance",
  col = col_fun,
  na_col = "white",
  
  # Use the imputed matrix for clustering
  clustering_distance_rows = function(x) as.dist(dist_matrix_imputed),
  clustering_distance_columns = function(x) as.dist(dist_matrix_imputed),
  
  # CNS-level styling
  rect_gp = gpar(col = "white", lwd = 0.1),
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

# For PDF (vector format, best for publications)
svg("distance_heatmap.svg", 
    width = 25, 
    height = 25, 
    pointsize = 8) # CNS standard: 8-10 pt font

draw(ht) # Your Heatmap object
dev.off() # Essential: closes the graphics device

# Define your genes of interest
genes_of_interest <- c("PER1", "NR1D2", "BMAL1", "DDIT4", "FKBP5","PER2","PER3","CRY1","CRY2","NR1D1",'CLOCK') # Replace with your actual IDs

# Create the base heatmap without row/column labels
ht <- Heatmap(
  dist_matrix_imputed,
  name = "Phase shift\ndeviation",
  col = col_fun,
  na_col = "white",
  clustering_distance_rows = function(x) as.dist(dist_matrix_imputed),
  clustering_distance_columns = function(x) as.dist(dist_matrix_imputed),
  rect_gp = gpar(col = "white", lwd = 0.1),
  show_row_names = FALSE,    # Turn off all row labels
  show_column_names = FALSE, # Turn off all column labels
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

# Add text annotations for specific genes
# Note: This requires knowing the positions after clustering
# You might need to run the clustering first to get positions

# Alternative: Use decorate_annotation for more control
ht <- ht + rowAnnotation(
  foo = anno_mark(
    at = which(rownames(dist_matrix_imputed) %in% genes_of_interest),
    labels = rownames(dist_matrix_imputed)[rownames(dist_matrix_imputed) %in% genes_of_interest],
    labels_gp = gpar(fontsize = 8),
    link_gp = gpar(lwd = 0.5),
    padding = unit(1, "mm")
  )
)
draw(ht) # Your Heatmap object

# analyse circadian from CYTOF data
srt.cytof<-readRDS('~/analysis/circadian/R/CYTOF.300k.reannotation.srt.rds')
srt.metacell<-readRDS('~/analysis/circadian/R/seacell.16individual.rds')

JTK.metacell<-readRDS('~/analysis/circadian/R/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
JTK.cytof<-readJTKFromMetaCells('~/analysis/circadian/R/CYTOF.meta2d.by.newtype/')
JTK.cytof<-JTK.cytof[JTK.cytof$PER>=20&JTK.cytof$PER<=28,]
saveRDS(JTK.cytof,"~/analysis/circadian/R/CYTOF.JTK.newresult.filtered.rds")
JTK.cytof<-readRDS('~/analysis/circadian/R/CYTOF.JTK.newresult.filtered.addmedianexp.addp2t.bytype.rds')
JTK.cytof$CycID<-CYTOF2RNA_FEATURES[JTK.cytof$CycID]
#JTK.cytof$celltype.raw<-JTK.cytof$celltype
#JTK.cytof$celltype<-CYTOF2RNA_CELLTYPES[JTK.cytof$celltype.raw]
JTK.cytof<-JTK.cytof[JTK.cytof$median_expression>0.1&JTK.cytof$fold_peak2trough>1,]
JTK.cytof<-JTK.cytof[!is.na(JTK.cytof$CycID),]
#srt.cytof@meta.data$RNA_type<-CYTOF2RNA_CELLTYPES[srt.cytof@meta.data$type]
#srt.cytof$type_raw<-srt.cytof$type
srt.cytof$type<-srt.cytof$manual.new

plotCytofbyIndividual(srt.cytof,cell.type = "CD8 TEM",features = "CXCR4",proportion = 1)

comparePhase<-function(JTK.result1,JTK.result2,feature){
  message(feature)
  individuals.used=intersect(JTK.result1$individual,JTK.result2$individual) %>% unique()
  features.used=intersect(JTK.result1$CycID,JTK.result2$CycID) %>% unique()
  JTK.result1=JTK.result1[JTK.result1$CycID %in% features.used,c("CycID","LAG","celltype","individual")]
  JTK.result2=JTK.result2[JTK.result2$CycID %in% features.used,c("CycID","LAG","celltype","individual")]
  JTK.result1$method="protein"
  JTK.result2$method="mRNA"
  JTK.result.all=rbind(JTK.result1,JTK.result2)
  JTK.result.all=JTK.result.all[JTK.result.all$CycID==feature&JTK.result.all$individual %in% individuals.used,]
  mRNA_LAG=JTK.result.all[JTK.result.all$method=="mRNA",] %>% group_by(individual) %>% summarise(median_LAG=median(LAG))
  JTK.result.all=left_join(JTK.result.all,mRNA_LAG,by="individual")
  JTK.result.all$relative_LAG=JTK.result.all$LAG-JTK.result.all$median_LAG
  print(JTK.result.all)
  JTK.result.all[JTK.result.all$method=="protein","relative_LAG"]=ifelse(JTK.result.all[JTK.result.all$method=="protein","relative_LAG"]<0,
                                                                         JTK.result.all[JTK.result.all$method=="protein","relative_LAG"]+24,
                                                                         JTK.result.all[JTK.result.all$method=="protein","relative_LAG"])
  JTK.result.all[JTK.result.all$method=="mRNA","relative_LAG"]=ifelse(JTK.result.all[JTK.result.all$method=="mRNA","relative_LAG"]< -12,
                                                                         JTK.result.all[JTK.result.all$method=="mRNA","relative_LAG"]+24,
                                                                         JTK.result.all[JTK.result.all$method=="mRNA","relative_LAG"])
  #print(JTK.result.all)
  ggplot(JTK.result.all, aes(x = method, y = relative_LAG, fill = individual)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.7,outlier.alpha = 0) +
    labs(title = "Time since mRNA peak",
         x = "Method",
         y = "Peaking Time (LAG)",
         fill = "Individual") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3")
}
comparePhase(JTK.cytof,JTK.metacell,"CXCR4")
comparePhase(JTK.cytof,JTK.metacell,"CD3D")

for (feature in unique(JTK.cytof$CycID)){
  comparePhase(JTK.cytof,JTK.metacell,feature)
  ggsave(paste0('~/analysis/circadian/R/CYTOFvsmRNA/',feature,'.pdf'),width = 4,height=3)
}

plotCytofbyIndividual(srt.cytof,cell.type = "CD14 Mono",features = "CXCR4",proportion = 0.10)

# merge velocity data from loom
test.mat<-getWholeDaySplicedData('HZD')
test.mat$CT09$spliced[1:3,1:3]

# reclustering and annotation of CYTOF
test.cytof<-readRDS('~/analysis/circadian/R/CYTOF.5individual.srt.rds')
test.cytof<-ScaleData(test.cytof)
test.cytof<-RunPCA(test.cytof,features=Features(test.cytof))
test.cytof.subset<-subset(test.cytof,cells=sample(Cells(test.cytof),300000,replace = F))
DimPlot(test.cytof)
test.cytof.subset<-FindNeighbors(test.cytof.subset)
test.cytof.subset<-FindClusters(test.cytof.subset,resolution = 2)
test.cytof.subset<-RunUMAP(test.cytof.subset,dims = 1:20)
DimPlot(test.cytof.subset,reduction = "umap",label = T)
markers<-FindAllMarkers(test.cytof.subset)
markers$HGNC<-CYTOF2RNA_FEATURES[markers$gene]
markers %>% group_by(cluster) %>% top_n(10,pct.1) %>% dplyr::filter(.,cluster==2)
FeaturePlot(test.cytof.subset,c("CD45RA","CD4","CD8a","TCR-g-d-plt","CD19-plt","CD57",
                                "CD14-plt","CD25","CD56-plt","IgD-plt","CD123-plt","CD16-plt"),order = T)

# annotate myloids
cytof.myloids<-subset(test.cytof.subset,subset=seurat_clusters %in% c(8,22,26,29,37))
cytof.myloids<-dimentionalReductionSubsetSrtCYTOF(cytof.myloids)
DimPlot(cytof.myloids,reduction = "umap",label = T)
FeaturePlot(cytof.myloids,c("CD14-plt","CD16-plt","CD123-plt","HLA-DR-plt"))
cytof.myloids<-TypeCluster(cytof.myloids,c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17),type="CD14 Mono",new.meta="manual")
cytof.myloids<-TypeCluster(cytof.myloids,c(0),type="CD16 Mono",new.meta="manual")
cytof.myloids<-TypeCluster(cytof.myloids,c(3),type="cDC2",new.meta="manual")
cytof.myloids<-TypeCluster(cytof.myloids,c(11),type="pDC",new.meta="manual")
DimPlot(cytof.myloids,reduction = "umap",label = T,group.by = "manual")
saveRDS(cytof.myloids,"~/analysis/circadian/R/CYTOF.myeloid.srt.rds")

# annotate B
cytof.b<-subset(test.cytof.subset,subset=seurat_clusters %in% c(12,33,36,41))
cytof.b<-dimentionalReductionSubsetSrtCYTOF(cytof.b)
DimPlot(cytof.b,reduction = "umap",label = T)
FeaturePlot(cytof.b,c("IgD-plt","CD27",'CD38-plt','CD19-plt'))
cytof.b<-TypeCluster(cytof.b,c(1,5,10,12,15),type="B naive",new.meta="manual")
cytof.b<-TypeCluster(cytof.b,c(0,2,4,7,8,9,11,14),type="B memory",new.meta="manual")
cytof.b<-TypeCluster(cytof.b,c(6),type="Plasmablast",new.meta="manual")
cytof.b<-TypeCluster(cytof.b,c(3,13),type="doublet",new.meta="manual")
DimPlot(cytof.b,reduction = "umap",label = T,group.by = "manual")
saveRDS(cytof.b,"~/analysis/circadian/R/CYTOF.b.srt.rds")

# annotate NK
cytof.nk<-subset(test.cytof.subset,subset=seurat_clusters %in% c(0,1,10,34,40,44))
cytof.nk<-dimentionalReductionSubsetSrtCYTOF(cytof.nk)
DimPlot(cytof.nk,reduction = "umap",label = T)
FeaturePlot(cytof.nk,c("CD56-plt","CD16-plt"))
cytof.nk<-TypeCluster(cytof.nk,c(0:6,8:15),type="NK",new.meta="manual")
cytof.nk<-TypeCluster(cytof.nk,c(7),type="NK_CD56ᵇʳⁱᵍʰᵗ",new.meta="manual")
DimPlot(cytof.nk,reduction = "umap",label = T,group.by = "manual")
saveRDS(cytof.nk,"~/analysis/circadian/R/CYTOF.nk.srt.rds")

cytof.t<-subset(test.cytof.subset,subset=seurat_clusters %in% c(7,45,32,47,39,14,42,38,27,16,5,
                                                                4,31,25,20,24,43,18,11,17,3,6,2,21,
                                                                23,15,19,13,46,9,28,35))
cytof.t<-dimentionalReductionSubsetSrtCYTOF(cytof.t,res = 2)

FeaturePlot(cytof.t,c("CD27","CD4","CD8a","CD45RA","CD57","CD197-plt","CD25","TCR-g-d-plt","CX3CR1"))
DimPlot(cytof.t,reduction = "umap",label = T)

cytof.t<-TypeCluster(cytof.t,c(6,11,13,19,29,32,33,44),type="CD4 Naive",new.meta="manual")
cytof.t<-TypeCluster(cytof.t,c(0,1,41,2,5,12,14,17,21,26,31,46),type="CD4 TCM",new.meta="manual")
cytof.t<-TypeCluster(cytof.t,c(7,11,36),type="doublet",new.meta="manual")
cytof.t<-TypeCluster(cytof.t,c(35),type="CD4 CTL",new.meta="manual")
cytof.t<-TypeCluster(cytof.t,c(27,43),type="Treg",new.meta="manual")
cytof.t<-TypeCluster(cytof.t,c(10,22,24,25,28,30,40,45),type="CD8 Naive",new.meta="manual")
cytof.t<-TypeCluster(cytof.t,c(20),type="CD8 TCM",new.meta="manual")
cytof.t<-TypeCluster(cytof.t,c(3,18,23),type="MAIT",new.meta="manual")
cytof.t<-TypeCluster(cytof.t,c(8,9,15,37,38,42,48,50,47),type="CD8 TEM",new.meta="manual")
cytof.t<-TypeCluster(cytof.t,c(39),type="dnT",new.meta="manual")
cytof.t<-TypeCluster(cytof.t,c(4,16,34,49),type="gdT",new.meta="manual")
DimPlot(cytof.t,reduction = "umap",label = T,group.by = "manual")
cytof.t$manual %>% table()
saveRDS(cytof.t,"~/analysis/circadian/R/CYTOF.t.srt.rds")

generateAnnotationFile('~/analysis/circadian/R/CYTOF.b.srt.rds',
                       '~/analysis/circadian/R/CYTOF.myeloid.srt.rds',
                       '~/analysis/circadian/R/CYTOF.nk.srt.rds',
                       '~/analysis/circadian/R/CYTOF.t.srt.rds',
                       cols = c("type","individual","CT","sample","manual"),
                       out.path = '~/analysis/circadian/R/CYTOF.annotation.manual.tsv')

new.cytof.annotation<-read.delim('~/analysis/circadian/R/CYTOF.annotation.manual.tsv')
annotation<-new.cytof.annotation$manual
names(annotation)<-new.cytof.annotation$cell_id
test.cytof.subset<-AddMetaData(test.cytof.subset,metadata=annotation,col.name = "manual.new")
DimPlot(test.cytof.subset,reduction = "umap",label = T,group.by = "manual.new")
test.cytof.subset.subset<-subset(test.cytof.subset,manual.new!="doublet")
DimPlot(test.cytof.subset.subset,reduction = "umap",label = T,group.by = "manual.new")
saveRDS(test.cytof.subset.subset,"~/analysis/circadian/R/CYTOF.300k.reannotation.srt.rds")

# analyse mouse skin data from tauFisher group
setwd('/dssg/home/acct-medll/medll/data/mouse_skin_circadian')
srts<-list()
for (sample in list.dirs(full.names = F)[-1]){
  this.mtx=Read10X(sample)
  this.srt=CreateSeuratObject(this.mtx)
  this.srt$sample=sample
  this.srt=RenameCells(this.srt,add.cell.id=sample)
  srts[[sample]]=this.srt
}
srts<-merge(srts[[1]],srts[2:length(srts)])
srts[['RNA']]<-JoinLayers(srts[['RNA']])
srts<-dimentionalReductionMergedSrt(srts)
srts<-integrateMergedSrt(srts)
srts@meta.data$ZT<-getField(srts@meta.data$sample,"_",2)

DimPlot(srts)
FeaturePlot(srts,c("Ptprc","Cd14","Fcgr3","Cd3d","Vim","Ca12"))
saveRDS(srts,'~/analysis/circadian/R/mouseSkin.srt.integrated.rds')

srts<-FindClusters(srts, resolution = 0.5)
DimPlot(srts,label = T)
markers.table<-FindAllMarkers(srts,min.pct = 0.7,min.diff.pct = 0.5,logfc.threshold = log2(2))
srts<-TypeCluster(srts,cluster = c(0,2,3,4,5,6,13,15),type = "Fibroblast")
srts<-TypeCluster(srts,cluster = c(1,8,12,16),type = "Erythroid cell")
srts<-TypeCluster(srts,cluster = c(7,10,11,17),type = "Macrophage")
srts<-TypeCluster(srts,cluster = c(19),type = "Monocyte")
srts<-TypeCluster(srts,cluster = c(20),type = "Monocyte proliferating")
srts<-TypeCluster(srts,cluster = c(9,14),type = "Muscle cell")
srts<-TypeCluster(srts,cluster = c(18),type = "TNK cell")
srts<-TypeCluster(srts,cluster = c(22),type = "Mast cell")
srts<-TypeCluster(srts,cluster = c(21),type = "debris")
DimPlot(srts,label = T,group.by = "type")
saveRDS(srts,'~/analysis/circadian/R/mouseSkin.srt.integrated.annotated.rds')

srts<-readRDS('~/analysis/circadian/R/mouseSkin.srt.integrated.annotated.rds')

srts<-NormalizeData(srts,normalization.method = "RC",scale.factor = 1000000)
srts@meta.data$individual<-getField(srts@meta.data$sample,"_",1)
mouse_skin_pseudobulk<-srt2bulkMatrix(srts,c("individual","ZT","type"),layer = "data")
mouse_skin_pseudobulk<-mouse_skin_pseudobulk[as.vector(colSums(mouse_skin_pseudobulk)!=0)]
# Calculate scaling factors for each sample
scaling_factors<-colSums(mouse_skin_pseudobulk)/1000000
# Normalize the data
normalized_mouse_skin_pseudobulk <- sweep(mouse_skin_pseudobulk, 2, scaling_factors, "/")
normalized_mouse_skin_pseudobulk<-rownames_to_column(normalized_mouse_skin_pseudobulk,var = "feature")
saveRDS(normalized_mouse_skin_pseudobulk,'~/analysis/circadian/R/mouseSkin.psedobulk.roughType.CPM.rds')
normalized_mouse_skin_pseudobulk<-readRDS('~/analysis/circadian/R/mouseSkin.psedobulk.roughType.CPM.rds')
normalized_mouse_skin_pseudobulk_long<-gather(normalized_mouse_skin_pseudobulk,key="observation",value="CPM",-feature)
normalized_mouse_skin_pseudobulk_long$ZT<-getField(normalized_mouse_skin_pseudobulk_long$observation,sep = "-",field = 2)
normalized_mouse_skin_pseudobulk_long$type<-getField(normalized_mouse_skin_pseudobulk_long$observation,sep = "-",field = 3)

interest_gene="Ddit4"
plotdata<-normalized_mouse_skin_pseudobulk_long[normalized_mouse_skin_pseudobulk_long$feature==interest_gene,]
ggplot(plotdata)+geom_boxplot(aes(x=ZT,y=CPM,color=type))+facet_wrap(~type,scales = "free")+
  scale_x_discrete(limits=c("ZT2","ZT6","ZT10","ZT14","ZT18","ZT22"))+NoLegend()+ggtitle(interest_gene)

interest_gene=c("Per2","Arntl")
interest_gene=c("Per2","Cry1")
interest_gene=c("Per2","Nr1d2")
interest_gene=c("Per2","Dbp")
plotdata<-normalized_mouse_skin_pseudobulk_long[normalized_mouse_skin_pseudobulk_long$feature%in%interest_gene,]
plotdata<-plotdata %>% group_by(feature,type,ZT) %>% summarise(median_CPM=median(CPM))
plotdata<-plotdata %>% group_by(feature,type) %>% mutate(z_score=(median_CPM-mean(median_CPM))/sd(median_CPM))
ggplot(plotdata)+geom_point(aes(x=ZT,y=z_score,color=feature))+
  geom_line(aes(x=ZT,y=z_score,group=feature,color=feature))+facet_wrap(~type)+
  scale_x_discrete(limits=c("ZT2","ZT6","ZT10","ZT14","ZT18","ZT22"))+
  ggtitle(paste0(interest_gene,collapse = "_"))+theme(axis.text.x = element_text(angle=60,hjust=1))

ggplot(plotdata)+geom_point(aes(x=ZT,y=median_CPM,color=feature))+
  geom_line(aes(x=ZT,y=median_CPM,group=feature,color=feature))+facet_wrap(~type,scales = "free_y")+
  scale_x_discrete(limits=c("ZT2","ZT6","ZT10","ZT14","ZT18","ZT22"))+
  ggtitle(paste0(interest_gene,collapse = "_"))+theme(axis.text.x = element_text(angle=60,hjust=1))

# View other metadata
meta.plotdata<-srts@meta.data %>% group_by(ZT,type) %>% summarise(n.cell=n())
ggplot(meta.plotdata)+
  geom_bar(aes(x=ZT,y=n.cell),stat="identity")+
  facet_wrap(~type,scales = "free_y")+theme(axis.text.x=element_text(angle=60,hjust=1))+
  ggtitle("number of cells")+scale_x_discrete(limits=c("ZT2","ZT6","ZT10","ZT14","ZT18","ZT22"))

meta.plotdata<-srts@meta.data %>% group_by(ZT,type) %>% summarise(median.n.feature=median(nFeature_RNA))
ggplot(meta.plotdata)+
  geom_bar(aes(x=ZT,y=median.n.feature),stat="identity")+
  facet_wrap(~type,scales = "free_y")+theme(axis.text.x=element_text(angle=60,hjust=1))+
  ggtitle("median nFeature")+scale_x_discrete(limits=c("ZT2","ZT6","ZT10","ZT14","ZT18","ZT22"))

meta.plotdata<-srts@meta.data %>% group_by(ZT,type) %>% summarise(median.n.count=median(nCount_RNA))
ggplot(meta.plotdata)+
  geom_bar(aes(x=ZT,y=median.n.count),stat="identity")+
  facet_wrap(~type,scales = "free_y")+theme(axis.text.x=element_text(angle=60,hjust=1))+
  ggtitle("median nCount")+scale_x_discrete(limits=c("ZT2","ZT6","ZT10","ZT14","ZT18","ZT22"))

# View raw count
mouse_skin_pseudobulk<-readRDS('~/analysis/circadian/R/mouseSkin.psedobulk.roughType.rds')
mouse_skin_pseudobulk<-rownames_to_column(mouse_skin_pseudobulk,var = "feature")
mouse_skin_pseudobulk_long<-gather(mouse_skin_pseudobulk,key="observation",value="raw_count",-feature)
mouse_skin_pseudobulk_long$individual<-getField(mouse_skin_pseudobulk_long$observation,sep = "-",field = 1)
mouse_skin_pseudobulk_long$ZT<-getField(mouse_skin_pseudobulk_long$observation,sep = "-",field = 2)
mouse_skin_pseudobulk_long$type<-getField(mouse_skin_pseudobulk_long$observation,sep = "-",field = 3)
interest_gene<-c("Per2","Arntl")
plotdata2<-mouse_skin_pseudobulk_long[mouse_skin_pseudobulk_long$feature %in% interest_gene,]
plotdata2<-plotdata2 %>% group_by(feature,type,ZT) %>% summarise(median_raw_count=median(raw_count))

ggplot(plotdata2)+geom_point(aes(x=ZT,y=median_raw_count,color=feature))+
  geom_line(aes(x=ZT,y=median_raw_count,color=feature,group=feature))+
  facet_wrap(~type,scales = "free_y")+
  ggtitle(paste0(paste(interest_gene,collapse = "_")))+
  theme(axis.text.x=element_text(angle=60,hjust=1))+
  scale_x_discrete(limits=c("ZT2","ZT6","ZT10","ZT14","ZT18","ZT22"))

### human data
test<-readRDS('~/analysis/circadian/R/16individual.srt.annotated.rds')
test<-NormalizeData(test,normalization.method = "RC",scale.factor = 1000000)
test$type<-test$predicted.celltype.l2
human_pbmc_pseudobulk<-srt2bulkMatrix(test,c("individual","CT","type"),layer = "data")

scaling_factors<-colSums(human_pbmc_pseudobulk)/1000000
normalized_human_pbmc_pseudobulk<-sweep(human_pbmc_pseudobulk, 2, scaling_factors, "/")
saveRDS(normalized_human_pbmc_pseudobulk,'~/analysis/circadian/R/16individual.pseudobulk.bytype.CPM.rds')
normalized_human_pbmc_pseudobulk<-readRDS('~/analysis/circadian/R/16individual.pseudobulk.bytype.CPM.rds')

normalized_human_pbmc_pseudobulk<-rownames_to_column(normalized_human_pbmc_pseudobulk,var = "feature")
normalized_human_pbmc_pseudobulk_long<-gather(normalized_human_pbmc_pseudobulk,key="observation",value="CPM",-feature)
normalized_human_pbmc_pseudobulk_long$CT<-getField(normalized_human_pbmc_pseudobulk_long$observation,sep = "-",field = 2)
normalized_human_pbmc_pseudobulk_long$type<-getField(normalized_human_pbmc_pseudobulk_long$observation,sep = "-",field = 3)
plotdata2<-normalized_human_pbmc_pseudobulk_long[normalized_human_pbmc_pseudobulk_long$feature=="DDIT4",]
ggplot(plotdata2)+geom_boxplot(aes(x=CT,y=CPM),outliers = F)+
  scale_x_discrete(limits=c("CT09","CT13","CT17","CT21","CT25","CT29"))+facet_wrap(~type,scales = "free")
ggplot(plotdata2[plotdata2$type=="CD14_Mono",])+
  geom_boxplot(aes(x=CT,y=CPM),outliers = F)+scale_x_discrete(limits=c("CT11","CT15","CT19","CT23","CT27","CT31","CT35"))


# use a rough annotation
test@meta.data$type<-CELL_TYPES2[test@meta.data$predicted.celltype.l2]
test<-NormalizeData(test,normalization.method = "RC",scale.factor = 1000000)
human_pbmc_pseudobulk<-srt2bulkMatrix(test,c("individual","CT","type"),layer = "data")
scaling_factors<-colSums(human_pbmc_pseudobulk)/1000000
normalized_human_pbmc_pseudobulk<-sweep(human_pbmc_pseudobulk, 2, scaling_factors, "/")
saveRDS(normalized_human_pbmc_pseudobulk,'~/analysis/circadian/R/16individual.pseudobulk.byroughtype.CPM.rds')

normalized_human_pbmc_pseudobulk<-readRDS('~/analysis/circadian/R/16individual.pseudobulk.byroughtype.CPM.rds')

normalized_human_pbmc_pseudobulk<-rownames_to_column(normalized_human_pbmc_pseudobulk,var = "feature")
normalized_human_pbmc_pseudobulk_long<-gather(normalized_human_pbmc_pseudobulk,key="observation",value="CPM",-feature)
normalized_human_pbmc_pseudobulk_long$individual<-getField(normalized_human_pbmc_pseudobulk_long$observation,sep = "-",field = 1)
normalized_human_pbmc_pseudobulk_long$CT<-getField(normalized_human_pbmc_pseudobulk_long$observation,sep = "-",field = 2)
normalized_human_pbmc_pseudobulk_long$type<-getField(normalized_human_pbmc_pseudobulk_long$observation,sep = "-",field = 3)

interest_gene<-"CRY1"
plotdata2<-normalized_human_pbmc_pseudobulk_long[normalized_human_pbmc_pseudobulk_long$feature==interest_gene,]
plotdata2<-plotdata2 %>% group_by(type,individual) %>% mutate(z_score=(log2(CPM)-mean(log2(1+CPM)))/sd(log2(1+CPM)))

# y-axis: CPM
ggplot(plotdata2)+geom_boxplot(aes(x=CT,y=CPM),outliers = F)+
  scale_x_discrete(limits=c("CT09","CT13","CT17","CT21","CT25","CT29"))+facet_wrap(~type,scales = "free")+
  NoLegend()+ggtitle(interest_gene)

ggplot(plotdata2)+geom_boxplot(aes(x=CT,y=CPM),outliers = F)+
  scale_x_discrete(limits=CT_TIME_ORDER_FINE)+facet_wrap(~type,scales = "free")+
  NoLegend()+ggtitle(interest_gene)

ggplot(plotdata2)+geom_point(aes(x=CT,y=CPM,color=individual))+
  geom_line(aes(x=CT,y=CPM,color=individual,group=individual))+
  scale_x_discrete(limits=CT_TIME_ORDER_FINE)+facet_wrap(~type,scales = "free")+
  scale_color_d3("category20")+NoLegend()+ggtitle(interest_gene)

ggplot(plotdata2[plotdata2$type=="Monocyte",])+geom_point(aes(x=CT,y=CPM,color=individual))+
  geom_line(aes(x=CT,y=CPM,color=individual,group=individual))+
  scale_x_discrete(limits=CT_TIME_ORDER_FINE)+facet_wrap(~individual,scales = "free")+
  scale_color_d3("category20")+NoLegend()+ggtitle(paste0(interest_gene,"-NK"))

# y-axis: z-score
ggplot(plotdata2)+geom_boxplot(aes(x=CT,y=z_score),outliers = F)+
  scale_x_discrete(limits=CT_TIME_ORDER_FINE)+facet_wrap(~type,scales = "free")+
  NoLegend()+ggtitle(interest_gene)

ggplot(plotdata2)+geom_point(aes(x=CT,y=z_score,color=individual))+
  geom_line(aes(x=CT,y=z_score,color=individual,group=individual))+
  scale_x_discrete(limits=CT_TIME_ORDER_FINE)+facet_wrap(~type,scales = "free")+
  scale_color_d3("category20")+NoLegend()+ggtitle(interest_gene)

ggplot(plotdata2[plotdata2$type=="Monocyte",])+geom_point(aes(x=CT,y=z_score,color=individual))+
  geom_line(aes(x=CT,y=z_score,color=individual,group=individual))+
  scale_x_discrete(limits=CT_TIME_ORDER_FINE)+facet_wrap(~individual,scales = "free")+
  scale_color_d3("category20")+NoLegend()+ggtitle(paste0(interest_gene,"-Monocyte"))

# put significant circadian gene in same plot

normalized_human_pbmc_pseudobulk<-readRDS('~/analysis/circadian/R/16individual.pseudobulk.byroughtype.CPM.rds')

normalized_human_pbmc_pseudobulk<-rownames_to_column(normalized_human_pbmc_pseudobulk,var = "feature")
normalized_human_pbmc_pseudobulk_long<-gather(normalized_human_pbmc_pseudobulk,key="observation",value="CPM",-feature)
normalized_human_pbmc_pseudobulk_long$individual<-getField(normalized_human_pbmc_pseudobulk_long$observation,sep = "-",field = 1)
normalized_human_pbmc_pseudobulk_long$CT<-getField(normalized_human_pbmc_pseudobulk_long$observation,sep = "-",field = 2)
normalized_human_pbmc_pseudobulk_long$type<-getField(normalized_human_pbmc_pseudobulk_long$observation,sep = "-",field = 3)

interest_gene<-c("PER2","NR1D2","BMAL1","DBP","CRY1")
interest_gene<-c("PER2","CRY1")
interest_gene<-c("PER2","DBP")
interest_gene<-c("PER2","NR1D2")
interest_gene<-c("PER2","BMAL1")
plotdata2<-normalized_human_pbmc_pseudobulk_long[normalized_human_pbmc_pseudobulk_long$feature %in% interest_gene,]
plotdata2<-plotdata2 %>% group_by(type,individual,feature) %>% mutate(z_score=(log2(CPM)-mean(log2(1+CPM)))/sd(log2(1+CPM)))

ggplot(plotdata2[plotdata2$type=="Monocyte",])+geom_point(aes(x=CT,y=z_score,color=feature))+
  geom_line(aes(x=CT,y=z_score,color=feature,group=feature))+
  scale_x_discrete(limits=CT_TIME_ORDER_FINE)+facet_wrap(~individual)+
  scale_color_d3("category20")+ggtitle(paste0(paste(interest_gene,collapse = "_"),"-Monocyte"))+
  theme(axis.text.x=element_text(angle=60,hjust=1))

ggplot(plotdata2[plotdata2$type=="Monocyte",])+geom_point(aes(x=CT,y=CPM,color=feature))+
  geom_line(aes(x=CT,y=CPM,color=feature,group=feature))+
  scale_x_discrete(limits=CT_TIME_ORDER_FINE)+facet_wrap(~individual,scales = "free_y")+
  scale_color_d3("category20")+ggtitle(paste0(paste(interest_gene,collapse = "_"),"-Monocyte"))+
  theme(axis.text.x=element_text(angle=60,hjust=1))

# View raw count
normalized_human_pbmc_pseudobulk<-readRDS('~/analysis/circadian/R/16individual.pseudobulk.byroughtype.rds')
normalized_human_pbmc_pseudobulk<-rownames_to_column(normalized_human_pbmc_pseudobulk,var = "feature")
normalized_human_pbmc_pseudobulk_long<-gather(normalized_human_pbmc_pseudobulk,key="observation",value="raw_count",-feature)
normalized_human_pbmc_pseudobulk_long$individual<-getField(normalized_human_pbmc_pseudobulk_long$observation,sep = "-",field = 1)
normalized_human_pbmc_pseudobulk_long$CT<-getField(normalized_human_pbmc_pseudobulk_long$observation,sep = "-",field = 2)
normalized_human_pbmc_pseudobulk_long$type<-getField(normalized_human_pbmc_pseudobulk_long$observation,sep = "-",field = 3)
interest_gene<-c("PER2","BMAL1")
plotdata2<-normalized_human_pbmc_pseudobulk_long[normalized_human_pbmc_pseudobulk_long$feature %in% interest_gene,]

ggplot(plotdata2[plotdata2$type=="Monocyte",])+geom_point(aes(x=CT,y=raw_count,color=feature))+
  geom_line(aes(x=CT,y=raw_count,color=feature,group=feature))+
  scale_x_discrete(limits=CT_TIME_ORDER_FINE)+facet_wrap(~individual,scales = "free_y")+
  scale_color_d3("category20")+ggtitle(paste0(paste(interest_gene,collapse = "_"),"-Monocyte"))+
  theme(axis.text.x=element_text(angle=60,hjust=1))

# View other metadata
meta.plotdata<-test@meta.data %>% group_by(CT,individual,type) %>% summarise(n.cell=n())
ggplot(meta.plotdata[meta.plotdata$type=="Monocyte",])+
  geom_bar(aes(x=CT,y=n.cell),stat="identity")+
  facet_wrap(~individual,scales = "free")+theme(axis.text.x=element_text(angle=60,hjust=1))+
  ggtitle("number of Monocyte")

meta.plotdata<-test@meta.data %>% group_by(CT,individual,type) %>% summarise(median.n.feature=median(nFeature_RNA))
ggplot(meta.plotdata[meta.plotdata$type=="Monocyte",])+
  geom_bar(aes(x=CT,y=median.n.feature),stat="identity")+
  facet_wrap(~individual,scales = "free")+theme(axis.text.x=element_text(angle=60,hjust=1))+
  ggtitle("median nFeature Monocyte")

meta.plotdata<-test@meta.data %>% group_by(CT,individual,type) %>% summarise(median.n.count=median(nCount_RNA))
ggplot(meta.plotdata[meta.plotdata$type=="Monocyte",])+
  geom_bar(aes(x=CT,y=median.n.count),stat="identity")+
  facet_wrap(~individual,scales = "free")+theme(axis.text.x=element_text(angle=60,hjust=1))+
  ggtitle("median nCount Monocyte")
# investigaet if age affect circadian rhythm
## 1. whether age affect intensity of circadian rhythm
JTK_result_all<-readRDS('~/analysis/circadian/R/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
JTK_result_all<-JTK_result_all %>% group_by(CycID,individual) %>% summarise(median_fp2t=median(fold_peak2trough))
JTK_result_all$age<-AGES[JTK_result_all$individual]
JTK_result_all$age_group<-ifelse(JTK_result_all$age<40,"young","aged")
JTK_result_all<-JTK_result_all %>% group_by(CycID,age_group) %>% mutate(n.obervation=n())
JTK_result_filtered<-JTK_result_all %>% dplyr::filter(n.obervation>=3,median_fp2t<3)

# Create the boxplot with statistical annotations
ggplot(JTK_result_filtered, aes(x = age_group, y = median_fp2t)) +
  geom_boxplot(aes(fill = age_group), width = 0.6, alpha = 0.8, outliers = F) +
  stat_compare_means(
    method = "wilcox.test", # or "wilcox.test" for non-parametric
    comparisons = list(c("young", "aged")), # adjust based on your groups
    label = "p.signif", # shows asterisks
    tip.length = 0.01,
    vjust = 0.5) + scale_fill_manual(values = c("#1F77B4", "#FF7F0E"))+
  labs(x = "Age Group",y = "Fold Change (Peak to Trough)",title = "") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", size = 12),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )


JTK_result_filtered<-JTK_result_all %>% dplyr::filter(n.obervation>=3)
plt.list<-list()
for(gene in c("BMAL1","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2","DBP")){
  plt.list[[gene]]=ggplot(JTK_result_filtered[JTK_result_filtered$CycID==gene,], aes(x = age_group, y = median_fp2t)) +
    geom_boxplot(aes(fill = age_group), width = 0.6, alpha = 0.8, outliers = F) +
    stat_compare_means(
      method = "wilcox.test", # or "wilcox.test" for non-parametric
      comparisons = list(c("young", "aged")), # adjust based on your groups
      label = "p.signif", # shows asterisks
      tip.length = 0.01,
      vjust = 0.1) + scale_fill_manual(values = c("#1F77B4", "#FF7F0E"))+scale_y_continuous(expand = c(0,0.5))+
    labs(x = "",y = "",title = gene) +
    theme_classic() +
    theme(
      text = element_text(family = "Arial", size = 11),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "none"
    )
}
ggarrange(plotlist = plt.list)

JTK_result_filtered<-JTK_result_all %>% dplyr::filter(n.obervation>=3)
plt.list<-list()
for(gene in c("DDIT4","FKBP5","CXCR4","TSC22D3","CPT1A","LINC01619","PTOV1","STEAP1B","SLC25A20")){
  plt.list[[gene]]=ggplot(JTK_result_filtered[JTK_result_filtered$CycID==gene,], aes(x = age_group, y = median_fp2t)) +
    geom_boxplot(aes(fill = age_group), width = 0.6, alpha = 0.8, outliers = F) +
    stat_compare_means(
      method = "wilcox.test", # or "wilcox.test" for non-parametric
      comparisons = list(c("young", "aged")), # adjust based on your groups
      label = "p.signif", # shows asterisks
      tip.length = 0.01,
      vjust = 0.1) + scale_fill_manual(values = c("#1F77B4", "#FF7F0E"))+scale_y_continuous(expand = c(0,1))+
    labs(x = "",y = "",title = gene) +
    theme_classic() +
    theme(
      text = element_text(family = "Arial", size = 11),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "none"
    )
}
ggarrange(plotlist = plt.list)

##2. whether younsters experience more phase shift
JTK_result_all<-readRDS('~/analysis/circadian/R/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
JTK_result_all<-JTK_result_all %>% group_by(CycID,individual) %>% summarise(median.phase=median(LAG))
JTK_result_all$age<-AGES[JTK_result_all$individual]
JTK_result_all$age_group<-ifelse(JTK_result_all$age<40,"young","aged")
JTK_result_all<-JTK_result_all %>% group_by(CycID,age_group) %>% mutate(n.obervation=n())
ggplot(JTK_result_all[JTK_result_all$CycID=="CXCR4",])+
  geom_boxplot(aes(x=age_group,y=median.phase))+geom_jitter(aes(x=age_group,y=median.phase))

### RNA velosity
## generate JTK rds file from spliced and unspliced file
JTK.spliced<-readJTK('~/analysis/circadian/R/velosity.spliced.meta2d.by.type/')
JTK.spliced<-JTK.spliced %>% dplyr::filter(PER>=20,PER<=28,ADJ.P<0.05)
saveRDS(JTK.spliced,'~/analysis/circadian/R/JTK.spliced.filtered.16individual.rds')

JTK.unspliced<-readJTK('~/analysis/circadian/R/velosity.unspliced.meta2d.by.type/')
JTK.unspliced<-JTK.unspliced %>% dplyr::filter(PER>=20,PER<=28,ADJ.P<0.05)
saveRDS(JTK.unspliced,'~/analysis/circadian/R/JTK.unspliced.filtered.16individual.rds')

# merge 7 pooled sample with 12 individual sample
srt.7mix<-readRDS('~/analysis/circadian/R/7mixed.unspliced.metacell.srt.rds')
srt.12individual<-readRDS('~/analysis/circadian/R/12individual.unspliced.metacell.srt.rds')
srt.total<-merge(srt.7mix,srt.12individual)
srt.total<-JoinLayers(srt.total)
srt.total
saveRDS(srt.total,'~/analysis/circadian/R/16individual.unspliced.srt.rds')
JTK<-readRDS('~/analysis/circadian/R/JTK.spliced.filtered.16individual.addmedianexp.rds')

# read ready spliced/unspliced JTK files
spliced<-readRDS('~/analysis/circadian/R/JTK.spliced.filtered.16individual.addmedianexp.addp2t.rds')
unspliced<-readRDS('~/analysis/circadian/R/JTK.unspliced.filtered.16individual.addmedianexp.addp2t.rds')

# read spliced metacell file
spliced.srt<-readRDS('~/analysis/core_data/16individual.spliced.srt.rds')

spliced$CycID_Long<-paste0(spliced$CycID,"_",spliced$individual,"_",spliced$celltype)
unspliced$CycID_Long<-paste0(unspliced$CycID,"_",unspliced$individual,"_",unspliced$celltype)
shared<-intersect(spliced$CycID_Long,unspliced$CycID_Long)
spliced.filtered<-spliced[spliced$CycID_Long%in%shared,]
unspliced.filtered<-unspliced[unspliced$CycID_Long%in%shared,]

feature.by.lag<-(spliced.filtered[spliced.filtered$celltype=="CD14-Mono"&
                      spliced.filtered$individual=="KD"&spliced.filtered$fold_peak2trough>1.5,] %>% arrange(LAG))$CycID
plotGeneExpressionHeatmapByCT(spliced.srt,"CD14-Mono",individual = "KD",gene.list=feature.by.lag)




# read ready regular JTK file
regular<-readRDS('~/analysis/circadian/R/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')


### calculate correlation of FKBP5 and cortisol
test<-readRDS('~/analysis/core_data/16individual.srt.annotated.rds')
test@meta.data$type<-CELL_TYPES2[test@meta.data$predicted.celltype.l2]
test@meta.data$type<-test@meta.data$predicted.celltype.l2
plotPseudobulk(test,"NK_CD56bright",c("FKBP5"),bool.relative = T,bool.debug = T,
               char.cortisol.file = "~/analysis/core_data/ELISA_cortisol.rds")
data<-runtime.plotPseudobulk.mat.final[runtime.plotPseudobulk.mat.final$individual%in%c("JJC","KD","LYH","ZYR"),]
data<-data[c("features","CT","individual","values")]
data.fkbp5<-data[data$features=="SDHD",c("CT","individual","values")]
data.cortisol<-data[data$features=="cortisol",c("CT","individual","values")]
colnames(data.fkbp5)[3]<-"SDHD"
colnames(data.cortisol)[3]<-"cortisol"
data.merged<-left_join(data.fkbp5,data.cortisol,by=c("CT","individual"))
cotest<-cor.test(data.merged$SDHD,data.merged$cortisol,method = "spearman")

ggplot(data.merged)+geom_point(aes(x=SDHD,y=cortisol,color=individual))+
  geom_text(x=-1.2,y=1.5,label=paste0("R=",round(cotest$estimate,2)," p=",round(cotest$p.value,2)))+
  theme_half_open()

#hormone_test<-readRDS("~/analysis/core_data/Correlation.features2hormone.rds")
hormone_test_byindividual<-readRDS("~/analysis/R/Correlation.features2hormone.byindividual.rds")
#hormone_test<-hormone_test[!grepl("ENSG",hormone_test$features)&hormone_test$cor2cortisol>0.3,]
hormone_test_byindividual<-hormone_test_byindividual[!grepl("ENSG",hormone_test_byindividual$features)&
                                                       hormone_test_byindividual$p.value.cortisol<0.05&hormone_test_byindividual$cor2cortisol>0.3,]
hormone_test_cortisol_summarise<-hormone_test_byindividual %>% group_by(celltype,individual) %>% summarise(n.correlated2cortisol=n())
JTK<-readRDS('~/analysis/core_data/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
JTK_summarise<-JTK[JTK$individual %in% INDIVIDUALS_BATCH2,] %>% group_by(celltype,individual) %>% summarise(n.oscillating=n())
hormone_test_cortisol_summarise<-left_join(hormone_test_cortisol_summarise,JTK_summarise,by=c("celltype","individual"))
hormone_test_cortisol_summarise$ratio<-hormone_test_cortisol_summarise$n.correlated2cortisol/hormone_test_cortisol_summarise$n.oscillating
hormone_test_cortisol_summarise$celltype_rough<-CELL_TYPES2[hormone_test_cortisol_summarise$celltype]
ggplot(hormone_test_cortisol_summarise[!is.na(hormone_test_cortisol_summarise$celltype),])+
  geom_boxplot(aes(x=celltype,y=ratio,fill=celltype_rough))+theme_half_open()+
  theme(axis.text.x = element_text(angle=60,hjust=1))

### plot metacell to compare with pseudobuk
test.seacell<-readRDS('~/analysis/core_data/seacell.16individual.TF.rds')
test.seacell<-NormalizeData(test.seacell,normalization.method = "RC",scale.factor = 1000000)
plotMetaCellByIndividual(test.seacell,cell.type = "CD14 Mono",feature = "CRY1",
    layer = "data",merge = F,rm.outliers = T,float.alpha = 0.1,bool.maskzero = F)
test<-readRDS('~/analysis/core_data/16individual.srt.annotated.rds')
test$type<-test$predicted.celltype.l2
test<-NormalizeData(test,normalization.method = "RC",scale.factor = 1000000)
plotPseudobulk(test,char.cell.type = "CD14 Mono",vector.features = c("PER2","CRY1"))

# compare expression between each subsets
data.raw<-readRDS('~/analysis/circadian/R/16individual.pseudobulk.bytype.CPM.rds')
data.raw<-readRDS('~/analysis/circadian/R/16individual.pseudobulk.bytype.CPHC.rds')
data<-data.raw[CIRCADIAN_GENES_MAIN,]
#data<-data.raw[c("RORA", "RORB", "RORC","NFIL3","CSNK1D","CSNK1E","GSK3B","FBXL3"),]
data<-rownames_to_column(data,"feature")
#data<-gather(data,key="observation",value="CPHC",-feature)
data<-gather(data,key="observation",value="CPM",-feature)
data$type<-getField(data$observation,"-",3)
data$type_top_level<-CELL_TYPES2[data$type]
data<-data[!data$type%in%c("Platelet","Eryth","HSPC"),]
#data$type<-orderCellTypes(data$type,reverse = T)
data<-as.data.frame(data)
#data<-data[!is.na(data$CPHC),]
data<-data[!is.na(data$CPM),]
#data <- data %>%
#  group_by(feature, type) %>%  # group by whatever factors you need
#  mutate(
#    Q1 = quantile(CPHC, 0.25, na.rm = TRUE),
#    Q3 = quantile(CPHC, 0.75, na.rm = TRUE),
#    IQR = Q3 - Q1
#  ) %>%
#  filter(CPHC > Q1 - 1.5 * IQR & CPHC < Q3 + 1.5 * IQR) %>%
#  ungroup()
data <- data %>%
  group_by(feature, type) %>%  # group by whatever factors you need
  mutate(
    Q1 = quantile(CPM, 0.25, na.rm = TRUE),
    Q3 = quantile(CPM, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1
  ) %>%
  filter(CPM > Q1 - 1.5 * IQR & CPM < Q3 + 1.5 * IQR) %>%
  ungroup()
#data<-data %>% group_by(feature,type) %>% mutate(median_CPHC=round(median(CPHC)),max_CPHC=round(max(CPHC)))
#data<-data %>% group_by(feature) %>% arrange(.,median_CPHC,.by_group = T)
data<-data %>% group_by(feature,type) %>% mutate(median_CPM=round(median(CPM)),max_CPM=round(max(CPM)))
data<-data %>% group_by(feature) %>% arrange(.,median_CPM,.by_group = T)
plotlist<-list()
for (feature in CIRCADIAN_GENES_MAIN) {
  plotlist[[feature]]=ggplot(data[data$feature==feature,])+
    #geom_boxplot(aes(x=fct_inorder(type),y=CPHC,fill=type_top_level),outliers = F)+
    #geom_text(aes(x=fct_inorder(type),y=max_CPHC,label=median_CPHC),size=3,hjust=-0.5)+
    geom_boxplot(aes(x=fct_inorder(type),y=CPM,fill=type_top_level),outliers = F)+
    geom_text(aes(x=fct_inorder(type),y=max_CPM,label=median_CPM),size=3,hjust=-0.5)+
    scale_y_continuous(limits = c(0,180))+
    facet_wrap(~feature,scale="free",nrow = 2,ncol=5)+
    coord_flip()+xlab("")+scale_fill_bmj()+theme_minimal_vgrid(font_size = 9)+NoLegend()+
    theme(strip.background = element_rect(fill = "white", color = "black"))
}
ggarrange(plotlist = plotlist,ncol = 5,nrow = 2)


# try to get cells that express target gene
cell.data<-LayerData(test,layer = "count",features=c("CXCR4","DDIT4")) %>% as.data.frame()
cell.data<-rownames_to_column(cell.data,"feature")
cell.data<-gather(cell.data,key="observation",value="UMI_count",-feature)
cell.data.HQ<-cell.data[cell.data$UMI_count>=3,]
cell.data.HQ<-cell.data.HQ %>% group_by(observation) %>% mutate(valid.gene=n())
cell.data.HQ<-cell.data.HQ[cell.data.HQ$valid.gene==2,]
cell.data.HQ<-cell.data.HQ %>% spread(.,key=feature,value=UMI_count)
cell.data.HQ$individual<-getField(cell.data.HQ$observation,"_",2)
cell.data.HQ$CT<-getField(cell.data.HQ$observation,"_",3)

meta<-test@meta.data[,c("predicted.celltype.l2","nCount_RNA")]
meta<-rownames_to_column(meta,"observation")
cell.data.HQ<-left_join(cell.data.HQ,meta)
ggplot(cell.data.HQ[cell.data.HQ$predicted.celltype.l2=="CD14 Mono",])+
  geom_boxplot(aes(x=CT,y=DDIT4/CXCR4),outliers = F)+facet_wrap(~individual,scales = "free_y")

# calculate count per hundreds cell (CPHC) and save it to rds
mat<-srt2bulkMatrix(test,split.by=c("individual","CT","type"),bool.use.CPHC = T)
saveRDS(mat,'~/analysis/circadian/R/16individual.pseudobulk.bytype.CPHC.rds')

# correlation between genes and hormones
result<-readRDS('~/analysis/circadian/R/Correlation.features2hormone.rds')
result<-result %>% group_by(features,celltype) %>% mutate(n.individual=n(),mean.cor2melatonin=mean(cor2melatonin),
                                                          mean.cor2cortisol=mean(cor2cortisol),
                                                          mean.logp.melatonin=mean(-log10(p.value.melatonin)),
                                                          mean.logp.cortisol=mean(-log10(p.value.cortisol)))
result<-dplyr::filter(result,n.individual>=3)
result %>% arrange(desc(mean.cor2melatonin))
plotPseudobulk(test,"CD4 TCM","CAMK4",bool.relative = T,char.melatonin.file = '~/analysis/circadian/R/ELISA_melatonin.rds')

ggplot(unique(result[c("features","mean.cor2melatonin","mean.logp.melatonin")]))+
  geom_point(aes(x=mean.cor2melatonin,y=mean.logp.melatonin),alpha=0.5)+
  geom_text_repel(aes(x=mean.cor2melatonin,y=mean.logp.melatonin,label=features))+
  theme_half_open()

ggplot(unique(result[c("features","mean.cor2cortisol","mean.logp.cortisol")]))+
  geom_point(aes(x=mean.cor2cortisol,y=mean.logp.cortisol),alpha=0.5)+
  geom_text_repel(aes(x=mean.cor2cortisol,y=mean.logp.cortisol,label=features))+
  theme_half_open()

# VDJ analysis, BCR
VDJ.B<-readRDS('analysis/core_data/VDJ.B.clone_types.filtered.rds')
this.individual<-INDIVIDUALS_BATCH2[4]
VDJ.B.clone.at.CT<-VDJ.B[c("individual","CT","clone_id")] %>% unique()
VDJ.B.clone.at.CT<-VDJ.B.clone.at.CT %>% group_by(clone_id,individual) %>% mutate(n.occur=n())
VDJ.B.clone.at.CT$n.occur %>% table() %>% as.data.frame()

VDJ.B.clone.at.CT<-VDJ.B.clone.at.CT[VDJ.B.clone.at.CT$n.occur==2&VDJ.B.clone.at.CT$individual==this.individual,]
# Count occurrences per clone per CT
count_data <- VDJ.B.clone.at.CT %>%
  count(clone_id, CT) %>%
  group_by(clone_id) %>%
  mutate(percent = n / sum(n) * 100)


seqential_occur<-data.frame("clone_id"=character(0),"is.sequential"=character(0))
for (id in unique(count_data$clone_id)) {
  this.data=count_data[count_data$clone_id==id,]
  ct.time=this.data$CT %>% gsub("CT","",.) %>% as.numeric()
  ct.time.min=min(ct.time)
  if((ct.time.min+4) %in% ct.time){
    seqential_occur=rbind(seqential_occur,data.frame("clone_id"=id,"is.sequential"="yes"))
  }else{
    seqential_occur=rbind(seqential_occur,data.frame("clone_id"=id,"is.sequential"="no"))
  }
}

binom.test(x = table(seqential_occur$is.sequential)["yes"], n = nrow(seqential_occur), p = 2/7, alternative = "greater")

plot1<-ggplot(VDJ.B.clone.at.CT, aes(x = clone_id, y = CT, fill = CT)) +
  geom_tile(color="black") + 
  labs(title = "Clone Distribution Across Circadian Times",
       x = paste0("Clone ID of individual ",INDIVIDUAL_MASKING[this.individual]), y = "Circadian time", fill = "Circadian Time") +
  coord_flip() + scale_y_discrete(expand = c(0,0))+
  scale_fill_manual(values=c("#FFFFFF","#D4D4D4","#AAAAAA","#555555","#AAAAAA","#D4D4D4","#FFFFFF"))+
  theme(panel.background = element_rect(fill="white",color="black",linetype = 1),axis.text.y=element_blank(),axis.ticks.y = element_blank())

plot2<-ggplot(seqential_occur)+geom_tile(aes(x="is.sequential",y=clone_id,fill=is.sequential),color="black")+
  scale_fill_manual(values = c("white","red"))+theme_nothing()+theme(plot.margin = margin(0,1,0,0,"cm"))+
  scale_y_discrete(expand = c(0,0))+scale_x_discrete(expand = c(0,0))

ggarrange(plot1,plot2,common.legend = T,align = "hv",widths = c(3,1),legend = "bottom")

# VDJ analysis, TCR
VDJ.T<-readRDS('~/analysis/core_data/VDJ.T.clone_types.filtered.rds')
this.individual<-INDIVIDUALS_BATCH1[12]
VDJ.T.clone.at.CT<-VDJ.T[c("individual","CT","clone_id")] %>% unique()
VDJ.T.clone.at.CT<-VDJ.T.clone.at.CT %>% group_by(clone_id,individual) %>% mutate(n.occur=n())
VDJ.T.clone.at.CT$n.occur %>% table() %>% as.data.frame()

VDJ.T.clone.at.CT<-VDJ.T.clone.at.CT[VDJ.T.clone.at.CT$n.occur==2&VDJ.T.clone.at.CT$individual==this.individual,]
# Count occurrences per clone per CT
count_data <- VDJ.T.clone.at.CT %>%
  count(clone_id, CT) %>%
  group_by(clone_id) %>%
  mutate(percent = n / sum(n) * 100)

seqential_occur<-data.frame("clone_id"=character(0),"is.sequential"=character(0))
for (id in unique(count_data$clone_id)) {
  this.data=count_data[count_data$clone_id==id,]
  ct.time=this.data$CT %>% gsub("CT","",.) %>% as.numeric()
  ct.time.min=min(ct.time)
  if((ct.time.min+4) %in% ct.time){
    seqential_occur=rbind(seqential_occur,data.frame("clone_id"=id,"is.sequential"="yes"))
  }else{
    seqential_occur=rbind(seqential_occur,data.frame("clone_id"=id,"is.sequential"="no"))
  }
}

binom.test(x = table(seqential_occur$is.sequential)["yes"], n = nrow(seqential_occur), p = 1/3, alternative = "greater")
order.clone_id<-(seqential_occur %>% arrange(is.sequential))$clone_id
plot1<-ggplot(VDJ.T.clone.at.CT, aes(x = clone_id, y = CT, fill = CT)) +
  geom_tile(color="black") + 
  labs(title = "Clone Distribution Across Circadian Times",
       x = paste0("Clone ID of individual ",INDIVIDUAL_MASKING[this.individual]), y = "Circadian time", fill = "Circadian Time") +
  scale_x_discrete(limits=order.clone_id)+scale_y_discrete(expand = c(0,0))+
  scale_fill_manual(values=c("#FFFFFF","#D4D4D4","#AAAAAA","#555555","#AAAAAA","#D4D4D4","#FFFFFF"))+
  theme(panel.background = element_rect(fill="white",color="black",linetype = 1),axis.text.x=element_blank(),axis.ticks.x = element_blank())

plot2<-ggplot(seqential_occur)+geom_tile(aes(x=clone_id,y="is.sequential",fill=is.sequential),color="black")+
  scale_fill_manual(values = c("white","red"))+theme_nothing()+theme(plot.margin = margin(0,0,0,0,"cm"))+
  scale_y_discrete(expand = c(0,0))+scale_x_discrete(expand = c(0,0),limits=order.clone_id)

ggarrange(plot1,plot2,common.legend = T,ncol = 1,align = "hv",heights = c(3,1),legend = "bottom")

plot.list<-list()
for(this.individual in c(INDIVIDUALS_BATCH1,INDIVIDUALS_BATCH2)){
  this.data=VDJ.T[VDJ.T$individual==this.individual,]
  this.data=this.data[c("CT","clone_id")] %>% unique()
  this.data=this.data %>% group_by(clone_id) %>% mutate(n.occur=n())
  plotdata=this.data$n.occur %>% table() %>% as.data.frame()
  colnames(plotdata)[1]<-"occurrence"
  plotdata=plotdata %>%
    arrange(desc(Freq)) %>%
    mutate(
      percentage = Freq / sum(Freq) * 100,
      label_position = cumsum(Freq) - Freq / 2,  # Position for labels
      label = paste0(occurrence, " (", round(percentage, 1), "%)")
    )
  plot.list[[this.individual]]=ggplot(plotdata, aes(x = "", y = Freq, fill = label)) +
    geom_bar(stat = "identity", width = 1, linewidth = 0.5,color="black") +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(subtitle=INDIVIDUAL_MASKING[this.individual]) +
    scale_fill_aaas(name="occurrence")+
    theme(plot.title = element_text(hjust=0.5),text=element_text(size=14,color="black"))
}
ggarrange(plotlist = plot.list,ncol=4,nrow=4)

# CYTOF data to confirm rhythmicity of CXCR4
cytof<-readRDS('~/analysis/core_data/CYTOF.300k.reannotation.srt.rds')
test<-readRDS('~/analysis/core_data/16individual.srt.annotated.rds')
test.metacell<-readRDS('~/analysis/core_data/seacell.16individual.TF.rds')
test$type<-test$predicted.celltype.l2
test<-NormalizeData(test,normalization.method = "RC",scale.factor = 1000000)
test.metacell<-NormalizeData(test.metacell,normalization.method = "RC",scale.factor = 1000000)

JTK_CYTOF<-readRDS('~/analysis/core_data/CYTOF.JTK.newresult.filtered.addmedianexp.addp2t.bytype.rds')
JTK_metacell<-readRDS('~/analysis/core_data/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')


this.celltype<-"TNK"

JTK_CYTOF[JTK_CYTOF$celltype==this.celltype&JTK_CYTOF$ADJ.P<0.000001,"CycID"] %>% table() %>% sort()

this.feature<-as.vector(CYTOF2RNA_FEATURES["CD45-plt"])
#this.feature<-"PER2"
JTK_metacell[JTK_metacell$CycID==this.feature&JTK_metacell$celltype==this.celltype,]
JTK_CYTOF[JTK_CYTOF$CycID=="CD184-plt"&JTK_CYTOF$celltype==this.celltype,]

plot1<-plotPseudobulk(test,this.celltype,this.feature,c("KD","JJC","LYH","ZYR"))+
  ggtitle(paste0("RNA: pseudobulk"))+NoLegend()+xlab("")

plot2<-plotMetaCellByIndividual(test.metacell,this.celltype,this.feature,
  individuals = c("TF_KD","TF_JJC","TF_LYH","TF_ZYR"),
  layer="count",float.alpha = 0.07)+xlab("")+
  ggtitle(paste0("RNA: metacell"))

plot3<-plotRandomPseudobulkByCTByIndividual(srt = cytof,cell.type = this.celltype,
  individuals = c("TF_KD","TF_JJC","TF_LYH","TF_ZYR"),
  features=c(RNA2CYTOF[this.feature]),normalize.data=F,proportion = 1,
  split.by="individual",plot="raw",float.alpha = 1)+NoLegend()+
  ggtitle(paste0("CYTOF: pseudobulk"))

plot4<-plotCytofbyIndividual(cytof,this.celltype,this.feature,proportion = 1,
  individuals = c("TF_KD","TF_JJC","TF_LYH","TF_ZYR"))+xlab("")+
  ggtitle(paste0("CYTOF: raw"))

(plot1|plot2)/(plot3|plot4)+patchwork::plot_annotation(paste0(this.celltype,": ",this.feature))

this.celltype<-"CD16 Mono"
this.feature<-"CCR4"
plot1<-plotPseudobulk(test,this.celltype,this.feature)+
  ggtitle(paste0("RNA: pseudobulk"))+NoLegend()+xlab("")

plot2<-plotMetaCellByIndividual(test.metacell,this.celltype,this.feature,
                                layer="count",float.alpha = 0.07)+xlab("")+
  ggtitle(paste0("RNA: metacell"))
(plot1|plot2)+patchwork::plot_annotation(paste0(this.celltype,": ",this.feature))

cytof@meta.data$type<-CELL_TYPES2[cytof@meta.data$manual.new]
test@meta.data$type<-CELL_TYPES2[test@meta.data$predicted.celltype.l2]
test.metacell@meta.data$type<-CELL_TYPES2[test.metacell@meta.data$predicted.celltype.l2.main]

# merge 12 individual velosity data
mat<-getMultipleWholeDaySplicedData(INDIVIDUALS_BATCH1)
mat2=list()
for(sample_id in sort(names(BATCH_JJC))){
  sample_loom_path=paste0(sample_id,"/outs/per_sample_outs/",sample_id,"/count/velocyto/")
  sample_loom_file_name=list.files(sample_loom_path)[1]
  sample_loom_full_path=paste0(sample_loom_path,sample_loom_file_name)
  message("loom file: ")
  message(sample_loom_full_path)
  loomMatrix=read.loom.matrices(sample_loom_full_path)
  #loomMatrix=subsetLoomMatrix(loomMatrix,cell.preffix=i,cells=cells)
  colnames(loomMatrix$spliced)=paste(sample_id,colnames(loomMatrix$spliced),sep='_')
  colnames(loomMatrix$unspliced)=paste(sample_id,colnames(loomMatrix$unspliced),sep='_')
  colnames(loomMatrix$ambiguous)=paste(sample_id,colnames(loomMatrix$ambiguous),sep='_')
  mat2[[sample_id]]=loomMatrix
}
saveRDS(mat,"~/analysis/circadian/R/scVelocyto.12individual.rds")
saveRDS(mat2,"~/analysis/circadian/R/scVelocyto.4PooledIndividual.rds")

mat<-readRDS("~/analysis/R/scVelocyto.12individual.rds")
srts<-loom2SeuratObject("~/analysis/R/scVelocyto.12individual.rds",'~/analysis/R/cell.annotations.manual1_2_NI.predicted2.modforSeacells.tsv')

# compare phase of circadian genes in different cell type
KD_CD14_Monocyte<-read.delim('FineTypeLevelPseudobulk/meta2d_KD_CD14_Mono.txt')
KD_CD16_Monocyte<-read.delim('FineTypeLevelPseudobulk/meta2d_KD_B_intermediate.txt')

shared.CG<-intersect(KD_CD14_Monocyte[KD_CD14_Monocyte$meta2d_pvalue<0.05,"CycID"],
          KD_CD16_Monocyte[KD_CD16_Monocyte$meta2d_pvalue<0.05,"CycID"])
plotdata<-left_join(KD_CD14_Monocyte[KD_CD14_Monocyte$CycID %in% shared.CG,c("CycID","meta2d_phase")],
                    KD_CD16_Monocyte[KD_CD16_Monocyte$CycID %in% shared.CG,c("CycID","meta2d_phase")],by="CycID")
colnames(plotdata)[2:3]<-c("CD14 Monocyte","B memory")
plotdata2<-gather(plotdata,key="type",value="phase",-CycID)
plot1<-ggplot(plotdata)+geom_point(aes(x= get("CD14 Monocyte"),y=get("B memory")))+
  xlab("CD14 Monocyte")+ylab("B memory")+
  geom_abline(slope = 1,intercept = 0)+theme_minimal(base_size = 13)
plot2<-ggplot(plotdata2)+geom_density(aes(x=phase,fill=type),alpha=0.3)+
  theme_minimal(base_size = 13)+theme(legend.position=c(0.8,0.9))+
  guides(fill=guide_legend(title=NULL))
ggarrange(plot1,plot2,nrow=2,align="hv")

# impose metacell on pseudobulk
test.metacell<-readRDS('~/analysis/core_data/seacell.16individual.TF.rds')
test.table<-readRDS('~/analysis/core_data/16individual.pseudobulk.bytype.CPM.rds')
test.metacell<-NormalizeData(test.metacell,normalization.method = "RC",scale.factor = 1000000)
imposeMetaCellOnPseudoBulk(srt.metacell = test.metacell,char.cell.type="CD14 Mono",char.feature="GAPDH",
                           df.pseudobulk = test.table)

#
plotPseudobulk(S4.srt = NULL,char.cell.type = "CD14 Mono",vector.features = "LGALS3BP",char.data.from.mat = test.table)

# View pseudo bulk level 
test.meta2d<-read.delim('~/analysis/R/FineTypeLevelPseudobulk/meta2d_KD_CD14_Mono.txt')
test.meta2d<-test.meta2d[!is.na(test.meta2d$meta2d_pvalue),]
test.meta2d[test.meta2d$CycID %in% CIRCADIAN_GENES_MAIN,]
test.meta2d.filtered<-test.meta2d[test.meta2d$meta2d_pvalue<0.1,]
meta2d.all<-NULL
for (file in list.files('~/analysis/R/FineTypeLevelPseudobulk/',pattern = "^meta2d*")) {
  message(file)
  this.file=file.path('~/analysis/R/FineTypeLevelPseudobulk',file)
  this.data=read.delim(this.file)
  this.data=this.data[!is.na(this.data$meta2d_pvalue),]
  if(nrow(this.data)==0){
    next
  }
  this.data$individual=getField(file,"_",2)
  this.data$celltype=getField(file,"_",-c(1,2)) %>% gsub(".txt","",.,fixed = T)
  if(is.null(meta2d.all)){
    meta2d.all=this.data
  }else{
    meta2d.all=rbind(this.data,meta2d.all)
  }
}
saveRDS(meta2d.all,'~/analysis/core_data/meta2d.pseudobulkByCelltypeByIndividual.rds')
meta2d.filtered<-meta2d.all[meta2d.all$meta2d_pvalue<0.1&meta2d.all$meta2d_rAMP>=0.25,]
meta2d.filtered<-meta2d.filtered %>% group_by(CycID,celltype) %>% mutate(n.subset.occurence=n())
meta2d.filtered<-meta2d.filtered[meta2d.filtered$n.subset.occurence>2,]
meta2d.filtered$CycID %>% table() %>% sort() %>% tail()

plotPseudobulk(S4.srt = NULL,char.cell.type = "B intermediate",vector.features = "IGLV3-21",char.data.from.mat = test.table)

# calculate heavy chain and light chain expression of TCR and BCR
test.table<-test.table[,!is.nan(colSums(test.table))]
i=1
all_chain=NULL
names_chain<-c("BCR_light_lambda","BCR_light_kappa","BCR_heavy","TCR_alpha","TCR_beta","TCR_gamma","TCR_delta")
for(this.pattern in c("IGL[VDJC]","IGK[VDJC]","IGH[AGDVJ]","TRA[VDJC]","TRB[VDJC]","TRG[VDJC]","TRD[VDJC]")){
  this.data=test.table[grepl(this.pattern,rownames(test.table)),] %>% colSums() %>% as.data.frame()
  colnames(this.data)[1]=names_chain[i]
  print(head(this.data))
  if(is.null(all_chain)){
    all_chain=this.data
  }else{
    all_chain=mergeByRownames(all_chain,this.data)
    print(head(all_chain))
  }
  i=i+1
}
all_chain<-as.data.frame(t(all_chain))
plotPseudobulk(S4.srt = NULL,char.cell.type = "B intermediate",
               vector.features = "BCR_light_kappa",char.data.from.mat = all_chain)

# merge cells to obtain raw count of metacells
metadata<-readRDS('~/analysis/core_data/seacell.0.04.16individual.metadata.rds')
srt<-readRDS('~/analysis/core_data/16individual.srt.annotated.rds')
srt<-NormalizeData(srt,normalization.method = "RC",scale.factor = 1000000)
mat<-LayerData(srt,layer = "data")
available.sc.cells<-colnames(mat)
rm(srt)
gc()
i=1
metacell.mat=NULL
for(this.metacell in metadata$SEACell %>% unique()){
  if(i<=17300){
    i=i+1
    next
  }
  this.sc.cells=metadata[metadata$SEACell==this.metacell,"index"]
  this.sc.cells=intersect(this.sc.cells,available.sc.cells)
  #message(this.metacell)
  this.mat=NULL
  if(length(this.sc.cells)==0){
    next
  }else if(length(this.sc.cells)==1){
    this.mat=as.data.frame(mat[,this.sc.cells])
  }else{
    this.mat=as.data.frame(rowSums(mat[,this.sc.cells]))
  }
  colnames(this.mat)=this.metacell
  if(is.null(metacell.mat)){
    metacell.mat=this.mat
  }else{
    metacell.mat=mergeByRownames(metacell.mat,this.mat)
  }
  i=i+1
  if(i%%100==0){
    message(paste0(i,"/",length(metadata$SEACell %>% unique())))
  }
}
saveRDS(metacell.mat,'~/analysis/core_data/seacell.0.04.16individual.matrix.part2.rds')

mat.part1<-readRDS('~/analysis/core_data/seacell.0.04.16individual.matrix.part1.rds')
mat.part2<-readRDS('~/analysis/core_data/seacell.0.04.16individual.matrix.part2.rds')
overlapped<-intersect(colnames(mat.part1),colnames(mat.part2))
mat.part2<-mat.part2[,!(colnames(mat.part2) %in% overlapped)]
all.equal(rownames(mat.part2),rownames(mat.part1))
mat.whole<-mergeByRownames(mat.part1,mat.part2)
rm(mat.part1,mat.part2)
srt.metacell<-CreateSeuratObject(mat.whole)
srt.metacell@meta.data$individual<-getField(rownames(srt.metacell@meta.data),"_",1)
srt.metacell@meta.data$CT<-getField(rownames(srt.metacell@meta.data),"_",2)
srt.metacell@meta.data$type<-getField(rownames(srt.metacell@meta.data),"_",3)
srt.metacell@meta.data$type<-gsub("-"," ",srt.metacell@meta.data$type)
srt.metacell<-RenameCells(srt.metacell,add.cell.id = "TF")
saveRDS(srt.metacell,'~/analysis/core_data/seacell.0.04.16individual.RawCount.srt.rds')
srt.metacell<-readRDS('~/analysis/core_data/seacell.0.04.16individual.RawCount.srt.rds')
srt.metacell<-NormalizeData(srt.metacell,normalization.method = "RC",scale.factor = 1000000)

test.table<-readRDS('~/analysis/core_data/16individual.pseudobulk.bytype.CPM.rds')
imposeMetaCellOnPseudoBulk(srt.metacell = srt.metacell,char.cell.type="CD14 Mono",char.feature="PER2",
                           df.pseudobulk = test.table)

meta2d.result<-readRDS('~/analysis/core_data/meta2d.pseudobulkByCelltypeByIndividual.rds')
meta2d.result[meta2d.result$meta2d_pvalue<=0.01&meta2d.result$CycID=="PER2"&meta2d.result$celltype=="CD14_Mono",]
plotPseudobulk(char.cell.type = "CD14 Mono",vector.features="PER2",
               char.data.from.mat=test.table,bool.return.plotdata = T,df.meta2d=meta2d.result,
               bool.merge = T)

velocyto<-readRDS('~/analysis/R/scVelocyto.12individual.rds')
# spliced rna data to srt
spliced.srt.list<-list()
i=1
for (individual in names(velocyto)) {
  for(this.CT in names(velocyto[[individual]])){
    mat=velocyto[[individual]][[this.CT]]$unspliced
    mat=mat[!duplicated(rownames(mat)), ]
    this.srt=CreateSeuratObject(mat)
    spliced.srt.list[[i]]=this.srt
    i=i+1
  }
}
spliced.srt<-merge(spliced.srt.list[[1]],spliced.srt.list[2:length(spliced.srt.list)])
spliced.srt[["RNA"]]<-JoinLayers(spliced.srt[["RNA"]])
spliced.srt<-RenameCells(spliced.srt,add.cell.id = "TF")
saveRDS(spliced.srt,"~/analysis/R/scVelocyto.12individual.unspliced.srt.rds")
head(Cells(spliced.srt))
rm(spliced.srt.list,spliced.srt)
rm(velocyto)
gc(reset = T)
velocyto<-readRDS('~/analysis/R/scVelocyto.4PooledIndividual.rds')
test<-readRDS('~/analysis/core_data/4multiplexedindividual.srt.annotated.rds')
test$raw_id<-paste0(test$sample,"_",test$orig.ident)
mappings<-as.data.frame(test$raw_id)
rm(test)
gc()
colnames(mappings)[1]="cell_id_old"
mappings=rownames_to_column(mappings,"cell_id_new")
mapping_value=mappings$cell_id_new
mapping_key=mappings$cell_id_old
names(mapping_value)=mapping_key
spliced.srt.list<-list()
i=1
for(sample in names(velocyto)){
  mat=as.data.frame(velocyto[[sample]]$unspliced)
  colnames(mat)=colnames(mat) %>% gsub("_sample_alignments_.*:","_",.) %>% gsub("x$","-1",.)
  print(head(colnames(mat)))
  #print(mat[1:3,1:3])
  mat=mat[colnames(mat) %in% mappings$cell_id_old]
  
  colnames(mat)=as.vector(mapping_value[colnames(mat)])
  mat=mat[!duplicated(rownames(mat)), ]
  this.srt=CreateSeuratObject(mat)
  spliced.srt.list[[i]]=this.srt
  i=i+1
  gc()
}
spliced.srt<-merge(spliced.srt.list[[1]],spliced.srt.list[2:length(spliced.srt.list)])
spliced.srt[["RNA"]]<-JoinLayers(spliced.srt[["RNA"]])
spliced.srt2<-readRDS('~/analysis/R/scVelocyto.12individual.unspliced.srt.rds')
spliced.srt<-merge(spliced.srt,spliced.srt2)
spliced.srt[["RNA"]]<-JoinLayers(spliced.srt[["RNA"]])
saveRDS(spliced.srt,'~/analysis/core_data/scVelo.16individual.unspliced.srt.rds')

test<-readRDS('~/analysis/core_data/16individual.srt.annotated.rds')
valid.cells<-Cells(test)
meta.data<-test@meta.data[,c("individual","CT","manual.level1","manual.level2","predicted.celltype.l2","manual_NI")]
rm(test)
gc()
srt<-readRDS('~/analysis/core_data/scVelo.16individual.unspliced.srt.rds')
srt<-subset(srt,cells=valid.cells)
srt<-AddMetaData(srt,metadata = meta.data)
srt$type<-srt$predicted.celltype.l2
srt<-NormalizeData(srt,normalization.method = "RC",scale.factor = 1000000)
plotPseudobulk(srt,"CD14 Mono",c("PER2"))
saveRDS(srt,'~/analysis/core_data/scVelo.16individual.unspliced.srt.rds')

srt<-readRDS('~/analysis/core_data/scVelo.16individual.spliced.srt.rds')
mat<-srt2bulkMatrix(srt,c("individual","CT","type"),layer = "data")
saveRDS(mat,'~/analysis/core_data/16individual.spliced.pseudobulk.bytype.CPM.rds')

# gene enrichment analysis for genes shared by most cell type/each celltype
meta2d.result0<-readRDS('~/analysis/core_data/meta2d.pseudobulkByCelltypeByIndividual.rds')
meta2d.result<-meta2d.result0[meta2d.result0$meta2d_pvalue<=0.05&!grepl("ENSG|LINC",meta2d.result0$CycID)&meta2d.result0$meta2d_rAMP>0.25,]
meta2d.result<-meta2d.result %>% group_by(CycID,celltype) %>% mutate(n.individual=n())
meta2d.result<-meta2d.result[meta2d.result$n.individual>=2,]
meta2d.result<-meta2d.result[c("CycID",'celltype')] %>% unique()
genes.shared<-meta2d.result$CycID %>% table() %>% sort() %>% tail(100) %>% names()
result<-enrichGObyHGNC(genes.shared)
plotdata<-result@result[result@result$pvalue<0.01,] %>% arrange(desc(pvalue))
plotdata<-plotdata[!duplicated(plotdata$geneID),]
ggplot(plotdata)+geom_point(aes(x=fct_inorder(Description),color=-log10(pvalue),
                             y=-log10(pvalue),size=Count))+coord_flip()+
  xlab("")+ylim(c(1.95,3))+theme_light(base_size = 13)

# compare phase of different celltype
meta2d.result0<-readRDS('~/analysis/core_data/meta2d.pseudobulkByCelltypeByIndividual.rds')
meta2d.result<-meta2d.result0[meta2d.result0$meta2d_pvalue<=0.1&!grepl("ENSG|LINC",meta2d.result0$CycID)&meta2d.result0$meta2d_rAMP>0.25,]
meta2d.result<-meta2d.result %>% group_by(CycID,celltype) %>% mutate(n.individual=n())
meta2d.result<-meta2d.result[meta2d.result$n.individual>=2,]
meta2d.result.CD14<-meta2d.result[meta2d.result$celltype %in% c("CD14_Mono"),
                                      c("CycID","individual","meta2d_phase")]
colnames(meta2d.result.CD14)[3]="CD14_Mono"
meta2d.result.CD16<-meta2d.result[meta2d.result$celltype %in% c("CD16_Mono"),
                                  c("CycID","individual","meta2d_phase")]
colnames(meta2d.result.CD16)[3]="CD16_Mono"
meta2d.result.CD14CD16<-left_join(meta2d.result.CD14,meta2d.result.CD16,by=c("CycID","individual"))
meta2d.result.CD14CD16<-meta2d.result.CD14CD16[!is.na(meta2d.result.CD14CD16$CD16_Mono),]
ggplot(meta2d.result.CD14CD16)+geom_point(aes(x=CD14_Mono,y=CD16_Mono))+geom_abline(slope = 1,intercept = 0)+
  geom_abline(slope = 1,intercept = 1,linetype=2)+geom_abline(slope = 1,intercept = -1,linetype=2)+
  geom_abline(slope = 1,intercept = 4,linetype=2)+geom_abline(slope = 1,intercept = -4,linetype=2)+
  geom_abline(slope = 1,intercept = 12,linetype=2)+geom_abline(slope = 1,intercept = -12,linetype=2)+
  geom_point(data=meta2d.result.CD14CD16[meta2d.result.CD14CD16$CycID %in% CIRCADIAN_GENES_MAIN,],
                  aes(x=CD14_Mono,y=CD16_Mono),shape=21,color="red",stroke=1)+
  geom_text_repel(data=meta2d.result.CD14CD16[meta2d.result.CD14CD16$CycID %in% CIRCADIAN_GENES_MAIN,],
            aes(x=CD14_Mono,y=CD16_Mono,label=CycID),hjust=-0.5)+theme_classic2()
comparePhaseBetweenCellTypes(meta2d.result0,"CD14_Mono","CD16_Mono",float.pvalue = 0.1,float.relative2AMP = 0.3)

plotPeakAmongIndividuals(meta2d.result,mode = "meta2d",gene = "CXCR4")

# using meta2d result to test tauFisher
meta2d.result0<-readRDS('~/analysis/core_data/meta2d.pseudobulkByCelltypeByIndividual.rds')
meta2d.result<-meta2d.result0[meta2d.result0$meta2d_pvalue<=0.1&!grepl("ENSG|LINC",meta2d.result0$CycID)&meta2d.result0$meta2d_rAMP>0.25,]
bulk_by_type<-readRDS("~/analysis/core_data/16individual.pseudobulk.bytype.CPM.rds")
top.genes<-CIRCADIAN_GENES_MAIN
predicted.data<-predictIndividualwithTauFisher(data=bulk_by_type,individual = c("KD","JJC","HZD","ZYX","WLG","XSP"),celltype = "CD14_Mono",
                               predictor.genes = top.genes,bool.include.trained = T,bool.return.pred = T)

# calculate dispersion of inter-individual phase shift between genes within individual
for (used.celltype in c("CD4_TCM","CD8_TEM","NK","B_memory","B_naive")) {
  plotdata=NULL
  for (this.individual in unique(meta2d.result$individual)) {
    message(this.individual)
    this.meta2d=meta2d.result[meta2d.result$individual==this.individual&meta2d.result$celltype==used.celltype,c("CycID","meta2d_phase")]
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

plotdata<-readRDS('~/analysis/R/lag.diffs/lag.diff.Treg.rds')
plotdata.filtered<-plotdata[plotdata$n.pair>6,]
plotdata.filtered<-plotdata.filtered %>% group_by(CycID,CycID2) %>% mutate(sd.lag.diff=sd(lag_diff))
plotdata.filtered<-plotdata.filtered[c("CycID","CycID2","sd.lag.diff")] %>% unique()
plotdata.mat<-spread(plotdata.filtered,key="CycID2",value="sd.lag.diff")
plotdata.mat<-as.data.frame(plotdata.mat)
plotdata.mat<-column_to_rownames(plotdata.mat,"CycID")

df<-plotdata.filtered
# First, create the distance matrix as before
all_cycids <- unique(c(df$CycID, df$CycID2))
dist_matrix <- matrix(NA, nrow = length(all_cycids), ncol = length(all_cycids),
                      dimnames = list(all_cycids, all_cycids))

for(i in 1:nrow(df)) {
  dist_matrix[df$CycID[i], df$CycID2[i]] <- df$sd.lag.diff[i]
  dist_matrix[df$CycID2[i], df$CycID[i]] <- df$sd.lag.diff[i]
}
diag(dist_matrix) <- 0

# Impute missing values with the mean distance
mean_distance <- mean(dist_matrix, na.rm = TRUE)
dist_matrix_imputed <- dist_matrix
dist_matrix_imputed[is.na(dist_matrix_imputed)] <- mean_distance

# Or use a more sophisticated imputation: maximum distance + small buffer
max_distance <- max(dist_matrix, na.rm = TRUE)
dist_matrix_imputed <- dist_matrix
dist_matrix_imputed[is.na(dist_matrix_imputed)] <- max_distance * 1.1

# Now create the heatmap with imputed matrix
col_fun <- colorRamp2(c(0, max(df$sd.lag.diff, na.rm = TRUE)), 
                      c("#003366", "#FFEE99"))

ht <- Heatmap(
  dist_matrix_imputed,
  name = "Distance",
  col = col_fun,
  na_col = "white",
  
  # Use the imputed matrix for clustering
  clustering_distance_rows = function(x) as.dist(dist_matrix_imputed),
  clustering_distance_columns = function(x) as.dist(dist_matrix_imputed),
  
  # CNS-level styling
  rect_gp = gpar(col = "white", lwd = 0.1),
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

# For PDF (vector format, best for publications)
svg("distance_heatmap.svg", 
    width = 25, 
    height = 25, 
    pointsize = 8) # CNS standard: 8-10 pt font

draw(ht) # Your Heatmap object
dev.off() # Essential: closes the graphics device

# Define your genes of interest
genes_of_interest<-c("PER1", "NR1D2", "BMAL1", "DDIT4", "FKBP5","PER2","PER3","CRY1","CRY2","NR1D1",'CLOCK') # Replace with your actual IDs

# Create the base heatmap without row/column labels
ht <- Heatmap(
  dist_matrix_imputed,
  name = "Phase shift\ndeviation",
  col = col_fun,
  na_col = "white",
  clustering_distance_rows = function(x) as.dist(dist_matrix_imputed),
  clustering_distance_columns = function(x) as.dist(dist_matrix_imputed),
  rect_gp = gpar(col = "white", lwd = 0.1),
  show_row_names = FALSE,    # Turn off all row labels
  show_column_names = FALSE, # Turn off all column labels
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

# Add text annotations for specific genes
# Note: This requires knowing the positions after clustering
# You might need to run the clustering first to get positions

# Alternative: Use decorate_annotation for more control
ht <- ht + rowAnnotation(
  foo = anno_mark(
    at = which(rownames(dist_matrix_imputed) %in% genes_of_interest),
    labels = rownames(dist_matrix_imputed)[rownames(dist_matrix_imputed) %in% genes_of_interest],
    labels_gp = gpar(fontsize = 8),
    link_gp = gpar(lwd = 0.5),
    padding = unit(1, "mm")
  )
)
draw(ht) # Your Heatmap object

# differential expression analysis of individuals with good rhythmicity and bad rhythmicity
mat<-readRDS("~/analysis/core_data/16individual.pseudobulk.bytype.rds")
mat<-mat[!is.na(colSums(mat))]
mat<-mat[(colSums(mat))!=0]
saveRDS(mat,"~/analysis/core_data/16individual.pseudobulk.bytype.rds")
mat<-readRDS("~/analysis/core_data/16individual.pseudobulk.bytype.rds")
for (celltype in unique(getField(colnames(mat),"-",3))) {
  mat.test<-mat[,grepl(celltype,colnames(mat))]
  mat.test<-rownames_to_column(mat.test,var = "feature")
  mat.test<-gather(mat.test,key="observation",value="value",-feature)
  mat.test$individual<-getField(mat.test$observation,"-",1)
  mat.test<-mat.test %>% group_by(feature,individual) %>% summarise(value=sum(value))
  mat.test<-spread(mat.test,key="individual",value="value")
  mat.test<-column_to_rownames(mat.test,var = "feature")
  mat.test<-mat.test[,colSums(mat.test)!=0]
  mat.test<-mat.test[rowMin(mat.test %>% as.matrix())>=3,]
  mat.test<-mat.test[!(grepl("^AC[0-9].*",rownames(mat.test))|grepl("^AL[0-9].*",
                                                                    rownames(mat.test))|grepl("^FP[0-9].*",rownames(mat.test))),]
  condition<-INDIVIDUAL2CIRCADIAN_GROUP2[colnames(mat.test)]
  colData <- data.frame(row.names=colnames(mat.test), condition)
  
  #dds<-DESeqDataSetFromMatrix(countData = mat.test, colData = colData, design = ~ condition+dornor)
  dds<-DESeqDataSetFromMatrix(countData = mat.test, colData = colData, design = ~ condition)
  dds.res <- DESeq(dds, fitType = 'mean')
  
  result<-results(dds.res) %>% as.data.frame()
  result$significant<-ifelse((-log10(result$padj))>=2&abs(result$log2FoldChange)>1,"yes","no")
  result$gene<-rownames(result)
  result<-result[!grepl("ENSG",result$gene),]
  result<-result[!is.na(result$padj),]
  gene.list=result[result$significant=="yes","gene"]
  if(length(gene.list)>10){
    message(length(gene.list))
    print(gene.list)
    enrich.res=enrichGObyHGNC(gene.list,backgroud.genes=rownames(mat.test),ont = "BP")
    saveRDS(enrich.res,paste0("~/analysis/R/DEseq2/",celltype,"_normal_vs_disordered.enrichGO.rds"))
  }
  
  
  ggplot(result)+geom_point(aes(x=log2FoldChange,y=-log10(padj),color=significant))+
    geom_text_repel(data=result[result$significant=="yes",],aes(x=log2FoldChange,y=-log10(padj),label=gene))+
    geom_hline(yintercept=2,linetype=2)+
    geom_vline(xintercept=1,linetype=2)+
    geom_vline(xintercept=-1,linetype=2)+
    scale_color_manual(values = c("grey","red"))+theme_classic2()+NoLegend()+
    ggtitle(paste0("Normal vs Disordered, ",celltype))+theme(plot.title = element_text(hjust=0.5))
  ggsave(paste0("~/analysis/R/DEseq2/",celltype,"_normal_vs_disordered.pdf"),width = 6,height = 7)
  saveRDS(result,paste0("~/analysis/R/DEseq2/",celltype,"_normal_vs_disordered.rds"))
  saveRDS(dds.res,paste0("~/analysis/R/DEseq2/",celltype,"_normal_vs_disordered.DESeq2.obj.rds"))
}

test.deseq<-readRDS('~/analysis/R/DEseq2/CD14_Mono_normal_vs_disordered.DESeq2.obj.rds')
test.enrichGO<-readRDS('~/analysis/R/DEseq2/CD14_Mono_normal_vs_disordered.enrichGO.rds')
DESeq2DotPlot(test.deseq)
enrichGODotPlot(test.enrichGO)

# RNA velosity
meta2d.all<-NULL
for (file in list.files('~/analysis/R/FineTypeLevelUnsplicedPseudobulk/',pattern = "^meta2d*")) {
  message(file)
  this.file=file.path('~/analysis/R/FineTypeLevelUnsplicedPseudobulk',file)
  this.data=read.delim(this.file)
  this.data=this.data[!is.na(this.data$meta2d_pvalue),]
  if(nrow(this.data)==0){
    next
  }
  this.data$individual=getField(file,"_",2)
  this.data$celltype=getField(file,"_",-c(1,2)) %>% gsub(".txt","",.,fixed = T)
  if(is.null(meta2d.all)){
    meta2d.all=this.data
  }else{
    meta2d.all=rbind(this.data,meta2d.all)
  }
}
saveRDS(meta2d.all,'~/analysis/core_data/meta2d.pseudobulkByCelltypeByIndividual.unspliced.rds')

spliced.meta2d0<-readRDS('~/analysis/core_data/meta2d.pseudobulkByCelltypeByIndividual.spliced.rds')
spliced.meta2d<-filterMeta2dTable(spliced.meta2d0)

unspliced.meta2d0<-readRDS('~/analysis/core_data/meta2d.pseudobulkByCelltypeByIndividual.unspliced.rds')
unspliced.meta2d<-filterMeta2dTable(unspliced.meta2d0)

holistic.meta2d0<-readRDS('~/analysis/core_data/meta2d.pseudobulkByCelltypeByIndividual.rds')
holistic.meta2d<-filterMeta2dTable(holistic.meta2d0)

spliced.genes<-unique(spliced.meta2d$CycID)
unspliced.genes<-unique(unspliced.meta2d$CycID)
holistic.genes<-unique(holistic.meta2d$CycID)
my_list <- list(spliced = spliced.genes, unspliced = unspliced.genes, holistic = holistic.genes)

plotVenn(my_list)

spliced.summa<-spliced.meta2d %>% group_by(CycID,individual) %>% summarise(.,phase=median(meta2d_phase))
colnames(spliced.summa)[3]<-"spliced"
unspliced.summa<-unspliced.meta2d %>% group_by(CycID,individual) %>% summarise(.,phase=median(meta2d_phase))
colnames(unspliced.summa)[3]<-"unspliced"
phases<-left_join(spliced.summa,unspliced.summa,by=c("CycID","individual"))
phases<-phases[!is.na(phases$spliced)&!is.na(phases$unspliced),]
phases$diff<-phases$unspliced-phases$spliced
phases$spliced<-ifelse(phases$diff>12,phases$spliced+24,phases$spliced)
phases$spliced<-ifelse(phases$diff< -12,phases$spliced-24,phases$spliced)
phases$diff<-NULL
phases<-gather(phases,key="status",value="phase",-CycID,-individual)


ggplot(phases,aes(x=status,y=phase,fill=status))+geom_boxplot()+
  scale_x_discrete(limits=c("unspliced","spliced"))+
  stat_compare_means(method = "wilcox.test",comparisons = list(c("unspliced", "spliced")),label = "p.signif")+
  theme_half_open()+theme(axis.text.x=element_text(angle=60,hjust=1))+NoLegend()+ylab("peaking time")+xlab("")

# generate matrix for detection of protein rhythmicity
srt.cytof<-readRDS('~/analysis/core_data/CYTOF.300k.reannotation.srt.rds')
srt.cytof$type<-srt.cytof$manual.new
mat<-srt2bulkMatrix(srt.cytof,split.by = c("individual","CT","type"),bool.use.CPHC = T)
saveRDS(mat,"~/analysis/core_data/CYTOF.5individual.pseudobulk.mat.rds")

meta2d.all<-NULL
for (file in list.files('~/analysis/R/FineTypeLevelCYTOFPseudobulk/',pattern = "^meta2d*")) {
  message(file)
  this.file=file.path('~/analysis/R/FineTypeLevelCYTOFPseudobulk',file)
  this.data=read.delim(this.file)
  this.data=this.data[!is.na(this.data$meta2d_pvalue),]
  if(nrow(this.data)==0){
    next
  }
  this.data$individual=getField(file,"_",2)
  this.data$celltype=getField(file,"_",-c(1,2)) %>% gsub(".txt","",.,fixed = T)
  if(is.null(meta2d.all)){
    meta2d.all=this.data
  }else{
    meta2d.all=rbind(this.data,meta2d.all)
  }
}
saveRDS(meta2d.all,'~/analysis/core_data/meta2d.pseudobulkByCelltypeByIndividual.cytof.rds')


# explore cortisol regulated genes
meta2d.result0<-readRDS('~/analysis/core_data/JTK.result.filtered.16individual.addp2t.addfold2bkg.addmedianexp.bytype.rds')
meta2d.result<-filterMeta2dTable(meta2d.result,float.pvalue.s1 = 0.05,float.pvalue.s2 = 0.05,int.min.individual = 2,mode = "JTK")
cortisol.experiment.DE.genes0<-read.delim('~/analysis/core_data/GSE148639_processed.txt',header = T,skip = 9)
cortisol.experiment.DE.genes0<-cortisol.experiment.DE.genes0[1:4]
cortisol.experiment.DE.genes<-dplyr::filter(cortisol.experiment.DE.genes0,p.value<0.05&abs(Fold.change)>1&Gene.Symbol!="---")
gene_vector<-cortisol.experiment.DE.genes$Gene.Symbol
cortisol.regulated.genes<-unlist(strsplit(gene_vector, " /// "))
all.result=NULL
for(this.individual in unique(meta2d.result$individual)){
  this.meta2d.result=meta2d.result[meta2d.result$individual==this.individual,]
  for(this.celltype in unique(this.meta2d.result$celltype)){
    this.features=(this.meta2d.result[this.meta2d.result$celltype==this.celltype,] %>% arrange(ADJ.P) %>% head(50))$CycID
    if(this.celltype=="CD4 TCM"){
      print(this.features[this.features %in% cortisol.regulated.genes])
    }
    ratio.cortisol.regulated=sum(this.features %in% cortisol.regulated.genes)/length(this.features)
    this.result=data.frame("individual"=this.individual,"celltype"=this.celltype,
                           "ratio"=ratio.cortisol.regulated,"n.oscillating"=length(this.features))
    if(is.null(all.result)){
      all.result=this.result
    }else{
      all.result=rbind(all.result,this.result)
    }
  }
}
 