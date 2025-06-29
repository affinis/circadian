
#1. top level annotation pipeline
i=1
sample_path<-paste0(PATIENTS[9],"_",i,"/outs/per_sample_outs/",PATIENTS[9],"_",i,"/count/sample_filtered_feature_bc_matrix/")
rds_path<-paste0("~/data/analysis/",PATIENTS[9],"_",i,".levelTop.rds")
srt<-seuratWrap1(sample_path,min.features = 500)
srt<-seuratWrap2(srt)
srt<-seuratWrap3(srt,res = 0.01)

FeaturePlot(srt,c("GNLY","CD3D","JCHAIN","MS4A1","CD14"))+DimPlot(srt,label = T)
markers<-FindAllMarkers(srt,logfc.threshold = log(2),min.diff.pct = 0.4)

srt<-TypeCluster(srt,c(0,1,3,4,5,7,8,11),"TNK")
srt<-TypeCluster(srt,c(2),"Myeloid")
srt<-TypeCluster(srt,c(6,9,10),"B")
DimPlot(srt,label = T,group.by = "type")
saveRDS(srt,rds_path)

#1.1 using individual HZD as an example
srt<-seuratWrap1('TFSH190500A_HZD_1/outs/per_sample_outs/TFSH190500A_HZD_1/count/sample_filtered_feature_bc_matrix/')
srt<-seuratWrap2(srt)
srt<-seuratWrap3(srt,res = 0.01)
DimPlot(srt,label = T)
FeaturePlot(srt,c("nCount_RNA","CD3D","JCHAIN","MS4A1","CD14","CD34"),max=10000)
markers<-FindAllMarkers(srt,logfc.threshold = log(2),min.diff.pct = 0.4)

srt<-TypeCluster(srt,c(0,1,4,5,6,7),"TNK")
srt<-TypeCluster(srt,2,"Myeloid")
srt<-TypeCluster(srt,c(3,8),"B")
srt<-TypeCluster(srt,c(9,11),"Plasma")
srt<-TypeCluster(srt,10,"Stem")
DimPlot(srt,label = T,group.by = "type")
saveRDS(srt,"~/data/analysis/TFSH190500A_HZD_1.levelTop.rds")


system('gzip TFSH190500A_HZD_2/outs/per_sample_outs/TFSH190500A_HZD_2/count/sample_filtered_feature_bc_matrix/barcodes.tsv')
srt<-seuratWrap1('TFSH190500A_HZD_2/outs/per_sample_outs/TFSH190500A_HZD_2/count/sample_filtered_feature_bc_matrix/')
srt<-seuratWrap2(srt)
srt<-seuratWrap3(srt,res = 0.01)
DimPlot(srt,label = T)
FeaturePlot(srt,c("nCount_RNA","CD3D","JCHAIN","MS4A1","CD14","CD34"),max=10000)
markers<-FindAllMarkers(srt,logfc.threshold = log(2),min.diff.pct = 0.4)

srt<-TypeCluster(srt,c(0,1,4,5,6,7),"TNK")
srt<-TypeCluster(srt,9,"Plasma")
srt<-TypeCluster(srt,c(3,8),"B")
srt<-TypeCluster(srt,c(2,10,12,13),"Myeloid")
srt<-TypeCluster(srt,11,"Stem")
DimPlot(srt,label = T,group.by = "type")
saveRDS(srt,"~/data/analysis/TFSH190500A_HZD_2.levelTop.rds")

system('gzip TFSH190500A_HZD_3/outs/per_sample_outs/TFSH190500A_HZD_3/count/sample_filtered_feature_bc_matrix/barcodes.tsv')
srt<-seuratWrap1('TFSH190500A_HZD_3/outs/per_sample_outs/TFSH190500A_HZD_3/count/sample_filtered_feature_bc_matrix/')
srt<-seuratWrap2(srt)
srt<-seuratWrap3(srt,res = 0.01)
DimPlot(srt,label = T)
FeaturePlot(srt,c("nCount_RNA","CD3D","JCHAIN","MS4A1","CD14","CD34"),max=10000)
markers<-FindAllMarkers(srt,logfc.threshold = log(2),min.diff.pct = 0.4)

srt<-TypeCluster(srt,c(0,1,3,5,6,8,10),"TNK")
srt<-TypeCluster(srt,c(2,11),"Myeloid")
srt<-TypeCluster(srt,4,"B")
srt<-TypeCluster(srt,7,"Plasma")
srt<-TypeCluster(srt,9,"Stem")
DimPlot(srt,label = T,group.by = "type")
saveRDS(srt,"~/data/analysis/TFSH190500A_HZD_3.levelTop.rds")

system('gzip TFSH190500A_HZD_4/outs/per_sample_outs/TFSH190500A_HZD_4/count/sample_filtered_feature_bc_matrix/barcodes.tsv')
srt<-seuratWrap1('TFSH190500A_HZD_4/outs/per_sample_outs/TFSH190500A_HZD_4/count/sample_filtered_feature_bc_matrix/')
srt<-seuratWrap2(srt)
srt<-seuratWrap3(srt,res = 0.01)
DimPlot(srt,label = T)
FeaturePlot(srt,c("nCount_RNA","CD3D","JCHAIN","MS4A1","CD14","CD34"),max=10000)
markers<-FindAllMarkers(srt,logfc.threshold = log(2),min.diff.pct = 0.4)

srt<-TypeCluster(srt,c(0,1,4,5,6,7),"TNK")
srt<-TypeCluster(srt,8,"Plasma")
srt<-TypeCluster(srt,3,"B")
srt<-TypeCluster(srt,c(2,9),"Myeloid")
srt<-TypeCluster(srt,10,"Stem")
DimPlot(srt,label = T,group.by = "type")
saveRDS(srt,"~/data/analysis/TFSH190500A_HZD_4.levelTop.rds")


system('gzip TFSH190500A_HZD_5/outs/per_sample_outs/TFSH190500A_HZD_5/count/sample_filtered_feature_bc_matrix/barcodes.tsv')
srt<-seuratWrap1('TFSH190500A_HZD_5/outs/per_sample_outs/TFSH190500A_HZD_5/count/sample_filtered_feature_bc_matrix/')
srt<-seuratWrap2(srt)
srt<-seuratWrap3(srt,res = 0.01)
DimPlot(srt,label = T)
FeaturePlot(srt,c("nCount_RNA","CD3D","JCHAIN","MS4A1","CD14","CD34"),max=10000)
markers<-FindAllMarkers(srt,logfc.threshold = log(2),min.diff.pct = 0.4)

srt<-TypeCluster(srt,c(0,1,4,5,6,7),"TNK")
srt<-TypeCluster(srt,c(8,9),"Plasma")
srt<-TypeCluster(srt,3,"B")
srt<-TypeCluster(srt,c(2),"Myeloid")
srt<-TypeCluster(srt,10,"Stem")
DimPlot(srt,label = T,group.by = "type")
saveRDS(srt,"~/data/analysis/TFSH190500A_HZD_5.levelTop.rds")

system('gzip TFSH190500A_HZD_6/outs/per_sample_outs/TFSH190500A_HZD_6/count/sample_filtered_feature_bc_matrix/barcodes.tsv')
srt<-seuratWrap1('TFSH190500A_HZD_6/outs/per_sample_outs/TFSH190500A_HZD_6/count/sample_filtered_feature_bc_matrix/')
srt<-seuratWrap2(srt)
srt<-seuratWrap3(srt,res = 0.01)
DimPlot(srt,label = T)
FeaturePlot(srt,c("nCount_RNA","CD3D","JCHAIN","MS4A1","CD14","CD34"),max=10000)
markers<-FindAllMarkers(srt,logfc.threshold = log(2),min.diff.pct = 0.4)

srt<-TypeCluster(srt,c(0,1,3,5,6,8),"TNK")
srt<-TypeCluster(srt,c(7,12),"Plasma")
srt<-TypeCluster(srt,4,"B")
srt<-TypeCluster(srt,c(2,10,11),"Myeloid")
srt<-TypeCluster(srt,9,"Stem")
DimPlot(srt,label = T,group.by = "type")
saveRDS(srt,"~/data/analysis/TFSH190500A_HZD_6.levelTop.rds")

#2. transform mouse circadian genes to human's
mouse2human<-readRDS('../list/mouse2human.rds')
circadian_gene_mm<-read.delim('../list/mouse_circadian_genes.tsv',header = F)
colnames(circadian_gene_mm)<-c("gene_name","gene_id","pos","description","CT_RNAseq","CT_HDA")
mapping<-mouse2human[,c(1,3)]
colnames(mapping)<-c("gene_name","gene_name_human")
circadian_gene_mm<-left_join(circadian_gene_mm,mapping,by="gene_name")
circadian_gene_hs<-circadian_gene_mm[,c("CT_RNAseq","gene_name_human")]
circadian_gene_hs<-circadian_gene_hs[!is.na(circadian_gene_hs$gene_name_human),]
colnames(circadian_gene_hs)<-c("CT","gene_name")

test0<-getWholeDaySplicedData(PATIENTS[1])
test<-mergeSplicedData(test0,normalize = T,num.scale = 100000,normalize_by_splicing_stat = T)
all_data=NULL
i=1
for (this.gene in circadian_gene_hs$gene_name) {
  message(this.gene)
  if(!(this.gene %in% rownames(test$spliced))){
    i=i+1
    next
  }
  subdata=plotSingleSplicingCounts(test,this.gene,return.data = T,summarise = "sum",normalize = F)
  subdata$scaled_count=NA
  subdata$total_count=log(subdata$total_count+1)
  subdata[subdata$status=="spliced","scaled_count"]=(subdata[subdata$status=="spliced","total_count"] %>% scale(.))[,1]
  subdata[subdata$status=="unspliced","scaled_count"]=(subdata[subdata$status=="unspliced","total_count"] %>% scale(.))[,1]
  subdata$gene=this.gene
  if(is.null(all_data)){
    all_data=subdata
    message("New assign")
  }else{
    message("rbinding")
    all_data=rbind(all_data,subdata)
  }
  message(paste0(i,"/",length(circadian_gene_hs$gene_name)))
  i=i+1
}
all_data %>% saveRDS(.,"TFSH190500A_HZD.circadianInPub.rds")

all_data$status<-paste0(all_data$status,all_data$time)
colnames(all_data)[6]<-"feature"
all_data<-left_join(all_data,circadian_gene_hs,by="feature") %>% arrange(.,by=phase)
circadian_gene_hs<-arrange(circadian_gene_hs,by=phase)


colnames(circadian_gene_hs)<-c("phase","feature")
wanted_genes<-intersect(circadian_gene_hs$feature,rownames(test0[[1]]$spliced))
df<-NULL
for(i in 1:length(time_point)){
  this.total=sum(test0[[i]]$spliced)+sum(test0[[i]]$unspliced)+sum(test0[[i]]$ambiguous)
  df.this.spliced=test0[[i]]$spliced[wanted_genes,] %>% rowSums() %>% as.data.frame() %>% `colnames<-`(.,paste0("spliced_",time_point[i]))
  #this.spliced.total=sum(test0[[i]]$spliced)
  #df.this.spliced[,paste0("spliced_",time_point[i])]=df.this.spliced[,paste0("spliced_",time_point[i])]/this.spliced.total
  df.this.spliced[,paste0("spliced_",time_point[i])]=df.this.spliced[,paste0("spliced_",time_point[i])]/this.total
  
  df.this.unspliced=test0[[i]]$unspliced[wanted_genes,] %>% rowSums() %>% as.data.frame() %>% `colnames<-`(.,paste0("unspliced_",time_point[i]))
  #this.unspliced.total=sum(test0[[i]]$unspliced)
  #df.this.unspliced[,paste0("unspliced_",time_point[i])]=df.this.unspliced[,paste0("unspliced_",time_point[i])]/this.unspliced.total
  df.this.unspliced[,paste0("unspliced_",time_point[i])]=df.this.unspliced[,paste0("unspliced_",time_point[i])]/this.total
  
  df.this=cbind(df.this.spliced,df.this.unspliced)
  if(is.null(df)){
    df=df.this
  }else{
    df=cbind(df,df.this)
  }
}

df<-df*100000
df<-df %>% as.matrix() %>% proportions(.,2)
df<-df*100000
df<-as.data.frame(df)

if(T){
  spliced.order=c("spliced_1","spliced_5","spliced_9",
             "spliced_13","spliced_17","spliced_21")
  unspliced.order=c("unspliced_1","unspliced_5","unspliced_9",
                  "unspliced_13","unspliced_17","unspliced_21")
  df.spliced<-df[,spliced.order]
  df.unspliced<-df[,unspliced.order]

  df.spliced<-scale(t(df.spliced)) %>% as.data.frame()
  df.unspliced<-scale(t(df.unspliced)) %>% as.data.frame()
  
  order<-c("spliced_1","unspliced_1","spliced_5","unspliced_5","spliced_9","unspliced_9",
           "spliced_13","unspliced_13","spliced_17","unspliced_17","spliced_21","unspliced_21")
  df<-rbind(df.spliced,df.unspliced)[order,]
  df<-t(df) %>% as.data.frame()
}



df$feature<-rownames(df)
df<-left_join(df,circadian_gene_hs,by="feature")
rownames(df)<-df$feature
df[is.na(df)]<-0
ncmdplot.mod(df,1:12,"phase",minors.threshold = NULL,show.obname = F,cluster.by.y = T,
            features.sort = order)



#3. try to merge seurat object directly and take a look at core circadian gene expression
test<-mergeSamplesOneDay(PATIENTS[1])
plotGeneExpressionByCT(test,cell.type = "B",alpha = 0.1)
plotGeneExpressionByCT(test,cell.type = "Myeloid",alpha = 0.1,k=6)
plotGeneExpressionByCT(test,cell.type = "TNK",alpha = 0.1,k=6)
plotGeneExpressionByCT(test,cell.type = "Plasma",alpha = 0.1)


data=data.frame("features"=character(0),"time"=character(0),"cell_type"=character(0),"expression"=numeric(0))
for(celltype in c("B","T","Myeloid","NK","Plasma")){
  thisdata=fetchMergedDataOneCellType(test,gene.list=c("CLOCK","ARNTL","PER1","PER2","CRY1","CRY2","NR1D1","NR1D2","RORC","DBP"),cell.type = celltype)
  thisdata=thisdata %>% group_by(.,features,time,cell_type) %>% summarise(.,expression=median(values))
  data=rbind(data,thisdata)
}
timesave<-data$time
data<-data %>% group_by(.,features,cell_type) %>% summarise(.,scaled_data=scale(expression))
data$time=timesave
data$scaled_data<-data$scaled_data %>% as.vector()
data$time<-factor(data$time,levels=c("CT1","CT5","CT9","CT13","CT17","CT21"))
ggplot(data)+geom_raster(aes(x=time,y=cell_type,fill=scaled_data))+
  scale_fill_gradientn(colors = rev(brewer.pal(11,"RdBu")))+
  facet_wrap(vars(features),ncol=11)+theme(axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle=45,hjust=1))+xlab("")

#4. code from RNA velocyto, test velocyto
testdata<-getWholeDaySplicedData("TFSH190500A_HZD","B")
test<-mergeSplicedData(testdata)

emat<-test$spliced
emat<-emat[,colSums(emat)>=1e3]
emat<-emat[!duplicated(rownames(emat)),]
r<-Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)
r$adjustVariance(plot=T,do.par=T,gam.k=10)
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine')
r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)
r$getEmbedding(type='PCA',embeddingType='UMAP',perplexity=50,verbose=T)
r$getEmbedding(type='PCA',embeddingType='largeVis',perplexity=50,verbose=T)
r$getEmbedding(type='PCA',embeddingType='FR',perplexity=50,verbose=T)
r$getEmbedding(type='PCA',embeddingType='UMAP_graph',perplexity=50,verbose=T)

r$plotEmbedding(type='PCA',embeddingType='largeVis',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters')
r$plotEmbedding(type='PCA',embeddingType='UMAP',colors=r$counts[,"MS4A1"],main='MS4A1')
r$plotEmbedding(type='PCA',embeddingType='UMAP',colors=r$counts[,"CD3D"],main='CD3D')
r$plotEmbedding(type='PCA',embeddingType='largeVis',colors=r$counts[,"CD3D"],main='CD3D')

emat <- test$spliced
nmat <- test$unspliced
emat <- emat[,rownames(r$counts)]
nmat <- nmat[,rownames(r$counts)] # restrict to cells that passed p2 filter
# take cluster labels
cluster.label<-r$clusters$PCA[[1]]
values_label<-names(cluster.label) %>% strsplit(.,"-") %>% lapply(.,`[`,2) %>% unlist()%>% as.vector()
names(values_label)<-names(cluster.label)
values_label<-factor(values_label,levels = CT_TIME)
cell.colors <- sccore::fac2col(values_label)
# take embedding
emb <- r$embeddings$PCA$UMAP_graph
cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
shared<-intersect(rownames(emat),rownames(nmat))
shared.cell<-intersect(colnames(emat),colnames(nmat))
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat[shared,shared.cell],nmat[shared,shared.cell],deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)
show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),
                               cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)


#5. plot cell ratio versus time
srts<-getSamplesOneDay("TFSH190500A_HZD")
plotCellRatioVersusCT(srts)+NoLegend()

#no integration is needed, merge is enough, as the plots following show 
test<-mergeSamplesOneDay(PATIENTS[1])
test<-NormalizeData(test)
test<-FindVariableFeatures(test)
test<-ScaleData(test)
test<-RunPCA(test)
test<-FindNeighbors(test, dims = 1:30, reduction = "pca")
test<-FindClusters(test, resolution = 2, cluster.name = "unintegrated_clusters")
test<-RunUMAP(test, dims = 1:30, reduction = "pca", reduction.name = "umap")
test<-RunTSNE(test,reduction = "pca",dims = 1:30)
dimplot_publication(test, reduction = "tsne", group.by="type")
dimplot_publication(test, reduction = "tsne", group.by="CT",colors = pal_cosmic()(8))
Idents(test)<-test$CT
DimPlot(test,cells.highlight=CellsByIdentities(test,"CT21"),reduction = "tsne")+theme_dimplot()+ggtitle("CT21")+NoLegend()

#6. try to create dimensional reduction with circadian genes
test<-mergeSamplesOneDay(PATIENTS[1])
test<-NormalizeData(test)
test<-ScaleData(test)
test<-FindVariableFeatures(test)
test<-RunPCA(test,reduction.name="pca.circadian",features=circadian_gene_hs$gene_name)
ElbowPlot(test,reduction = "pca.circadian")
VizDimLoadings(test,reduction = "pca.circadian",dims = 1:4)+theme_publication(base.size = 10)
test<-FindNeighbors(test, dims = 1:10, reduction = "pca.circadian")
test<-FindClusters(test, resolution = 0.1, cluster.name = "pca.circadian.clusters")
test<-RunUMAP(test, dims = 1:10, reduction = "pca.circadian", reduction.name = "umap.circadian")
DimPlot(test,reduction="umap.circadian")
dimplot_publication(test,reduction.name="umap.circadian",reduction="umap",group.by = "CT")
dimplot_publication(test,reduction.name="umap.circadian",reduction="umap",group.by = "type",colors = pal_aaas()(8))
p<-DimPlot(test)
occilating_cells<-CellSelector(p)

testdata<-getWholeDaySplicedData("TFSH190500A_HZD")
testdata<-mergeSplicedData(testdata)
subtestdata<-subsetLoomMatrix(testdata,cells=occilating_cells,genes=circadian_gene_hs$gene_name)
#these cells did not pass filtering 
viewRNAVelosity(subtestdata)

#
test$CT<-factor(test$CT,levels=CT_TIME_ORDER)
Idents(test)<-test$CT
dimplot_publication(subtest,group.by = "CT",reduction.name="umap.circadian",reduction ="umap")
allTNK<-subset(test,type=="TNK")
allTNK@meta.data$isOccilating<-ifelse(rownames(allTNK@meta.data) %in% occilating_cells,"occilating","non-occilating")
VlnPlot(allTNK,c("nCount_RNA","nFeature_RNA","percent.mt"),group.by="isOccilating",alpha = 0)
allTNK.j<-JoinLayers(allTNK)
occilating_vs_non<-FindMarkers(allTNK.j,ident.1="occilating",ident.2="non-occilating",group.by="isOccilating")

subtest<-subset(allTNK,isOccilating=="occilating")
subtest2<-subset(allTNK,isOccilating=="non-occilating")

dimplot_publication(subtest,reduction.name="umap.circadian",reduction="umap",group.by = "CT")
plotGeneExpressionByCT(subtest,cell.type = "TNK",k=6)
occilating_genes<-subtest@reductions$pca.circadian@feature.loadings %>% as.data.frame() %>% rownames()

#7. try to find markers by CT
test<-mergeSamplesOneDay(PATIENTS[1])
test<-JoinLayers(test)
test<-NormalizeData(test)
test<-ScaleData(test)
CTmarkers1<-FindMarkers(test,CT_TIME_ORDER[1],CT_TIME_ORDER[4],group.by="CT")
CTmarkers2<-FindMarkers(test,CT_TIME_ORDER[2],CT_TIME_ORDER[5],group.by="CT")
CTmarkers3<-FindMarkers(test,CT_TIME_ORDER[3],CT_TIME_ORDER[6],group.by="CT")
CTmarkers1_names<-CTmarkers1 %>% dplyr::filter(.,p_val_adj<0.01,abs(avg_log2FC)>1.5) %>% rownames()
CTmarkers2_names<-CTmarkers2 %>% dplyr::filter(.,p_val_adj<0.01,abs(avg_log2FC)>1.5) %>% rownames()
CTmarkers3_names<-CTmarkers3 %>% dplyr::filter(.,p_val_adj<0.01,abs(avg_log2FC)>1.5) %>% rownames()
plotGeneExpressionByCT(test,CTmarkers3_names[73:96],"B",k = 6,alpha = 0.01)
plotGeneExpressionByCT(test,CTmarkers3_names[73:96],"Myeloid",alpha = 0.01,k = 6)
plotGeneExpressionByCT(test,CTmarkers3_names[73:96],"TNK",k = 6,alpha = 0.01)


#8. integration of total 36 data and make annotation manually to level2
test<-mergeAllSamples()
saveRDS(test,"../analysis/all.samples.merged.rds")
test<-dimentionalReductionMergedSrt(test)
dimplot_publication(test,"patient")
saveRDS(test,"../analysis/all.samples.merged.red.rds")

test[["RNA"]] <- split(test[["RNA"]], f = test$patient)
test<-IntegrateLayers(object=test, method=CCAIntegration, orig.reduction="pca", new.reduction="integrated.cca",verbose = FALSE)
saveRDS(test,"../analysis/all.samples.integrated.unjoined.rds")

test[["RNA"]] <- JoinLayers(test[["RNA"]])
test <- FindNeighbors(test, reduction = "integrated.cca", dims = 1:30)
test <- FindClusters(test, resolution = 0.1)
test <- RunUMAP(test, dims = 1:30, reduction = "integrated.cca")
saveRDS(test,"../analysis/all.samples.integrated.umapran.rds")

test<-readRDS('../analysis/all.samples.integrated.umapran.rds')
dimplot_publication(test,group.by = "patient")
dimplot_publication(test,group.by = "seurat_clusters",label = T)
FeaturePlot(test,c("CD3D","GNLY","JCHAIN","MS4A1","FCN1","CD34","PPBP","CCR3","LILRA4"))
test<-TypeCluster(test,c(0,2,3,6,10,12),type="TNK",new.meta="type.top.level.manual")
test<-TypeCluster(test,c(1,5,8,11),type="Myeloid",new.meta="type.top.level.manual")
test<-TypeCluster(test,c(4,9),type="B",new.meta="type.top.level.manual")
test<-TypeCluster(test,c(7),type="Platelet",new.meta="type.top.level.manual")
dimplot_publication(test,group.by = "type.top.level.manual")
saveRDS(test,"../analysis/all.samples.integrated.toplevel.rds")

#TNK
test<-readRDS("../analysis/all.samples.integrated.toplevel.rds")
TNK.test<-subset(test,type.top.level.manual=="TNK")
remove(test)
gc()
TNK.test<-integrateSubset(TNK.test)
saveRDS(TNK.test,"../analysis/all.samples.TNK.integrated.rds")
TNK.test<-readRDS("../analysis/all.samples.TNK.integrated.rds")
TNK.test<-FindClusters(TNK.test,resolution=1)
TNK.test <- RunAzimuth(TNK.test, reference = "pbmcref")
dimplot_publication(TNK.test,label=T,reduction = "int",colors = generateColor(35,col.dist.min=0.1))
FeaturePlot(TNK.test,c("CD3D","GNLY","JCHAIN","MS4A1","FCN1","CD34","PPBP","CCR3","LILRA4"),order = T)
FeaturePlot(TNK.test,c("CD3D","CD4","CD8A","CCR7","GNLY","MKI67","SLC4A10","PTGER2","FOXP3"))
FeaturePlot(TNK.test,"nFeature_RNA",max.cutoff = 4000)
#BiocManager::install('multtest')
#install.packages('metap')
TNK.conserved<-FindAllMarkers(TNK.test,min.pct = 0.7,min.diff.pct = 0.5)
i=25
dplyr::filter(TNK.conserved,cluster==i) %>% head()
TNK.test@meta.data[TNK.test@meta.data$seurat_clusters==i,"predicted.celltype.l2"] %>% table()
TNK.test@meta.data[TNK.test@meta.data$seurat_clusters==i,"predicted.celltype.l2.score"] %>% density() %>% plot()
TNK.test<-TypeCluster(TNK.test,c(0,6,25,33),type="cTNK01_CD4+CCR7+Tn",new.meta="type.bottom.level.manual")
TNK.test<-TypeCluster(TNK.test,c(1),type="cTNK02_CD8+CCR7+Tn",new.meta="type.bottom.level.manual")
TNK.test<-TypeCluster(TNK.test,c(2,5,8,10,18,27,29),type="cTNK03_CD4+PTGER2+Tcm",new.meta="type.bottom.level.manual")
TNK.test<-TypeCluster(TNK.test,c(3,9,14,16,24,28,32),type="cTNK04_CD8+GZMK+Tem",new.meta="type.bottom.level.manual")
TNK.test<-TypeCluster(TNK.test,c(4,7,11,22),type="cTNK05_GNLY+NK",new.meta="type.bottom.level.manual")
TNK.test<-TypeCluster(TNK.test,c(12),type="cTNK06_gdT",new.meta="type.bottom.level.manual")
TNK.test<-TypeCluster(TNK.test,c(13),type="cTNK07_CD8+SLC4A10+MAIT",new.meta="type.bottom.level.manual")
TNK.test<-TypeCluster(TNK.test,c(15),type="cTNK08_CD4+NKG7+CTL",new.meta="type.bottom.level.manual")
TNK.test<-TypeCluster(TNK.test,c(17,20,23),type="debris",new.meta="type.bottom.level.manual")
TNK.test<-TypeCluster(TNK.test,c(21),type="cTNK09_CD4+FOXP3+Treg",new.meta="type.bottom.level.manual")
TNK.test<-TypeCluster(TNK.test,c(26),type="cTNK10_MKI67+proliferating",new.meta="type.bottom.level.manual")
TNK.test<-TypeCluster(TNK.test,c(30),type="cTNK11_LYST+dnT",new.meta="type.bottom.level.manual")
TNK.test<-TypeCluster(TNK.test,c(19),type="cTNK12_CD4+cTRBV+Tcm-like",new.meta="type.bottom.level.manual")
TNK.test<-TypeCluster(TNK.test,c(31),type="cTNK13_CD8+cTRBV+Tcm-like",new.meta="type.bottom.level.manual")

TNK.test<-TypeCluster(TNK.test,c(0:16,18:19,21,22,23),type="TNK",new.meta="type.top.level.manual")
TNK.test<-TypeCluster(TNK.test,c(17,20,23),type="debris",new.meta="type.top.level.manual")
dimplot_publication(TNK.test,label=T,reduction = "int",group.by = "type.top.level.manual")
dimplot_publication(TNK.test,label=T,reduction = "int",group.by = "type.bottom.level.manual")
saveRDS(TNK.test,"../analysis/all.samples.TNK.integrated.bottom.level.rds")

test<-readRDS("../analysis/all.samples.integrated.toplevel.rds")
Myeloid.test<-subset(test,type.top.level.manual=="Myeloid")
remove(test)
gc()
Myeloid.test<-integrateSubset(Myeloid.test)
saveRDS(Myeloid.test,"../analysis/all.samples.Myeloid.integrated.rds")
Myeloid.test<-readRDS("../analysis/all.samples.Myeloid.integrated.rds")
dimplot_publication(Myeloid.test,label=T,reduction = "int")
dimplot_publication(Myeloid.test,group.by = "patient",reduction = "int")
FeaturePlot(Myeloid.test,c("CD3D","GNLY","MS4A1","FCN1","PPBP","LILRA4","CD14","FCGR3A","nFeature_RNA"))
FeaturePlot(Myeloid.test,"FCER1A")
Myeloid.test <- RunAzimuth(Myeloid.test, reference = "pbmcref")
Myeloid.conserved<-FindAllMarkers(Myeloid.test,min.pct = 0.7,min.diff.pct = 0.5)
i=13
dplyr::filter(Myeloid.conserved,cluster==i) %>% head()
Myeloid.test@meta.data[Myeloid.test@meta.data$seurat_clusters==i,"predicted.celltype.l2"] %>% table()
Myeloid.test@meta.data[Myeloid.test@meta.data$seurat_clusters==i,"predicted.celltype.l2.score"] %>% density() %>% plot()
Myeloid.test<-TypeCluster(Myeloid.test,c(0,1,5),type="cM01_CD14+FCGR3A-Monocyte",new.meta="type.bottom.level.manual")
Myeloid.test<-TypeCluster(Myeloid.test,c(2),type="cM02_CD14-FCGR3A+Monocyte",new.meta="type.bottom.level.manual")
Myeloid.test<-TypeCluster(Myeloid.test,c(3,6,7,10),type="doublet",new.meta="type.bottom.level.manual")
Myeloid.test<-TypeCluster(Myeloid.test,c(4),type="cM04_FCER1A+cDC2",new.meta="type.bottom.level.manual")
Myeloid.test<-TypeCluster(Myeloid.test,c(8),type="debris",new.meta="type.bottom.level.manual")
Myeloid.test<-TypeCluster(Myeloid.test,c(9,13),type="cM05_LILRA4+pDC",new.meta="type.bottom.level.manual")
Myeloid.test<-TypeCluster(Myeloid.test,c(11),type="cS01_CD34+HSPC",new.meta="type.bottom.level.manual")
Myeloid.test<-TypeCluster(Myeloid.test,c(12),type="cM03_CLEC9A+cDC1",new.meta="type.bottom.level.manual")
dimplot_publication(Myeloid.test,label=T,reduction = "int",group.by = "type.bottom.level.manual")
Myeloid.test<-TypeCluster(Myeloid.test,c(0,1,2,4,5,9,12,13),type="Myeloid",new.meta="type.top.level.manual")
Myeloid.test<-TypeCluster(Myeloid.test,c(11),type="Stem",new.meta="type.top.level.manual")
Myeloid.test<-TypeCluster(Myeloid.test,c(8),type="debris",new.meta="type.top.level.manual")
Myeloid.test<-TypeCluster(Myeloid.test,c(3,6,7,10),type="doublet",new.meta="type.top.level.manual")
dimplot_publication(Myeloid.test,label=T,reduction = "int",group.by = "type.top.level.manual")
saveRDS(Myeloid.test,"../analysis/all.samples.Myeloid.integrated.bottom.level.rds")

test<-readRDS("../analysis/all.samples.integrated.toplevel.rds")
B.test<-subset(test,type.top.level.manual=="B")
remove(test)
gc()
B.test<-integrateSubset(B.test)
saveRDS(B.test,"../analysis/all.samples.B.integrated.rds")
dimplot_publication(B.test,label=T,reduction = "int")
dimplot_publication(B.test,group.by = "patient",reduction = "int")
FeaturePlot(B.test,c("CD3D","GNLY","JCHAIN","MS4A1","FCN1","CD34","PPBP","CCR3","LILRA4","CD14"))
FeaturePlot(B.test,"JCHAIN")
B.conserved<-FindAllMarkers(B.test,min.pct = 0.7,min.diff.pct = 0.5)
B.test <- RunAzimuth(B.test, reference = "pbmcref")
i=1
dplyr::filter(B.conserved,cluster==i) %>% head()
B.test@meta.data[B.test@meta.data$seurat_clusters==i,"predicted.celltype.l2"] %>% table()
B.test<-TypeCluster(B.test,c(3),type="cB01_TCL1A+B_naive",new.meta="type.bottom.level.manual")
B.test<-TypeCluster(B.test,c(0,1,2),type="cB03_IgG+IgA+B_memory",new.meta="type.bottom.level.manual")
B.test<-TypeCluster(B.test,c(4,9,11),type="cB04_Plasma",new.meta="type.bottom.level.manual")
B.test<-TypeCluster(B.test,c(10),type="cB02_EGR1+B_intermediate",new.meta="type.bottom.level.manual")
B.test<-TypeCluster(B.test,c(5,6,7),type="doublet",new.meta="type.bottom.level.manual")
B.test<-TypeCluster(B.test,c(5,6,7),type="doublet",new.meta="type.top.level.manual")
B.test<-TypeCluster(B.test,c(8),type="debris",new.meta="type.bottom.level.manual")
B.test<-TypeCluster(B.test,c(8),type="debris",new.meta="type.top.level.manual")
dimplot_publication(B.test,label=T,reduction = "int",group.by = "type.bottom.level.manual")
dimplot_publication(B.test,label=T,reduction = "int",group.by = "type.top.level.manual")
saveRDS(B.test,"../analysis/all.samples.B.integrated.bottom.level.rds")

generateAnnotationFile(TNK.file = '../analysis/all.samples.TNK.integrated.bottom.level.rds',
                       B.file = '../analysis/all.samples.B.integrated.bottom.level.rds',
                       Myeloid.file = '../analysis/all.samples.Myeloid.integrated.bottom.level.rds',
                       cols = c("type.top.level.manual","type.bottom.level.manual"),
                       new.colnames = c("manual.level1","manual.level2"))

#get only B/Myeloid/TNK cells
annotation.meta<-read.delim('../analysis/cell.annotations.tsv',header = T)
rownames(annotation.meta)<-annotation.meta$cell_id
test<-readRDS('../analysis/all.samples.integrated.toplevel.rds')
test<-AddMetaData(test,metadata=annotation.meta)
test<-subset(test,manual.level1 %in% c("TNK","B","Myeloid","Stem"))
dimplot_publication(test,group.by = "manual.level2",reduction = "int",colors = generateColor(24,alpha = 0.5))
test<-RunTSNE(test)
dimplot_publication(test,group.by = "manual.level2",reduction = "tsne",reduction.name = "tsne",colors = generateColor(24))
saveRDS(test,"../analysis/all.samples.oddsremoved.tsneran.rds")

test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.rds")
test<-RunAzimuth(test, reference = "pbmcref")
test$predicted.celltype.l1.5<-case_when(test$predicted.celltype.l2=="ASDC" ~ "ASDC",
                                        test$predicted.celltype.l2=="B intermediate" ~ "B memory",
                                        test$predicted.celltype.l2=="B memory" ~ "B memory",
                                        test$predicted.celltype.l2=="B naive" ~ "B naive",
                                        test$predicted.celltype.l2=="CD14 Mono" ~ "CD14 Mono",
                                        test$predicted.celltype.l2=="CD16 Mono" ~ "CD16 Mono",
                                        test$predicted.celltype.l2 %in% c("CD4 CTL","CD4 Naive","CD4 Proliferating","CD4 TCM","CD4 TEM","Treg") ~ "CD4 T",
                                        test$predicted.celltype.l2 %in% c("CD8 Naive","CD8 Proliferating","CD8 TCM","CD8 TEM") ~ "CD8 T",
                                        test$predicted.celltype.l2=="cDC1" ~ "cDC1",
                                        test$predicted.celltype.l2=="cDC2" ~ "cDC2",
                                        test$predicted.celltype.l2=="dnT" ~ "dnT",
                                        test$predicted.celltype.l2=="Eryth" ~ "Eryth",
                                        test$predicted.celltype.l2=="gdT" ~ "gdT",
                                        test$predicted.celltype.l2=="HSPC" ~ "HSPC",
                                        test$predicted.celltype.l2=="ILC" ~ "ILC",
                                        test$predicted.celltype.l2=="MAIT" ~ "MAIT",
                                        test$predicted.celltype.l2 %in% c("NK","NK Proliferating","NK_CD56bright") ~ "NK",
                                        test$predicted.celltype.l2=="pDC" ~ "pDC",
                                        test$predicted.celltype.l2=="Plasmablast" ~ "Plasma",
                                        test$predicted.celltype.l2=="Platelet" ~ "Platelet"
                                        )
saveRDS(test,"../analysis/all.samples.oddsremoved.tsneran.Azimuth.rds")

test<-subset(test,patient %in% HEALTH)
saveRDS(test,"../analysis/all.samples.oddsremoved.tsneran.Azimuth.healthonly.rds")

plotCellRatioVersusCT(list(test),type.col.name="manual.level2",colors = generateColor(24))+NoLegend()
plotdata<-test@meta.data[,c("cell_id","CT","patient","manual.level1","manual.level2")]
plotdata<-plotdata %>% group_by(.,manual.level2,CT,patient) %>% summarise(.,count=n())
plotdata$CT<-factor(plotdata$CT,levels = CT_TIME_ORDER)
total.count<-plotdata %>% group_by(CT,patient) %>% summarise(total=sum(count))
plotdata<-left_join(plotdata,total.count,by=c("CT","patient"))
plotdata$percentage<-plotdata$count*100/plotdata$total
plotdata$time_numeric<-plotdata$CT %>% gsub("CT","",.) %>% as.numeric()
ggplot(plotdata)+geom_point(aes(x=time_numeric,y=percentage,group=patient,color=patient))+geom_smooth(aes(x=time_numeric,y=percentage))+
  facet_wrap(~manual.level2,scales = "free_y")+scale_color_aaas()+theme_publication()+NoLegend()+xlab("CT")+ylab("percentage(%)")

#9. plot heatmap of core circadian gene expression versus time
plotIntegratedGeneExpressionHeatmapByCT(test)
plotIntegratedGeneExpressionHeatmapByCT(subset(test,patient %in% EMI_PATIENTS))
plotIntegratedGeneExpressionHeatmapByCT(subset(test,patient %in% HEALTH))
plotIntegratedGeneExpressionHeatmapByCT(subset(test,patient==HEALTH[3]))
plotIntegratedGeneExpressionHeatmapByCT(subset(test,patient!=HEALTH[3]),method = "mean")

test$type<-test$manual.level2
celltypes<-test$type %>% table %>% names()
test.sub<-subset(test,patient!=PATIENTS[4])
plotGeneExpressionByCT(test.sub,cell.type=celltypes[11],CT_field = 3,k=6,patient_filed = 1,alpha = 0.01)
ggplot(plotdata,aes(x=time_numeric,y=values,color=patient,group=patient))+
  geom_point(alpha=0.1)+geom_smooth(method="gam",formula=y~s(x,bs="cs",k=5))+scale_color_d3()+NoLegend()+facet_wrap(~features)

############################
# circadian gene prediction
############################
#1. trial of linear model
lm.test.data<-FetchData(subset(normal.test,manual.level2==celltypes[6]),vars="DNTTIP1")
lm.test.data<-dplyr::filter(lm.test.data,DNTTIP1>0)
lm.test.data$CT<-rownames(lm.test.data) %>% strsplit(.,"_") %>% lapply(.,`[`,3) %>% unlist() %>% gsub("CT","",.) %>% as.numeric()
colnames(lm.test.data)<-c("expression","CT")
lm.model<-geneIsOscillated(lm.test.data)
lm.model$rawdata
lm.model$residual
lm.model$normality_of_residual

#2. using metacycle to detect circadian genes
test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.rds")
celltypes<-test$manual.level2 %>% table() %>% names()
test$type<-test$manual.level2
for(i in 1:length(celltypes)){
  usecelltype<-celltypes[i]
  detectCircadianGenes(test,usecelltype,individuals=PATIENTS,out.dir = '../analysis/OSgenePredict_Highly_expressed_genes_3000/',HVarGene = F)
}


LS.result<-read.delim(paste0('../analysis/OSgenePredict/LSresult_',celltypes[14],".txt"))
JTK.result<-read.delim(paste0('../analysis/OSgenePredict/JTKresult_',celltypes[15],".txt"))
dplyr::filter(LS.result,p<0.2)
dplyr::filter(JTK.result,ADJ.P<0.05)
candidate_genes<-dplyr::filter(LS.result,p<0.3)$CycID
#test.data<-test %>% subset(.,patient %in% HEALTH[1:2])
plotGeneExpressionByCT(test.data,gene.list=candidate_genes,cell.type = celltypes[14],CT_field = 3,patient_filed = c(1,2),alpha=0.05)

#plot gene expression as clock
celltypes<-test$manual.level2 %>% table() %>% names()
plotCircadianGeneCount(0.2,celltypes,generateColor(25))+ylim(0,40)
plotCircadianGeneCount(0.2,celltypes,generateColor(25),out.dir = "../analysis/OSgenePredict_EMI/")+ylim(0,40)
plotCircadianGeneCount(0.2,celltypes,generateColor(25),out.dir = "../analysis/OSgenePredict_HEALTH/")+ylim(0,40)

plotPhaseofCellType(0.2,celltypes,generateColor(25),out.dir = "../analysis/OSgenePredict/")
plotPhaseofCellType(0.2,celltypes,generateColor(25),out.dir = "../analysis/OSgenePredict_HEALTH/",method = "LS")


test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.rds")
test$type<-test$manual.level2
celltypes<-test$type %>% table() %>% names()


celltypes<-test$type %>% table() %>% names()
normal.test<-subset(test,patient %in% HEALTH[1:2])
normal.test2<-subset(test,patient %in% HEALTH)
plotGeneExpressionByCT(test,cell.type=celltypes[8],CT_field=3,patient_filed=c(1,2),use.median = F,layer="count")+NoLegend()

plotClockByGenes(normal.test,cell.type=celltypes[2],gene.list=dplyr::filter(JTK.result,ADJ.P<0.05)$CycID[1:12])
plotCircadianGeneCount(0.05,celltypes,generateColor(25))+ylim(c(0,90))
plotCircadianGeneCount(0.05,celltypes,generateColor(25),out.dir = '../analysis/OSgenePredict_HEALTH/')+ylim(c(0,125))

plotClockByGenes(normal.test,cell.type=celltypes[15],gene.list=dplyr::filter(JTK.result,ADJ.P<0.05)$CycID[1:12])
plotPhaseofCellType(0.05,celltypes,generateColor(25),out.dir = "../analysis/OSgenePredict_HEALTH/",embeding = "box")

#take a look at tregs
test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.rds")
celltypes<-test$manual.level2 %>% table() %>% names()
for(i in 1:length(celltypes)){
  usecelltype<-celltypes[i]
  detectCircadianGenes(test,usecelltype,individuals=EMI_PATIENTS,features=CIRCADIAN_GENES_MAIN,out.dir = '../analysis/OSgenePredict_CIRCADIAN_GENES_PATIENT/')
}
VlnPlot(test,features = CIRCADIAN_GENES_MAIN[4],group.by = "manual.level2",cols =generateColor(25))+NoLegend()+
  theme(axis.text.x = element_text(angle=60,hjust=1),plot.margin = margin(l=1,unit = "cm"))+xlab("")
test.treg<-subset(test,type=="cTNK11_CD4+FOXP3+T_Treg")
test.treg<-subset(test.treg,patient %in% PATIENTS[1:2])
VlnPlot(test.treg,features = CIRCADIAN_GENES_MIAN,group.by = "type")
plotClockByGenes(test.treg,cell.type="cTNK11_CD4+FOXP3+T_Treg")
plotRhythmicityPvalue(celltypes=celltypes,colors = generateColor(25),method = "JTK",out.dir = "../analysis/OSgenePredict_CIRCADIAN_GENES_PATIENT/")

#celltype change
test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.rds")
saveRDS(test,"../analysis/all.samples.oddsremoved.tsneran.rds")
test$type<-test$manual.level2
DotPlot(test,SUBSET_MARKERS,group.by = "type")+theme_publication()
plotdata<-test@meta.data
plotdata<-plotdata %>% group_by(.,CT,patient,type) %>% summarise(.,ncell.patient=n())
plotdata2<-plotdata %>% group_by(.,patient,CT) %>% summarise(.,ncell.total=sum(ncell.patient))
plotdata_all<-left_join(plotdata,plotdata2,by=c("patient","CT"))
plotdata_all$percentage<-plotdata_all$ncell.patient*100/plotdata_all$ncell.total
plotdata_all$group<-ifelse(plotdata_all$patient %in% EMI_PATIENTS,"EMI","healthy")
ggplot(plotdata_all)+geom_boxplot(aes(x=CT,y=percentage,fill=group),color="black")+facet_wrap(~type,scales="free_y",ncol=6)+
  scale_color_d3()+scale_x_discrete(limits=CT_TIME_ORDER)

#take a look at MT1 and MT2 receptor
MTdata<-FetchData(test,vars=c("MTNR1A","MTNR1B"))

#reanalysis of spliced data
spliced.test<-getWholeDaySplicedData('TFSH190500A_HZD')
spliced.test.merge.raw<-mergeSplicedData(spliced.test)
plotSplicingCounts(spliced.test.merge.raw)

spliced.test.merge.norm<-mergeSplicedData(spliced.test,normalize=T)
plotSplicingCounts(spliced.test.merge.norm)

#take a look at Treg
Tregs<-Cells(subset(test,manual.level2=="cTNK11_CD4+FOXP3+T_Treg"&patient==PATIENTS[1]))
plotSplicingCounts(spliced.test.merge.raw,cells = Tregs)
plotSplicingCounts(spliced.test.merge.norm,cells = Tregs)

spliced.treg.merge.raw<-subsetLoomMatrix(spliced.test.merge.raw,cells = Tregs)
spliced.treg.merge.norm<-subsetLoomMatrix(spliced.test.merge.norm,cells = Tregs)

viewRNAVelosity(spliced.treg.merge.raw,min.transcripts.per.cell=1000,min.cells.per.gene=3)
viewRNAVelosity(spliced.treg.merge.norm,min.transcripts.per.cell=1000,min.cells.per.gene=3)

#try to use circadian genes discovered to do rna velosity
#not feasible
Tregs<-Cells(subset(test,manual.level2=="cTNK11_CD4+FOXP3+T_Treg"&patient==HEALTH[1]))
circadian.genes<-read.delim('../analysis/OSgenePredict/LSresult_cTNK11_CD4+FOXP3+T_Treg.txt') %>% dplyr::filter(.,p<0.2)
circadian.spliced.treg.raw<-subsetLoomMatrix(spliced.test.merge.raw,cells = Tregs,genes = circadian.genes$CycID)
viewRNAVelosity(circadian.spliced.treg.raw,min.transcripts.per.cell=1,min.cells.per.gene=1)
plotSplicingCounts(circadian.spliced.treg.raw,cells = Tregs,features = circadian.genes$CycID)


#circadian genes in PMID_38190520
test$type<-test$manual.level2
plotClockByGenes(subset(test,patient %in% HEALTH[1:3]),gene.list = CIRCADIAN_GENES_PMID_38190520,
                 cell.type="cTNK11_CD4+FOXP3+T_Treg",CT_field = 3,patient_filed = c(1,2),summa="sum")


#using mean to examine spliced data
splice_test<-getWholeDaySplicedData("TFSH190500A_HZD")
splice_test<-mergeSplicedData(splice_test)
plotSplicingCounts(splice_test)
plotSplicedVersusUnspliced(splice_test,summarise = "mean")

#view core circadian genes in Treg
Tregs<-Cells(subset(test,manual.level2=="cTNK11_CD4+FOXP3+T_Treg"&patient==HEALTH[1]))
plotSplicedVersusUnspliced(splice_test,summarise = "sum",cells =Tregs )
plotSplicedVersusUnspliced(splice_test,summarise = "mean",cells =Tregs )
circadian.genes<-read.delim('../analysis/OSgenePredict/LSresult_cTNK11_CD4+FOXP3+T_Treg.txt') %>% dplyr::filter(.,p<0.1)
plotSplicedVersusUnspliced(splice_test,summarise = "sum",cells =Tregs,features=circadian.genes$CycID)
plotSplicedVersusUnspliced(splice_test,summarise = "mean",cells =Tregs,features=circadian.genes$CycID)


#load all spliced data
all_splice0=list()
for(patient in PATIENTS){
  this.data=getWholeDaySplicedData(patient)
  all_splice0=c(all_splice0,this.data)
}

#saveRDS(all_splice0,file = "../analysis/all.sample.rnaVelosity.rds")
all_splice0<-readRDS("../analysis/all.sample.rnaVelosity.rds")
all_splice<-mergeSplicedData(all_splice0)


Tregs_health<-Cells(subset(test,manual.level2=="cTNK11_CD4+FOXP3+T_Treg"&patient %in% HEALTH))
#shared_Tregs_health<-intersect(Tregs_health,all_splice$spliced %>% colnames())
plotSplicingCounts(all_splice,summarise = "mean",cells = Tregs_health)
plotSplicedVersusUnspliced(all_splice,summarise = "mean",cells = Tregs_health)

celltypes<-test$manual.level2 %>% unique() %>% sort()
using_type<-20
cells_health<-Cells(subset(test,manual.level2==celltypes[using_type]&patient %in% HEALTH))
p1<-plotSplicingCounts(all_splice,summarise = "mean",cells = cells_health)+ggtitle(celltypes[using_type])
p2<-plotSplicedVersusUnspliced(all_splice,summarise = "mean",cells = cells_health)
ggarrange(p2,p1)


splice_table<-plotSplicedVersusUnspliced(all_splice,summarise ="mean",
              cells = cells_health,return.data=T,features = NULL)
splice_table<-splice_table %>% arrange(.,feature,time_numeric)
gene_CT_list<-split(splice_table,splice_table$feature)
results<-isPolygonParallel(gene_CT_list)
results<-results %>% unlist()
candidate_genes_by_rnaVelosity_test<-names(gene_CT_list)[results]

#saveRDS(candidate_genes_by_rnaVelosity,"../analysis/Treg.velocity.rhythmic.rds")


#splice_table<-plotSplicedVersusUnspliced(all_splice,summarise ="mean",cells = cells_health,return.data=T,features = NULL)
splice_table<-splice_table %>% arrange(.,feature,time_numeric)
gene_CT_list<-split(splice_table,splice_table$feature)
results<-isConvexParallel(gene_CT_list)
results<-results %>% unlist()
candidate_genes_by_rnaVelosity2<-names(gene_CT_list)[results]
candidate_genes_by_rnaVelosity<-intersect(candidate_genes_by_rnaVelosity,candidate_genes_by_rnaVelosity2)
p1<-plotSplicingCounts(all_splice,summarise="mean",cells = cells_health,
                       features=candidate_genes_by_rnaVelosity)+ggtitle(celltypes[using_type])
p2<-plotSplicedVersusUnspliced(all_splice,summarise = "mean",cells = cells_health,
                               features=candidate_genes_by_rnaVelosity)
ggarrange(p2,p1)

saveRDS(candidate_genes_by_rnaVelosity,"../analysis/Treg.velocity.rhythmic.convex.rds")


#using spliced data to find oscillating genes
test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.rds")
all_splice0<-readRDS("../analysis/all.sample.rnaVelosity.rds")
all_splice<-mergeSplicedData(all_splice0)
remove(all_splice0)
celltypes<-test$manual.level2 %>% unique() %>% sort()
for(celltype in celltypes[1:24]){
  all.patient.velosity.rhythm<-data.frame("type"=character(0),"individual"=character(0),"feature"=character(0))
  message(celltype)
  for(individual in PATIENTS){
    message(individual)
    message("fetching cell data")
    subsrt=subset(test,manual.level2==celltype&patient==individual)
    cells_this=Cells(subsrt)
    #subsrt[["RNA"]]=split(subsrt[["RNA"]], f = subsrt$CT)
    subsrt=NormalizeData(subsrt)
    #subsrt=FindVariableFeatures(subsrt,nfeatures=3000)
    genes_wanted=rowSums(subsrt) %>% sort() %>% tail(.,3000) %>% names()
    message("fetching splicing table")
    #splice_table=plotSplicedVersusUnspliced(all_splice,summarise ="mean",normalize=F,
    #                  cells = cells_this,return.data=T,features = VariableFeatures(subsrt))
    splice_table=plotSplicedVersusUnspliced(all_splice,summarise ="mean",normalize=T,
                      cells = cells_this,return.data=T,features = genes_wanted)
    splice_table=splice_table %>% arrange(.,feature,time_numeric)
    message("splicing table")
    gene_CT_list=split(splice_table,splice_table$feature)
    message("determining rhythmicity")
    results=isConvexParallel(gene_CT_list)
    results=results %>% unlist()
    candidate_genes_by_rnaVelosity=names(gene_CT_list)[results]
    message("done")
    this.data.frame=data.frame("type"=rep(celltype,length(candidate_genes_by_rnaVelosity)),
                               "individual"=rep(individual,length(candidate_genes_by_rnaVelosity)),
                               "feature"=candidate_genes_by_rnaVelosity)
    all.patient.velosity.rhythm=rbind(all.patient.velosity.rhythm,this.data.frame)
    gc(reset = T)
  }
  if(grepl("Temra/Teff",celltype)){
    celltype=gsub("Temra/Teff","Temra#Teff",celltype)
  }else{
    celltype=celltype
  }
  saveRDS(all.patient.velosity.rhythm,paste0("/lustre/home/acct-medll/medll/data/analysis/velosity_rhythm_highly_expressed_genes_3000_by_type_norm/",celltype,".rds"))
  gc(reset = T)
}


using_type<-20
spliceB<-readRDS('../analysis/velosity_rhythm_variable_genes_by_type/cTNK11_CD4+FOXP3+T_Treg.rds')
(spliceB$feature %>% table())[spliceB$feature %>% table()>=2] %>% length()
test_genes<-(spliceB$feature %>% table())[spliceB$feature %>% table()>=2] %>% names()
using_genes<-test_genes[1]
using_genes<-"NEAT1"
#celltypes<-test$manual.level2 %>% unique() %>% sort()

using_patients<-spliceB[spliceB$feature==using_genes,"individual"]
cells_health<-Cells(subset(test,manual.level2==celltypes[using_type]&patient %in% using_patients))
p1<-plotSingleSplicingCounts(all_splice,summarise = "mean",cells = cells_health,feature =using_genes )+ggtitle(paste0(celltypes[using_type],": ",using_genes))
p2<-plotSingleSplicedVersusUnspliced(all_splice,summarise = "mean",feature =using_genes,cells = cells_health)
ggarrange(p1,p2,nrow = 2)

#view genes oscillating with splicing data
all_table=NULL
for(file in list.files('../analysis/velosity_rhythm_variable_genes_3000_by_type_normalize_scale10000/') %>% grep("^c",.,value = T)){
  file_table=readRDS(paste0('../analysis/velosity_rhythm_variable_genes_3000_by_type_normalize_scale10000/',file))
  if(is.null(all_table)){
    all_table=file_table
  }else{
    all_table=rbind(all_table,file_table)
  }
}
all_table_join<-all_table %>% group_by(.,feature) %>% summarise(.,n.observation=n())
all_table_bytype<-all_table[c("type","feature")] %>% unique() %>% group_by(.,feature) %>% summarise(.,n.celltype=n())
all_table_byindividual<-all_table[c("individual","feature")] %>% unique() %>% group_by(.,feature) %>% summarise(.,n.individual=n())
all_table_join<-left_join(all_table_join,all_table_bytype,by="feature")
all_table_join<-left_join(all_table_join,all_table_byindividual,by="feature")
all_table_join$mean.observation<-all_table_join$n.observation/all_table_join$n.celltype


plot.table<-all_table[,c(1,3)] %>% unique()
ggplot(plot.table)+geom_bar(aes(x=type,fill=type),alpha=0.75)+
  scale_fill_manual(values=generateColor(24))+
  theme(axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  panel.background = element_rect(fill="white"))


plot.table2<-all_table %>% group_by(.,type,feature) %>% summarise(n.observaion=n()) %>% 
  dplyr::filter(n.observaion>=2)

ggplot(plot.table2)+geom_bar(aes(x=type,fill=type),alpha=0.75)+
  scale_fill_manual(values=generateColor(24))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill="white"))


#take a look at phase of enriched genes
genes<-all_table_join[all_table_join$n.individual>=2,"feature"] %>% unlist() %>% as.vector()
gene_ids<-bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID
out<-enrichGO(gene_ids,OrgDb = org.Hs.eg.db)
barplot(out)
interest_id<-out@result$geneID[1] %>% strsplit(.,"/") %>% unlist()
interest_name<-bitr(interest_id,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db)$SYMBOL
j=9
individuals<-all_table[all_table$type==celltypes[i]&all_table$feature==interest_name[j],"individual"]
cells<-Cells(subset(test,manual.level2==celltypes[i]&patient %in% individuals))
p1<-plotSingleSplicingCounts(all_splice,summarise = "mean",cells = "cTNK05_CD8+CCR7+T_CD8+Tn",feature="AAK1")
p2<-plotSingleSplicedVersusUnspliced(all_splice,summarise = "mean",feature =interest_name[j],cells = cells)
ggarrange(p1,p2,nrow = 2)

#GO gene enrichment analysis
phase_data<-NULL
for(celltype in celltypes){
  message(celltype)
  genes=all_table_bytypefeature[all_table_bytypefeature$type==celltype,"feature"] %>% unlist() %>% as.vector()
  gene_ids=bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID
  GOout=enrichGO(gene_ids,OrgDb = org.Hs.eg.db)
  if(is.null(GOout)){
    next
  }
  interest_id=GOout@result[GOout@result$p.adjust<=0.05,"geneID"] %>% strsplit(.,"/")
  interest_pathway=GOout@result[GOout@result$p.adjust<=0.05,"Description"]
  if(length(interest_id)==0){
    next
  }
  #interest_name=bitr(interest_id,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db)$SYMBOL
  for(i in 1:length(interest_pathway)){
    interest_name=bitr(interest_id[[i]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db)$SYMBOL
    for(feature in interest_name){
      individuals=all_table[all_table$type==celltype&all_table$feature==feature,"individual"]
      cells=Cells(subset(test,manual.level2==celltype&patient %in% individuals))
      data=fetchSingleSplicedCounts(all_splice,feature=feature,cells=cells)
      data_max=data %>% group_by(.,individual) %>% slice(which.max(mean_count))
      data_max$feature=feature
      data_max$celltype=celltype
      data_max$pathway=interest_pathway[i]
      data_max$label="max"
      data_min=data %>% group_by(.,individual) %>% slice(which.min(mean_count))
      data_min$feature=feature
      data_min$celltype=celltype
      data_min$pathway=interest_pathway[i]
      data_min$label="min"
      data=rbind(data_max,data_min)
      if(is.null(phase_data)){
        phase_data=data
      }else{
        phase_data=rbind(phase_data,data)
      }
    }
  }
}
phase_data$time<-factor(phase_data$time,levels=CT_TIME_ORDER)
saveRDS(phase_data,"../analysis/velosity_rhythm_variable_genes_by_type/phasedata_total.rds")
phase_data<-readRDS("../analysis/velosity_rhythm_variable_genes_by_type/phasedata_total.rds")

phase_data_sumar<-phase_data %>% group_by(.,time,celltype,pathway,label,individual) %>% summarise(n.observation=n())
using_type<-17

ggplot(phase_data_sumar[phase_data_sumar$celltype %in% celltypes[using_type],])+
  geom_point(aes(x=time,y=pathway,size=n.observation,shape=label))+facet_wrap(~individual)+
  scale_color_d3(palette="category20")+scale_shape_manual(values=c(16,1))+ggtitle(celltypes[using_type])

#take a look at proportion of core circadian genes
test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.Azimuth.rds")
counts.main<-FetchData(test,vars=c(CIRCADIAN_GENES_MAIN,"CXCR3","ITGAL","ITGB2","GZMB","MT-CO1"),layer = "count")
meta<-test@meta.data[,c("nCount_RNA","manual.level2")]
counts.main<-mergeByRownames(counts.main,meta)
counts.main$CLOCK<-counts.main$CLOCK*10000/counts.main$nCount_RNA
counts.main$ARNTL<-counts.main$ARNTL*10000/counts.main$nCount_RNA
counts.main$PER1<-counts.main$PER1*10000/counts.main$nCount_RNA
counts.main$PER2<-counts.main$PER2*10000/counts.main$nCount_RNA
counts.main$PER3<-counts.main$PER3*10000/counts.main$nCount_RNA
counts.main$CRY1<-counts.main$CRY1*10000/counts.main$nCount_RNA
counts.main$CRY2<-counts.main$CRY2*10000/counts.main$nCount_RNA
counts.main$NR1D1<-counts.main$NR1D1*10000/counts.main$nCount_RNA
counts.main$NR1D2<-counts.main$NR1D2*10000/counts.main$nCount_RNA

counts.main$CXCR3<-counts.main$CXCR3*10000/counts.main$nCount_RNA
counts.main$ITGAL<-counts.main$ITGAL*10000/counts.main$nCount_RNA
counts.main$ITGB2<-counts.main$ITGB2*10000/counts.main$nCount_RNA
counts.main$GZMB<-counts.main$GZMB*10000/counts.main$nCount_RNA
counts.main$`MT-CO1`<-counts.main$`MT-CO1`*10000/counts.main$nCount_RNA
counts.main$CT<-getField(rownames(counts.main),"_",3)
counts.main$individual<-getField(rownames(counts.main),"_",c(1,2))
counts.main<-counts.main %>% group_by(.,manual.level2,CT,individual) %>% summarise(.,mean_CLOCK=mean(CLOCK),
                mean_ARNTL=mean(ARNTL),mean_PER1=mean(PER1),mean_PER2=mean(PER2),
                mean_PER3=mean(PER3),mean_PER3=mean(PER3),mean_CRY1=mean(CRY1),
                mean_CRY2=mean(CRY2),mean_NR1D1=mean(NR1D1),mean_NR1D2=mean(NR1D2),
                mean_CXCR3=mean(CXCR3),mean_ITGAL=mean(ITGAL),mean_ITGB2=mean(ITGB2),mean_GZMB=mean(GZMB),
                "mean_MT-CO1"=mean(`MT-CO1`))
test_gene<-"MT-CO1"
ggplot(counts.main)+geom_point(aes(x=CT,y=get(paste0("mean_",test_gene)),color=manual.level2))+scale_x_discrete(limits=CT_TIME_ORDER)+
  facet_wrap(~individual)+scale_color_manual(values=generateColor(25))
ggplot(counts.main)+
  geom_boxplot(aes(x=CT,y=get(paste0("mean_",test_gene)),color=manual.level2))+
  scale_x_discrete(limits=CT_TIME_ORDER)+
  scale_color_manual(values=generateColor(25))+ylab("")+ggtitle(test_gene)

celltypes<-names(table(test$manual.level2))
JTK.result=NULL
for(celltype in celltypes){
  if(file.exists(paste0('../analysis/OSgenePredict_Highly_expressed_genes_3000/JTKresult_',celltype,".txt"))){
    JTK.result.part=read.delim(paste0('../analysis/OSgenePredict_Highly_expressed_genes_3000/JTKresult_',celltype,".txt"))
    JTK.result.part$celltype=celltype
  }else{
    next
  }
  if(is.null(JTK.result)){
    JTK.result=JTK.result.part
  }else{
    JTK.result=rbind(JTK.result,JTK.result.part)
  }
}
JTK.result.filtered<-JTK.result%>% filter(.,ADJ.P<0.05)
length(JTK.result.filtered$CycID)
all_table_join.filtered<-all_table_join %>% filter(.,n.observation>=2,n.celltype>=2)
nrow(all_table_join.filtered)

intersect(JTK.result.filtered$CycID,all_table_join.filtered$feature) %>% length()

#FIG. 1B
test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.Azimuth.rds")
test$type<-test$predicted.celltype.l1.5
dimplot_publication(test,reduction = "umap",reduction.name="umap",
                    colors=generateColor(30,col.dist.min=0.3,seed = 10),group.by = "type")

ggplot(test@meta.data)+geom_bar(aes(x=CT,fill=type),position="fill",color="black")+
  scale_fill_manual(values=generateColor(30,col.dist.min=0.3,seed=10))+
  scale_x_discrete(limits=rev(CT_TIME_ORDER))+scale_y_continuous(expand=c(0,0))+
  theme(axis.text = element_text(size=11),plot.margin=margin(0,0.5,0,0,"cm"))+
  NoLegend()+ylab("")+xlab("")+coord_flip()

#FIG. S1B
test$type<-test$predicted.celltype.l2
dimplot_publication(test,reduction = "tsne",reduction.name="tsne",
       colors=generateColor(30,col.dist.min=0.3,seed = 10),group.by = "type")
dimplot_publication(test,colors=generateColor(30,col.dist.min=0.3,seed = 10),group.by = "type")
ggplot(test@meta.data)+geom_bar(aes(x=CT,fill=type),position="fill",color="black")+
  scale_fill_manual(values=generateColor(30,col.dist.min=0.3,seed=10))+
  scale_x_discrete(limits=rev(CT_TIME_ORDER))+theme_minimal(base_size = 13)+
  NoLegend()+ylab("ratio")+xlab("")+coord_flip()

#FIG. 1C
meta<-test@meta.data
meta$CT<-factor(meta$CT,levels=CT_TIME_ORDER)
meta$type<-meta$predicted.celltype.l1.5
celltypes<-sort(unique(meta$type))
meta$type<-factor(meta$type,levels=rev(celltypes))


plotFIG1C1<-ggplot(meta)+geom_bar(aes(x=type,fill=CT),position="fill",color="black")+
  scale_fill_d3(alpha = 0.75)+scale_y_continuous(expand = c(0,0))+
  theme(plot.margin=margin(0.5,0,0,0,"cm"),
       axis.text.x=element_blank(),axis.ticks.x=element_blank(),
       axis.text.y=element_text(color="black"))+
  ylab("")+xlab("")+NoLegend()

plotFIG1C2<-DotPlot(test,rev(c("AXL","SIGLEC6",
               "MS4A1","CD27",
               "CD14","FCGR3A",
               "CD3D","CD4","CD8A","FOXP3",
               "CD1C","CLEC9A","HBA1",
               "TRGC1","TRDC","CD34","SLC4A10",
               "GNLY","GZMK","GZMH","LILRA4",
               "JCHAIN","PPBP")),group.by="type",
        cols = c("#507AAF","#BE5C37"))+scale_y_discrete(limits=sort(unique(test$type)))+
  scale_size_continuous(range = c(0,4))+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=9),
        plot.margin=margin(0,0,0,0,"cm"),axis.text.y=element_text(size=9),
        panel.border = element_rect(color = "black",linewidth = 1),
        legend.text = element_text(size=10),legend.title=element_text(size=10))+
  xlab("")+ylab("")+coord_flip()

ggarrange(plotFIG1C1,plotFIG1C2,nrow=2,align="v",heights = c(1,5))
  

#FIG. S1C
meta<-test@meta.data
meta$CT<-factor(meta$CT,levels=CT_TIME_ORDER)
celltypes<-sort(unique(meta$type))
meta$type<-factor(meta$type,levels=rev(celltypes))


ggplot(meta)+geom_bar(aes(x=type,fill=CT),position="fill",color="black")+
  scale_fill_aaas()+
  coord_flip()+theme(axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.y=element_text(size = 11))+ylab("")+xlab("")+NoLegend()
DotPlot(test,c("MS4A1","CD27","TCL1A","JCHAIN","EGR1",
               "CD14","FCGR3A","CD1C","LILRA4","CLEC9A",
               "CD3D","CD4","CCR7","CD8A","GNLY","GZMK","GZMH",
               "FCER1G","PTGER2","FGFBP2","TRGC1","TRDC",
               "FOXP3","SLC4A10","MKI67","MX1"),group.by="type",
        cols = c("#507AAF","#BE5C37"))+scale_y_discrete(limits=sort(unique(test$type),decreasing = T))+
  scale_size_continuous(range = c(0,4))+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=9),axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(color = "black",linewidth = 1),
        legend.title=element_text(size=11))+xlab("")+ylab("")

ggarrange(left,right,ncol = 2,align="h",widths = c(1:5),common.legend = T)

#FIG. S1A
sample_summary<-read.delim('../../script/circadian_sample_summary.tsv')
plotdataFIGS1A<-sample_summary[sample_summary$FEATURE %in% c("Estimated number of cells","Sequencing saturation",
                          "Median genes per cell","Confidently mapped to transcriptome"),]
rownames(plotdataFIGS1A)<-plotdataFIGS1A$FEATURE
plotdataFIGS1A$FEATURE<-NULL
plotdataFIGS1A<-t(plotdataFIGS1A)
plotdataFIGS1A<-TransformContinuousDf(plotdataFIGS1A,NULL)
plotdataFIGS1A$values<-gsub("[%,]","",plotdataFIGS1A$values) %>% as.numeric()

plotdataFIGS1A$individual<-getField(plotdataFIGS1A$observations,"_",c(1,2))
plotdataFIGS1A$time<-getField(plotdataFIGS1A$observations,"_",3)
plotdataFIGS1A$values<-as.numeric(plotdataFIGS1A$values)
ggplot(plotdataFIGS1A)+geom_bar(aes(x=observations,y=values,fill=individual),linewidth=0.2,color="black",stat="identity")+
  facet_wrap(~features,nrow = 4,scales = "free_y")+scale_fill_d3("category20")+ylab("")+xlab("")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "top")

#FIG. S1B
DimPlot(test,group.by = "type",reduction = "tsne")+theme_dimplot()+ggtitle("")
Idents(test)<-test$patient
DimPlot(test,cells.highlight=CellsByIdentities(test,idents=PATIENTS[6]),order = T,reduction = "tsne")+
  NoLegend()+theme_dimplot()

#FIG. 2A
plotdataFIG2A1<-test$predicted.celltype.l1.5 %>% table() %>% as.data.frame()
test$type<-test$predicted.celltype.l1.5
celltypes<-sort(unique(test$predicted.celltype.l1.5))
colnames(plotdataFIG2A1)<-c("cell_type","number_of_cell")
FIG2A1<-ggplot(plotdataFIG2A1)+geom_bar(aes(y=cell_type,x=log10(number_of_cell)),stat="identity",fill="blue",alpha=0.5,color="black")+
  geom_text(aes(y=cell_type,x=log10(number_of_cell)/2,label=number_of_cell),color="white")+
  scale_y_discrete(limits=rev(plotdataFIG2A1$cell_type))+scale_x_continuous(expand=c(0,0))+
  theme(panel.grid = element_blank(),panel.background=element_rect(fill="white"))+ylab("")+xlab("log10(# of cell)")

JTK_CYCL<-NULL
for(block_data in list.files("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5/")){
  block_path=paste0("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5/",block_data,"/")
  block=getJTK_CYCLEouts(block_path,celltypes=celltypes,filter.data = F)
  if(is.null(JTK_CYCL)){
    JTK_CYCL=block
  }else{
    JTK_CYCL=rbind(JTK_CYCL,block)
  }
}

plotdataFIG2A2<-test@meta.data[,c("nFeature_RNA","type")] %>% group_by(type) %>% summarise(median_genes=max(nFeature_RNA))
plotdataFIG2A2<-left_join(plotdataFIG2A2,(JTK_CYCL %>% dplyr::filter(.,ADJ.P<=0.01))$celltype %>% table() %>% as.data.frame() %>% `colnames<-`(.,c("type","rhythmic")))
plotdataFIG2A2$percent_rhythmic<-plotdataFIG2A2$rhythmic*100/plotdataFIG2A2$median_genes
plotdataFIG2A2$percent_rhythmic<-round(plotdataFIG2A2$percent_rhythmic,2)
plotdataFIG2A2[is.na(plotdataFIG2A2)]<-0

FIG2A2<-ggplot(plotdataFIG2A2)+geom_bar(aes(y=type,x=percent_rhythmic),stat="identity",fill="red",alpha=0.5,color="black")+
  geom_text(data=plotdataFIG2A2[plotdataFIG2A2$percent_rhythmic>=12,],aes(y=type,x=percent_rhythmic/2,label=percent_rhythmic),color="white")+
  geom_text(data=plotdataFIG2A2[plotdataFIG2A2$percent_rhythmic<12,],aes(y=type,x=percent_rhythmic+5,label=percent_rhythmic),color="black")+
  scale_y_discrete(limits=rev(plotdataFIG2A1$cell_type))+scale_x_continuous(expand=c(0,0))+
  theme(panel.grid = element_blank(),panel.background=element_rect(fill="white"),
        axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  ylab("")+xlab("percentage of\nrhythmic genes(%)")

plotdataFIG2A3<-JTK_CYCL[JTK_CYCL$CycID %in% CIRCADIAN_GENES_MAIN,]
plotdataFIG2A3$celltype<-factor(plotdataFIG2A3$celltype,levels=rev(celltypes))
FIG2A3<-ggplot(plotdataFIG2A3)+geom_point(aes(x=celltype,y=CycID,size=-log10(ADJ.P),color=log2(AMP)))+
  scale_color_gradientn(colours = c("#507AAF","#BE5C37"))+
  theme(axis.text.x = element_text(angle=45,hjust=1),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        panel.background = element_rect(fill="white"))+coord_flip()+
  xlab("")+ylab("rhythmicity of\ncore circadian genes")

ggarrange(FIG2A1,FIG2A2,FIG2A3,ncol=3,widths=c(1.5,1.1,2.7),align="h")


#FIG. S2A
plotdataFIGS2A1<-test$predicted.celltype.l2 %>% table() %>% as.data.frame()
celltypes<-test$predicted.celltype.l2 %>% unique() %>% sort()
test$type<-test$predicted.celltype.l2
colnames(plotdataFIGS2A1)<-c("cell_type","number_of_cell")
FIGS2A1<-ggplot(plotdataFIGS2A1)+geom_bar(aes(y=cell_type,x=log10(number_of_cell)),stat="identity",fill="blue",alpha=0.5,color="black")+
  geom_text(aes(y=cell_type,x=log10(number_of_cell)/2,label=number_of_cell),color="white")+
  scale_y_discrete(limits=rev(plotdataFIGS2A1$cell_type))+scale_x_continuous(expand=c(0,0))+
  theme(panel.grid = element_blank(),panel.background=element_rect(fill="white"))+ylab("")+xlab("log10(# of\ncell)")

JTK_CYCL<-NULL
for(block_data in list.files("../analysis/OSgenePredict_PseudoBulk_20k_all_genes/")){
  block_path=paste0("../analysis/OSgenePredict_PseudoBulk_20k_all_genes/",block_data,"/")
  block=getJTK_CYCLEouts(block_path,celltypes=celltypes,filter.data = F)
  if(is.null(JTK_CYCL)){
    JTK_CYCL=block
  }else{
    JTK_CYCL=rbind(JTK_CYCL,block)
  }
}

plotdataFIGS2A2<-test@meta.data[,c("nFeature_RNA","type")] %>% group_by(type) %>% summarise(median_genes=max(nFeature_RNA))
plotdataFIGS2A2<-left_join(plotdataFIGS2A2,(JTK_CYCL %>% dplyr::filter(.,ADJ.P<=0.01))$celltype %>% table() %>% as.data.frame() %>% `colnames<-`(.,c("type","rhythmic")))
plotdataFIGS2A2$percent_rhythmic<-plotdataFIGS2A2$rhythmic*100/plotdataFIGS2A2$median_genes
plotdataFIGS2A2$percent_rhythmic<-round(plotdataFIGS2A2$percent_rhythmic,2)
plotdataFIGS2A2[is.na(plotdataFIGS2A2)]<-0

FIGS2A2<-ggplot(plotdataFIGS2A2)+geom_bar(aes(y=type,x=percent_rhythmic),stat="identity",fill="red",alpha=0.5,color="black")+
  geom_text(data=plotdataFIGS2A2[plotdataFIGS2A2$percent_rhythmic>=20,],aes(y=type,x=percent_rhythmic/2,label=percent_rhythmic),color="white")+
  geom_text(data=plotdataFIGS2A2[plotdataFIGS2A2$percent_rhythmic<20,],aes(y=type,x=percent_rhythmic+10,label=percent_rhythmic),color="black")+
  scale_y_discrete(limits=rev(celltypes))+scale_x_continuous(expand=c(0,0))+
  theme(panel.grid = element_blank(),panel.background=element_rect(fill="white"),
        axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  ylab("")+xlab("percentage of\nrhythmic genes(%)")

plotdataFIGS2A3<-JTK_CYCL[JTK_CYCL$CycID %in% CIRCADIAN_GENES_MAIN,]
plotdataFIGS2A3$celltype<-factor(plotdataFIGS2A3$celltype,levels=rev(celltypes))
FIGS2A3<-ggplot(plotdataFIGS2A3)+geom_point(aes(x=celltype,y=CycID,size=-log10(ADJ.P),color=log2(AMP)))+
  scale_color_gradientn(colours = c("#507AAF","#BE5C37"))+
  theme(axis.text.x = element_text(angle=45,hjust=1),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        panel.background = element_rect(fill="white"))+coord_flip()+
  xlab("")+ylab("")

ggarrange(FIGS2A1,FIGS2A2,FIGS2A3,ncol=3,widths=c(1.75,1,3),align="h")

#FIG. 2B
plotdataFIG2B<-NULL
for(celltype in c("CD14 Mono","CD16 Mono","CD4 T")){
  plotdataFIG2B.part1=plotPseudobulkByCT(test,celltype,features=c("ARNTL","CRY1","NR1D1","PER1","TEF"),return.data = T,rep=2,proportion=0.5)
  plotdataFIG2B.part1$type=celltype
  plotdataFIG2B.part2=plotPseudobulkByCT(test,celltype,features=c("ARNTL","CRY1","NR1D1","PER1","TEF"),return.data = T,rep=2,proportion=0.5)
  plotdataFIG2B.part2$type=celltype
  plotdataFIG2B.part2$time=case_when(plotdataFIG2B.part2$time=="CT1"~"CT25",
                                     plotdataFIG2B.part2$time=="CT5"~"CT29",
                                     plotdataFIG2B.part2$time=="CT9"~"CT33",
                                     plotdataFIG2B.part2$time=="CT13"~"CT37",
                                     plotdataFIG2B.part2$time=="CT17"~"CT41",
                                     plotdataFIG2B.part2$time=="CT21"~"CT45")
  plotdataFIG2B.part=rbind(plotdataFIG2B.part1,plotdataFIG2B.part2)
  if(is.null(plotdataFIG2B)){
    plotdataFIG2B=plotdataFIG2B.part
  }else{
    plotdataFIG2B=rbind(plotdataFIG2B,plotdataFIG2B.part)
  }
}
summ.data<-plotdataFIG2B %>% group_by(feature,type) %>% summarise(min=min(count)+1)
plotdataFIG2B[plotdataFIG2B$count!=0,]
plotdataFIG2B<-left_join(plotdataFIG2B,summ.data,by=c("feature","type"))
plotdataFIG2B$relative_expression<-plotdataFIG2B$count/plotdataFIG2B$min
plotdataFIG2B.sup<-plotdataFIG2B %>% group_by(feature,type,time) %>% summarise(median_expression=median(relative_expression))

use_color<-generateColor(30,col.dist.min=0.3,seed = 10)
ggplot(plotdataFIG2B)+geom_boxplot(aes(x=time,y=relative_expression,color=type))+
  geom_point(data=plotdataFIG2B.sup,aes(x=time,y=median_expression,color=type))+
  geom_line(data=plotdataFIG2B.sup,aes(x=time,y=median_expression,group=type,color=type))+
  scale_color_manual(values=c(use_color[4],use_color[5],use_color[6],use_color[9]))+
  scale_x_discrete(limits=c(CT_TIME_ORDER,"CT25","CT29","CT33","CT37","CT41","CT45"))+
  facet_wrap(~feature,scales = "free_y",nrow=3)+xlab("")+
  theme(axis.text.x=element_text(angle=60,hjust=1,color="black",size=9),
        axis.text.y=element_text(color="black",size=9),legend.text=element_text(color="black",size=9),
        panel.background=element_rect(fill="white"),
        legend.position=c(0.75,0.1))

#FIG. 2C
JTK_CYCL.filtered<-NULL
celltypes<-test$type %>% unique()
for(block_data in list.files("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5/")){
  block_path=paste0("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5/",block_data,"/")
  block.filtered=getJTK_CYCLEouts(block_path,celltypes=celltypes,period = 24,padj = 0.01)
  if(is.null(JTK_CYCL.filtered)){
    JTK_CYCL.filtered=block.filtered
  }else{
    JTK_CYCL.filtered=rbind(JTK_CYCL.filtered,block.filtered)
  }
}
use_color<-generateColor(30,col.dist.min=0.3,seed = 10)[-11]
ggplot(JTK_CYCL.filtered[JTK_CYCL.filtered$celltype!="Eryth",])+
  geom_bar(aes(x=LAG,fill=celltype),alpha=0.5,color="black")+
  facet_wrap(~celltype,scales = "free_y",ncol=6)+
  scale_fill_manual(values=use_color)+NoLegend()+
  theme(axis.text=element_text(color="black"))+ylab("count of peaking genes")+
  xlab("ciradian time")

bulk_result<-read.delim('/lustre/home/acct-medll/medll/data/analysis/tauFisher_test/GSE113883_hs_full_adj_meta2d/meta2d_GSE113883_hs_full_adj.txt')
bulk_result_IDs<-bulk_result[bulk_result$JTK_pvalue<0.01,"CycID"]
ggplot(JTK_CYCL.filtered[JTK_CYCL.filtered$celltype!="Eryth"&(JTK_CYCL.filtered$CycID%in%bulk_result_IDs),])+
  geom_bar(aes(x=LAG,fill=celltype),alpha=0.5,color="black")+
  facet_wrap(~celltype,scales = "free_y",ncol=6)+
  scale_fill_manual(values=use_color)+NoLegend()+
  theme(axis.text=element_text(color="black"))+ylab("count of peaking genes")+
  xlab("ciradian time")

#FIG. 2D
plotPhaseofCellType(JTK_CYCL.filtered[JTK_CYCL.filtered$celltype!="Eryth",],p.cutoff = 0.01,
                    colors=generateColor(30,col.dist.min=0.3,seed = 10)[-11])+NoLegend()

#FIG. 2E
plotdataFIG2E<-NULL
celltypes<-sort(unique(JTK_CYCL.filtered$celltype))
use_color<-generateColor(30,col.dist.min=0.3,seed = 10)
for(i in 1:length(celltypes)){
  genes=JTK_CYCL.filtered[JTK_CYCL.filtered$celltype==celltypes[i],"CycID"]
  gene_ids=bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID
  out=enrichGO(gene_ids,OrgDb = org.Hs.eg.db)
  out=out@result[out@result$p.adjust<0.05,c("Count","Description","p.adjust")]
  if(nrow(out)!=0){
    out$type=celltypes[i]
  }else{
    next
  }
  if(is.null(plotdataFIG2E)){
    plotdataFIG2E=out
  }else{
    plotdataFIG2E=rbind(plotdataFIG2E,out)
  }
}
use_color<-use_color[celltypes %in% plotdataFIG2E$type]
plotdataFIG2E$Description<-getField(plotdataFIG2E$Description,sep = ",",field = 1)

ggplot(plotdataFIG2E)+geom_point(aes(x=Count,y=Description,size=-log10(p.adjust),color=type),alpha=0.75)+
  scale_color_manual(values=use_color)+
  theme(panel.grid=element_line(color="grey"),panel.background=element_rect(fill="white",color="black"))+
  xlab("count of genes")+ylab("")+guides(color=guide_legend(override.aes=list(size=3)))




#FIG. S2A
test$type<-test$manual.level1
DotPlot(test,CIRCADIAN_GENES_MAIN,group.by="type",scale = FALSE,
        cols = c("#507AAF","#BE5C37"))+scale_y_discrete(limits=unique(test$type))+
  scale_size_continuous(range = c(0,4))+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=9),axis.text.y=element_text(size=9),
        panel.border = element_rect(color = "black",linewidth = 1),
        legend.title=element_text(size=11))+xlab("")+ylab("")+coord_flip()

#FIG. S2B
counts.main<-FetchData(test,vars=c(CIRCADIAN_GENES_MAIN),layer = "count")
meta<-test@meta.data[,c("nCount_RNA","manual.level1")]
counts.main<-mergeByRownames(counts.main,meta)
counts.main$CLOCK<-counts.main$CLOCK*10000/counts.main$nCount_RNA
counts.main$ARNTL<-counts.main$ARNTL*10000/counts.main$nCount_RNA
counts.main$PER1<-counts.main$PER1*10000/counts.main$nCount_RNA
counts.main$PER2<-counts.main$PER2*10000/counts.main$nCount_RNA
counts.main$PER3<-counts.main$PER3*10000/counts.main$nCount_RNA
counts.main$CRY1<-counts.main$CRY1*10000/counts.main$nCount_RNA
counts.main$CRY2<-counts.main$CRY2*10000/counts.main$nCount_RNA
counts.main$NR1D1<-counts.main$NR1D1*10000/counts.main$nCount_RNA
counts.main$NR1D2<-counts.main$NR1D2*10000/counts.main$nCount_RNA
counts.main$CT<-getField(rownames(counts.main),"_",3)
counts.main$individual<-getField(rownames(counts.main),"_",c(1,2))
counts.main<-counts.main %>% group_by(.,manual.level1,CT,individual) %>% summarise(.,CLOCK=mean(CLOCK),
                    ARNTL=mean(ARNTL),PER1=mean(PER1),PER2=mean(PER2),
                    PER3=mean(PER3),CRY1=mean(CRY1),
                    CRY2=mean(CRY2),NR1D1=mean(NR1D1),NR1D2=mean(NR1D2))
counts.main<-counts.main %>% gather(.,"gene","mean_normalized_count",-c("manual.level1","individual","CT"))
counts.line<-counts.main %>% group_by(.,manual.level1,CT,gene) %>% summarise(.,median=median(mean_normalized_count))
ggplot(counts.main)+
  geom_boxplot(aes(x=CT,y=mean_normalized_count,color=gene))+
  geom_line(data=counts.line,aes(x=CT,y=median,group=gene,color=gene))+
  scale_x_discrete(limits=CT_TIME_ORDER)+facet_wrap(~manual.level1)+
  ylab("")+scale_color_d3()+theme_bw()+ylim(c(0,1))

#FIG. S2C
test$type<-test$predicted.celltype.l1.5
DotPlot(test,CIRCADIAN_GENES_MAIN,group.by="type",scale = FALSE,
        cols = c("#507AAF","#BE5C37"))+scale_y_discrete(limits=sort(unique(test$type)))+
  scale_size_continuous(range = c(0,4))+
  theme(axis.text.x=element_text(angle=60,hjust=1,size=9),axis.text.y=element_text(size=9),
        panel.border = element_rect(color = "black",linewidth = 1),
        legend.title=element_text(size=11))+xlab("")+ylab("")+coord_flip()
#FIG.S2D
counts.main<-FetchData(test,vars=c(CIRCADIAN_GENES_MAIN),layer = "count")
meta<-test@meta.data[,c("nCount_RNA","type")]
counts.main<-mergeByRownames(counts.main,meta)
counts.main$CLOCK<-counts.main$CLOCK*10000/counts.main$nCount_RNA
counts.main$ARNTL<-counts.main$ARNTL*10000/counts.main$nCount_RNA
counts.main$PER1<-counts.main$PER1*10000/counts.main$nCount_RNA
counts.main$PER2<-counts.main$PER2*10000/counts.main$nCount_RNA
counts.main$PER3<-counts.main$PER3*10000/counts.main$nCount_RNA
counts.main$CRY1<-counts.main$CRY1*10000/counts.main$nCount_RNA
counts.main$CRY2<-counts.main$CRY2*10000/counts.main$nCount_RNA
counts.main$NR1D1<-counts.main$NR1D1*10000/counts.main$nCount_RNA
counts.main$NR1D2<-counts.main$NR1D2*10000/counts.main$nCount_RNA
counts.main$CT<-getField(rownames(counts.main),"_",3)
counts.main$individual<-getField(rownames(counts.main),"_",c(1,2))
counts.main<-counts.main %>% group_by(.,type,CT,individual) %>% summarise(.,CLOCK=mean(CLOCK),
                                                                                   ARNTL=mean(ARNTL),PER1=mean(PER1),PER2=mean(PER2),
                                                                                   PER3=mean(PER3),CRY1=mean(CRY1),
                                                                                   CRY2=mean(CRY2),NR1D1=mean(NR1D1),NR1D2=mean(NR1D2))
counts.main<-counts.main %>% gather(.,"gene","mean_normalized_count",-c("type","individual","CT"))
counts.line<-counts.main %>% group_by(.,type,CT,gene) %>% summarise(.,median=median(mean_normalized_count))
top<-ggplot(counts.main[counts.main$type=="MAIT",])+
  geom_boxplot(aes(x=CT,y=mean_normalized_count,color=gene))+
  geom_line(data=counts.line[counts.line$type=="MAIT",],aes(x=CT,y=median,group=gene,color=gene))+
  scale_x_discrete(limits=CT_TIME_ORDER)+facet_wrap(~type)+xlab("")+
  ylab("")+scale_color_d3()+theme_bw()+ylim(c(0,1))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

bottom<-ggplot(counts.main[counts.main$type=="NK",])+
  geom_boxplot(aes(x=CT,y=mean_normalized_count,color=gene))+
  geom_line(data=counts.line[counts.line$type=="NK",],aes(x=CT,y=median,group=gene,color=gene))+
  scale_x_discrete(limits=CT_TIME_ORDER)+facet_wrap(~type)+
  ylab("")+scale_color_d3()+theme_bw()+ylim(c(0,1))

ggarrange(top,bottom,nrow = 2,common.legend = T,legend ="right")


#FIG. S2E
celltypes<-names(table(test$predicted.celltype.l2))
JTK.result<-getJTK_CYCLEouts("../analysis/OSgenePredict_PseudoBulkKnownMM_all_data/",celltypes = celltypes,padj = 1)
JTK.result.filtered<-getJTK_CYCLEouts("../analysis/OSgenePredict_PseudoBulkKnownMM_all_data/",celltypes = celltypes,padj = 0.001)
validcelltypes=JTK.result.filtered$celltype %>% unique()
top<-ggplot(JTK.result.filtered)+geom_bar(aes(x=celltype),fill="black",color="black")+
  scale_x_discrete(limits=validcelltypes,labels=getField(validcelltypes,sep = "_",1))+
  scale_y_continuous(expand=c(0,0),limits = c(0,50))+xlab("")+NoLegend()+ylab("count of \noscillating genes")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.background=element_rect(fill="white"),
        plot.margin = margin(0.5,0,0,0,"cm"))

share.data=NULL
for(type.this in validcelltypes){
  for(type.that in validcelltypes){
    n.gene.type.this=nrow(JTK.result.filtered[JTK.result.filtered$celltype==type.this,])
    n.share.gene=length(intersect(JTK.result.filtered[JTK.result.filtered$celltype==type.this,"CycID"],
                                  JTK.result.filtered[JTK.result.filtered$celltype==type.that,"CycID"]))
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
bottom<-ggplot(share.data)+geom_raster(aes(x=from.type,y=to.type,fill=shared.ratio))+
  scale_x_discrete(labels=getField(validcelltypes,sep = "_",1))+
  scale_fill_gradientn(colours = c("#283168","#009fc0","#f0f0cc","white"),values = c(0,0.5,0.9,1))+
  scale_y_discrete(limits=rev(validcelltypes),labels=rev(getField(validcelltypes,sep = "_",1)))+
  theme(axis.text.x = element_text(angle=60,hjust=1),plot.margin=margin(0,0,0,0,"cm"))+ylab("")+xlab("")

ggarrange(top,bottom,nrow=2,align = "v",heights = c(1,3),common.legend = T,legend = "right")

#FIG. S2F

MM_JTK<-getJTK_CYCLEouts('../analysis/OSgenePredict_PseudoBulkKnownMM/',celltypes = celltypes,filter.data = F)
MM_JTK.filtered<-getJTK_CYCLEouts('../analysis/OSgenePredict_PseudoBulkKnownMM/',celltypes = celltypes)
plotdata<-data.frame(group="test",ratio=length(MM_JTK.filtered$CycID %>% unique())/length(MM_JTK$CycID %>% unique()))

for(i in 1:6){
  path.char=paste0("../analysis/OSgenePredict_PseudoBulk352_Random_test/rep",i,"/")
  Ctrl_JTK<-getJTK_CYCLEouts(path.char,celltypes = celltypes,filter.data = F)
  Ctrl_JTK.filtered<-getJTK_CYCLEouts(path.char,celltypes = celltypes)
  plotdata=rbind(plotdata,data.frame(group="control",ratio=length(Ctrl_JTK.filtered$CycID %>% unique())/length(Ctrl_JTK$CycID %>% unique())))                               
}

ggplot(plotdata)+geom_boxplot(aes(x=group,y=ratio))+
  #stat_compare_means(aes(x=group,y=ratio))+
  theme_classic2()

#FIG. 3A
#compare bulk and scPseudoBulk oscillating genes
bulk_result<-read.delim('/lustre/home/acct-medll/medll/data/analysis/tauFisher_test/GSE113883_hs_full_adj_meta2d/meta2d_GSE113883_hs_full_adj.txt')
bulk_result_IDs<-bulk_result[bulk_result$JTK_pvalue<0.01,"CycID"]
JTK_CYCL<-NULL
for(block_data in list.files("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5_health/")){
  block_path=paste0("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5_health/",block_data,"/")
  block=getJTK_CYCLEouts(block_path,celltypes=celltypes,filter.data = F)
  if(is.null(JTK_CYCL)){
    JTK_CYCL=block
  }else{
    JTK_CYCL=rbind(JTK_CYCL,block)
  }
}
pseudo_bulk_IDs<-dplyr::filter(JTK_CYCL,ADJ.P<0.01)$CycID %>% unique()
JTK_CYCL_filter<-dplyr::filter(JTK_CYCL,ADJ.P<0.01)
shared<-intersect(pseudo_bulk_IDs,bulk_result_IDs)
vennplot<-VennDiagram::venn.diagram(list("SC"=pseudo_bulk_IDs,"bulk"=bulk_result_IDs),filename = NULL)
grid.draw(vennplot)
dev.off()
setdiff(bulk_result_IDs,shared)
LayerData(test,"count")[intersect(setdiff(bulk_result_IDs,shared),Features(test)),] %>% rowSums() %>% median()
LayerData(test,"count")[shared,] %>% rowSums() %>% median()
JTK_CYCL_filter[JTK_CYCL_filter$CycID%in%shared,"CycID"] %>% table() %>% density() %>% plot()

#these genes identified by polygon recognition did not show obvious pattern
all_table=NULL
for(file in list.files('../analysis/velosity_rhythm_highly_expressed_genes_3000_by_type/') %>% grep("^c",.,value = T)){
  file_table=readRDS(paste0('../analysis/velosity_rhythm_highly_expressed_genes_3000_by_type/',file))
  if(is.null(all_table)){
    all_table=file_table
  }else{
    all_table=rbind(all_table,file_table)
  }
}
all_table_join<-all_table %>% group_by(.,feature) %>% summarise(.,n.observation=n())
all_table_bytype<-all_table[c("type","feature")] %>% unique() %>% group_by(.,feature) %>% summarise(.,n.celltype=n())
all_table_byindividual<-all_table[c("individual","feature")] %>% unique() %>% group_by(.,feature) %>% summarise(.,n.individual=n())
all_table_join<-left_join(all_table_join,all_table_bytype,by="feature")
all_table_join<-left_join(all_table_join,all_table_byindividual,by="feature")
all_table_join$mean.observation<-all_table_join$n.observation/all_table_join$n.celltype
all_table_join.filtered<-all_table_join[all_table_join$n.individual>=6,]
#intersect(all_table_join.filtered$feature, JTK.result.filtered$CycID)
cell.meta<-test@meta.data
all_table_join.filtered %>% arrange(.,desc(n.observation))
for(gene in all_table_join.filtered$feature){
  n.celltype=all_table[all_table$feature==gene,"type"] %>% unique() %>% length()
  plotSplicingCountsByGene(all_splice,cell.meta,all_table[all_table$feature==gene,],normalize = T)
  ggsave(paste0("../analysis/velosity_rhythm_highly_expressed_genes_3000_by_type_norm/pdf/",gene,".pdf"),
         width=12,height=n.celltype*2)
}
plotSplicingCountsByGene(all_splice,cell.meta,all_table[all_table$feature=="STK38",],normalize = T)
plotSplicingCountsByGene(all_splice,cell.meta,all_table[all_table$feature=="RPLP0",],normalize = F)
all_table[all_table$feature=="UCP2",]
interest_genes<-all_table_join.filtered$feature

#still using JTK cycle, but using pseudobulk
test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.Azimuth.rds")
#test2<-readRDS("../analysis/36.samples.oddsremoved.tsneran.rds")
#celltypes<-test$manual.level2 %>% unique() %>% sort()
celltypes<-test$predicted.celltype.l2 %>% unique() %>% sort()
#celltypes2<-test2$manual.level2 %>% unique() %>% sort()
#test$type<-test$manual.level2
test$type<-test$predicted.celltype.l2
#test2$type<-test2$manual.level2
genes_wanted<-intersect(Features(test),circadian_gene_hs$gene_name)

#genes_wanted<-scan('../list/circadian_genes_human.lst',what = "character")
#genes_wanted<-intersect(Features(test),genes_wanted)
genes_wanted<-setdiff(Features(test),genes_wanted) %>% sample(.,356)
intersect(genes_wanted,circadian_gene_hs$gene_name)
for(i in 1:length(celltypes)){
  usecelltype=celltypes[i]
  if(file.exists(paste0('../analysis/OSgenePredict_PseudoBulkKnownMM_all_data_random356/rep3/',usecelltype,".txt"))){
    next
  }
  message(usecelltype)
  detectCircadianGenesByPseudoBulk(test,features=genes_wanted,
      usecelltype = usecelltype,out.dir = '../analysis/OSgenePredict_PseudoBulkKnownMM_all_data_random356/rep3/',
      rep=3,proportion=0.3)
}

testJTK<-getJTK_CYCLEouts('../analysis/OSgenePredict_PseudoBulkKnownMM_all_data/',celltypes,padj = 0.01)

#diffusionmap
test$type<-test$predicted.celltype.l1.5
testCells<-subset(test,type=="CD14 Mono")
testCells<-NormalizeData(testCells)
testCells<-FindVariableFeatures(testCells,nfeatures = 4000)
testCellCounts<-FetchData(testCells,vars=VariableFeatures(testCells))
housekeepers<-testCellCounts %>% colSums() %>% sort() %>% tail(.,50) %>% names()
testCellCounts<-testCellCounts %>% rownames_to_column(.,var="Cell")
#testCellCounts$CT<-as.numeric(testCells$CT %>% gsub("CT","",.))
testCellCounts<-as.ExpressionSet(as.data.frame(testCellCounts))
normalizations <- colMeans(exprs(testCellCounts)[housekeepers,])
guo_norm <- testCellCounts
exprs(guo_norm) <- exprs(guo_norm) - normalizations
dm <- DiffusionMap(guo_norm,n_pcs = 50)
plot(dm,1:2,pch = 20,legend_main = 'circadian time')
DCs<-dm@eigenvectors %>% as.data.frame()
DCs$CT<-testCells$CT %>% as.character()
ggplot(DCs)+geom_point(aes(x=DC1,y=DC2,color=CT))

#plot count of count for core circadian genes
countdata<-FetchData(test,vars=CIRCADIAN_GENES_MAIN,layer = "count")
countdata_n<-gather(countdata,key="gene",value="count")
countdata_n<-countdata_n %>% table() %>% as.data.frame()
countdata_n$total<-Cells(test) %>% length()
countdata_n$ratio<-round(countdata_n$Freq*100/countdata_n$total,digits = 2)
ggplot(countdata_n)+geom_bar(aes(x=count,y=log10(Freq+1)),stat = "identity")+
  geom_text(data=countdata_n[countdata_n$ratio>=0.01,],aes(x=count,y=log10(Freq+1)/2,label=ratio),color="white")+
  facet_wrap(~gene,nrow=3)+coord_flip()+xlab("UMI count")+ylab("log10(1+count of cell)")

countdata<-FetchData(subset(test,type=="MAIT"),vars=CIRCADIAN_GENES_MAIN,layer = "count")
countdata_n<-gather(countdata,key="gene",value="count")
countdata_n<-countdata_n %>% table() %>% as.data.frame()
countdata_n$total<-Cells(subset(test,type=="MAIT")) %>% length()
countdata_n$ratio<-round(countdata_n$Freq*100/countdata_n$total,digits = 2)
ggplot(countdata_n)+geom_bar(aes(x=count,y=log10(Freq+1)),stat = "identity")+
  geom_text(data=countdata_n[countdata_n$ratio>=0.01&log10(countdata_n$Freq+1)>=1.5,],aes(x=count,y=log10(Freq+1)/2,label=ratio),color="white")+
  geom_text(data=countdata_n[countdata_n$ratio>=0.01&log10(countdata_n$Freq+1)<1.5,],aes(x=count,y=log10(Freq+1)+1,label=ratio),color="black")+
  facet_wrap(~gene,nrow=3)+coord_flip()+xlab("UMI count")+ylab("log10(1+count of cell)")+ggtitle("MAIT")

plotPseudobulkByCT(test,cell.type="CD16 Mono",features=c("ARNTL","DBP","PER1","PER2","PER3","NR1D1","NR1D2","CRY1","CRY2"))

test$immune_cell<-ifelse(test$predicted.celltype.l1.5 %in% c("ASDC","B memory","B naive","CD14 Mono",
                              "CD16 Mono","CD4 T","CD8 T","cDC1","cDC2","dnT"),"immune","non-immune")
test$type<-test$predicted.celltype.l1.5
test$type<-test$immune_cell
plotPseudobulkByCT(test,cell.type="immune",
                   features=c("ARNTL","DBP","PER1","PER2","PER3","NR1D1","NR1D2","CRY1","CRY2"),
                   rep = 5,proportion = 0.2)
plotPseudobulkByCT(test,cell.type="CD8 T",
                   features=c("ARNTL","DBP","PER1","PER2","PER3","NR1D1","NR1D2","CRY1","CRY2"),
                   rep = 5,proportion = 0.2)

#plot cells ratio by CT
cells_CT<-test@meta.data[,c("CT","type","patient")]
cells_CT<-cells_CT %>% group_by(.,CT,type,patient) %>% summarise(.,count=n())
cells_CT_total<-cells_CT %>% group_by(CT,patient) %>% summarise(.,total=sum(count))
cells_CT<-left_join(cells_CT,cells_CT_total,by=c("CT","patient"))
cells_CT$ratio<-cells_CT$count/cells_CT$total
cells_CT$CT_numeric<-cells_CT$CT %>% gsub("CT","",.) %>% as.numeric()
cells_CT_clean<-cells_CT[cells_CT$type!="Eryth",]
cells_CT_median<-cells_CT_clean %>% group_by(.,CT_numeric,type) %>% summarise(.,median=median(ratio))
ggplot(cells_CT_clean)+geom_boxplot(aes(x=CT_numeric,y=ratio,group=CT_numeric))+
  geom_line(data=cells_CT_median,aes(x=CT_numeric,y=median))+
  facet_wrap(~type,nrow = 3,scales = "free_y")

cells_CT<-test@meta.data[,c("CT","type")]
cells_CT<-cells_CT %>% group_by(.,CT,type) %>% summarise(.,count=n())
cells_CT_total<-cells_CT %>% group_by(CT) %>% summarise(.,total=sum(count))
cells_CT<-left_join(cells_CT,cells_CT_total,by=c("CT"))
cells_CT$ratio<-cells_CT$count/cells_CT$total
cells_CT$CT_numeric<-cells_CT$CT %>% gsub("CT","",.) %>% as.numeric()
cells_CT_clean<-cells_CT[cells_CT$type!="Eryth",]

ggplot(cells_CT_clean)+geom_boxplot(aes(x=CT_numeric,y=ratio,group=CT_numeric))+
  geom_line(aes(x=CT_numeric,y=ratio))+
  facet_wrap(~type,nrow = 3,scales = "free_y")


test<-readRDS('../analysis/all.samples.oddsremoved.tsneran.Azimuth.healthonly.rds')
test$type<-test$predicted.celltype.l1.5
plotPseudobulkByCT(test,"CD14 Mono",features = c("IL10","IL6","TNF"),proportion = 0.25)
plotPseudobulkByCT(test,"CD14 Mono",features=CIRCADIAN_GENES_MAIN,proportion = 0.25)
plotPseudobulkByCTByIndividual(test,"CD14 Mono",features = c("PER1"),proportion = 1)

test<-readRDS('../analysis/all.samples.oddsremoved.tsneran.Azimuth.healthonly.SCTransform.rds')


#using tauFisher to predict cell internal time
#get circadian genes for one cell type and filter for top10 significant genes
JTK_CYCL<-NULL
for(block_data in list.files("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5/")){
  block_path=paste0("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5/",block_data,"/")
  block=getJTK_CYCLEouts(block_path,celltypes=celltypes,filter.data = T,padj = 0.01)
  if(is.null(JTK_CYCL)){
    JTK_CYCL=block
  }else{
    JTK_CYCL=rbind(JTK_CYCL,block)
  }
}

JTK_CYCL_bulk<-read.delim('/lustre/home/acct-medll/medll/data/analysis/tauFisher_test/GSE113883_hs_full_adj_meta2d/meta2d_GSE113883_hs_full_adj.txt')
JTK_CYCL_bulk<-dplyr::filter(JTK_CYCL_bulk,JTK_pvalue<0.01)

celltype=celltypes[7]
#strong_oscillating<-JTK_CYCL[JTK_CYCL$celltype==celltype&JTK_CYCL$CycID %in% JTK_CYCL_bulk$CycID,] %>% arrange(.,ADJ.P) %>% head(.,20)
strong_oscillating<-JTK_CYCL[JTK_CYCL$celltype==celltype,] %>% arrange(.,ADJ.P) %>% head(.,20)
strong_oscillating<-strong_oscillating$CycID
cc_genes<-toupper(c("Arntl", "Dbp", "Nr1d1", "Nr1d2", 
                     "Per1", "Per2", "Per3", "Cry1", "Cry2"))
strong_oscillating<-unique(c(strong_oscillating,cc_genes))
#plotPseudobulkByCT(test,"CD14 Mono",features = strong_oscillating,rep = 3,proportion = 0.1)
#get Pseudobulk data, relative to the bulk input in taufisher

Pseudobulk0<-plotPseudobulkByCT(test,celltype,features = strong_oscillating,return.data = T,rep = 3,proportion = 0.2)
Pseudobulk0$time_numeric<-Pseudobulk0$time %>% gsub("CT","",.) %>% as.numeric()
Pseudobulk0$time_numeric<-Pseudobulk0$time_numeric+24*(Pseudobulk0$downsample-1)
Pseudobulk0$downsample<-NULL
Pseudobulk0$time<-NULL

Pseudobulk<-spread(Pseudobulk0,key=feature,value=count) %>% as.data.frame()
rownames(Pseudobulk)<-Pseudobulk$time_numeric %>% as.character()
Pseudobulk$time_numeric<-NULL

bulk_log<-log2(Pseudobulk+1)
time_adj<-rownames(Pseudobulk) %>% as.numeric()
#run FDA
nrep = 3 # number of replicates
fda_expression = get_FDAcurves(dat=bulk_log, 
                               time=time_adj, 
                               numbasis=5) %>%
  dplyr::mutate(time_24 = fda_time - 24*floor(fda_time/24)) %>%
  dplyr::mutate(time_label = paste0(time_24, "_",
                                    rep((max(nrep)+1):(max(nrep)+max(nrep)),
                                        each = 24, length = nrow(.))))
new_fda_rownames = fda_expression$time_label

# Remove the unnecessary columns
fda_expression2 = fda_expression[, -c(1, ncol(fda_expression)-1, ncol(fda_expression))]
# Create the differences matrix and scale
fda_diff <- create_DiffMatrix(genes=strong_oscillating, dat=fda_expression2)
fda_diff_scaled = scale_DiffMatrix(diffs=fda_diff)

fda_mat = as.matrix(fda_diff_scaled)
#run PCA
# Set up train data for PCA
train=fda_mat
rownames(train)=new_fda_rownames
train_time=fda_expression$time_24

# PCA
X_PCA<-train
X_PCA<-as.data.frame(X_PCA)
X_PCA$CT<-train_time

pc<-stats::prcomp(X_PCA[, -ncol(X_PCA)], scale = FALSE)
ggplot(pc$x)+geom_point(aes(x=PC1,y=PC2))
#Run multinomial regression
ndims = 2
pc_data<-data.frame(pc$x[,1:ndims])

# Get the times and relevel them so the smallest CT is the reference level
pc_data$CT24<-as.numeric(vapply(stringr::str_split(rownames(pc_data), pattern = '_'),
                                '[', 1, FUN.VALUE = character(1) ))
pc_data$CT24_relevel<-stats::relevel(factor(pc_data$CT24),
                                     ref = as.character(min(train_time)))
mod <- nnet::multinom(CT24_relevel ~ PC1 + PC2, data = pc_data, trace=T)

#get single cell data
#test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.Azimuth.healthonly.rds")
all_predicted=NULL
for(i in 1:length(HEALTH)){
  sc_data<-plotPseudobulkByCTByIndividual(test,celltype,individuals = HEALTH[i],features = strong_oscillating,return.data = T,proportion = 1,layer = "count")
  sc_data<-sc_data[,c("time","count","feature")]
  sc_data<-sc_data %>% spread(key=feature,value=count) %>% as.data.frame()
  rownames(sc_data)<-sc_data$time
  sc_data$time<-NULL
  #good_cell<-colSums(sc_data) %>% sort() %>% tail(.,100) %>% names()
  #sc_data<-t(sc_data[,good_cell])
  sc_data_log<-log2(sc_data+1)

  #Calculate differences for each gene pair
  sc_data_diff <- create_DiffMatrix(genes=strong_oscillating, dat=sc_data_log)
  sc_data_diff_scaled <- scale_DiffMatrix(diffs=sc_data_diff)

  # Project data onto PCA space
  pc_pred <- stats::predict(pc, newdata = sc_data_diff_scaled)
  # Predict
  pred_vals = stats::predict(mod, newdata = pc_pred[,1:ndims])
  pred_vals = as.numeric(as.character(pred_vals))

  density(pred_vals) %>% plot()

  truth_time<-rownames(sc_data_diff_scaled) %>% gsub("CT","",.) %>% as.numeric()
  out=calc_error(truth=truth_time, pred=pred_vals)$differences
  out$individual=HEALTH[i]
  file_name=paste0("/lustre/home/acct-medll/medll/data/analysis/tauFhisher_predict/",HEALTH[i],"_",gsub(" ","_",celltype),".tsv")
  write_delim(out,file_name,delim = "\t",col_names = T)
  if(is.null(all_predicted)){
    all_predicted=out
  }else{
    all_predicted=rbind(all_predicted,out)
  }
}

ggplot(all_predicted)+geom_boxplot(aes(x=individual,y=differences,group=individual))+theme_publication()+ggtitle(celltype)


JTK_CYCL<-NULL
for(block_data in list.files("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5/")){
  block_path=paste0("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5/",block_data,"/")
  block=getJTK_CYCLEouts(block_path,celltypes=celltypes,filter.data = F)
  if(is.null(JTK_CYCL)){
    JTK_CYCL=block
  }else{
    JTK_CYCL=rbind(JTK_CYCL,block)
  }
}
JTK_CYCL<-JTK_CYCL[JTK_CYCL$ADJ.P<0.01,]

JTK_CYCL_bulk<-read.delim('/lustre/home/acct-medll/medll/data/analysis/tauFisher_test/GSE113883_hs_full_adj_meta2d/meta2d_GSE113883_hs_full_adj.txt')
JTK_CYCL_bulk<-dplyr::filter(JTK_CYCL_bulk,JTK_pvalue<0.01)
bulkshared<-LayerData(test,features=intersect(JTK_CYCL_bulk$CycID,unique(JTK_CYCL$CycID)),layer = "count") %>% rowSums()
bulkspec<-LayerData(test,features=setdiff(JTK_CYCL_bulk$CycID,unique(JTK_CYCL$CycID)),layer = "count") %>% rowSums
bulkshared<-as.data.frame(bulkshared)
bulkspec<-as.data.frame(bulkspec)
colnames(bulkshared)<-"count"
colnames(bulkspec)<-"count"
bulkshared$category<-"shared with scRNA-seq"
bulkspec$category<-"bulk specific"
genes_compare<-rbind(bulkshared,bulkspec)
genes_compare$log2<-log2(genes_compare$count+1)
ggplot(genes_compare)+geom_boxplot(aes(x=category,y=log2,fill=category),outliers = F)+
  stat_compare_means(aes(x=category,y=log2,group=category),vjust = 4)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.background=element_rect(fill = "white"),
        axis.line = element_line())+
  xlab("")+ylab("log2(UMI count)")

#cortisol data analysis
cortisol_data_test<-data.frame(concentrate=c(400,200,100,50,25,12.5,6.25,0),OD=c(0.089,0.145,0.236,0.440,0.687,1.098,1.414,1.910))
ggplot(cortisol_data_test)+geom_point(aes(x=log2(concentrate),y=log10(OD)))+
  theme(panel.background=element_rect(fill = "white"),
        axis.line = element_line())
lm.cortisol<-lm(log2(concentrate)~log10(OD),cortisol_data_test[-8,])
new<-data.frame("OD"=c(0.536,0.562,0.607,0.573))
data<-2^predict(lm.cortisol,new) %>% as.data.frame()
colnames(data)<-"concentrate"
data$individual<-c(1,2,1,2)
data$CT<-c("CT10","CT10","CT16","CT16")
ggplot(data)+geom_boxplot(aes(x=CT,y=concentrate,group=CT))+
  stat_compare_means(aes(x=CT,y=concentrate,group=CT),method = "t")+
  ylim(c(25,35))+
  theme(panel.background=element_rect(fill = "white"),
        axis.line = element_line())+
  xlab("")+ylab("concentrate (ng/ml)")

#vdj analysis
#read BCR data
test<-readRDS('../analysis/all.samples.oddsremoved.tsneran.Azimuth.rds')
test$type<-test$predicted.celltype.l1.5
BCR_data<-NULL
for(patient in PATIENTS){
  for(i in 1:6){
    file=file.path(paste0(patient,"_",i),"outs/per_sample_outs",paste0(patient,"_",i),"vdj_b/vdj_results/filtered_contig_igblast_db-pass.tsv")
    if(file.exists(file)){
      BCR_sample_data=readChangeoDb(file)
      BCR_sample_data$individual=patient
      BCR_sample_data$CT=CT_TIME[i]
      BCR_sample_data$cell_id=paste0(patient,"_",CT_TIME[i],"_",BCR_sample_data$cell_id)
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
dist_nearest<-distToNearest(dplyr::filter(BCR_prim,locus=="IGH"),nproc=1)

# find threshold for cloning automatically
threshold_output <- findThreshold(dist_nearest$dist_nearest,
                                  method = "gmm", model = "gamma-norm",
                                 cutoff = "user", spc = 0.995,progress = T)
#threshold=0.1218473
threshold<-threshold_output@threshold

# call clones using hierarchicalClones
results <- hierarchicalClones(BCR_prim, cell_id = 'cell_id',
                              threshold = threshold, only_heavy = FALSE,
                              split_light = TRUE, summarize_clones = FALSE)
results$sample<-getField(results$cell_id,"_",1:3)
results$sequence_id<-paste(results$sequence_id,results$sample,sep="_")
celltype.meta<-test$type %>% as.data.frame()
colnames(celltype.meta)[1]<-"cell_type"
celltype.meta<-rownames_to_column(celltype.meta,var="cell_id")
results<-left_join(results,celltype.meta)
results<-dplyr::filter(results,cell_type %in% c("B memory","B naive","Plasma"))
#saveRDS(results,"../analysis/BCR_filtered.clone_type.rds")

results_clone_num<-results[c("individual","CT","clone_id","cell_type","cell_id")] %>% unique()
results_clone_num<-balanceObsevations(results_clone_num,c("individual","cell_type"),"CT")
#one clone type seldom have replicates, compare diversity across time series is of less significance
results_clone_num$clone_id %>% unique() %>% length()
results_clone_num<-results_clone_num %>% group_by(.,individual,CT,cell_type) %>% summarise(.,clone_num=n())
results_clone_num<-results_clone_num %>% arrange(.,individual,cell_type)
#can not scale because clone number is same across time
results_clone_num$scaled_clone_num<-(results_clone_num %>% group_by(.,individual,cell_type) %>% reframe(.,scaled_clone_num=scale(clone_num)))$scaled_clone_num[,1]
ggplot(results_clone_num)+geom_boxplot(aes(x=CT,y=scaled_clone_number,group=CT))+
  scale_x_discrete(limits=CT_TIME_ORDER)+facet_wrap(~cell_type)


out<-(results[,c("clone_id","sample")] %>% unique())$clone_id %>% table() %>% sort()
shared_clones<-out[out>1] %>% names()
plotdata<-results[,c("sample","clone_id")] %>% unique() %>% dplyr::filter(.,clone_id %in% shared_clones)
ggplot(plotdata)+geom_point(aes(x=clone_id,y=sample))+theme(axis.text.x = element_text(angle=45,hjust=1))

plotdata2<-data.frame(clone_id=character(0),sample1=character(0),sample2=character(0))
for(clone in plotdata$clone_id %>% unique){
  samples=plotdata[plotdata$clone_id==clone,]$sample
  if(plotdata[plotdata$clone_id==clone,] %>% nrow()==2){
    plotdata2=rbind(plotdata2,data.frame(clone_id=clone,sample1=samples[1],sample2=samples[2]))
  }else{
    plotdata2=rbind(plotdata2,data.frame(clone_id=clone,sample1=samples[c(1,1)],sample2=samples[c(2,3)]))
  }
}
plotdata$individual<-getField(plotdata$sample,sep = "_",field = 1:2)
plotdata$CT<-getField(plotdata$sample,sep = "_",field = 3)
plotdata2$individual1<-getField(plotdata2$sample1,sep = "_",field = 1:2)
plotdata2$individual2<-getField(plotdata2$sample2,sep = "_",field = 1:2)
plotdata2$CT1<-getField(plotdata2$sample1,sep = "_",field = 3)
plotdata2$CT2<-getField(plotdata2$sample2,sep = "_",field = 3)
patient_id<-11
ggplot(plotdata[plotdata$individual==PATIENTS[patient_id],])+geom_point(aes(x=clone_id,y=CT))+theme(axis.text.x = element_text(angle=45,hjust=1))+
  geom_segment(data=plotdata2[plotdata2$individual1==PATIENTS[patient_id]&plotdata2$individual2==PATIENTS[patient_id],],aes(x=clone_id,xend=clone_id,y=CT1,yend=CT2))+
  scale_y_discrete(limits=rev(CT_TIME_ORDER))

ggplot(BCR_prim[grepl("IGH",BCR_prim$locus),])+geom_bar(aes(x=CT,fill=c_call),position="fill")+
  scale_x_discrete(limits=CT_TIME_ORDER)+scale_fill_aaas()

results<-readRDS('../analysis/BCR_filtered.clone_type.rds')
results$CT<-factor(results$CT,levels = CT_TIME_ORDER)
plotdata<-results[,c("CT","clone_id","cell_type")] %>% unique()
plotdata<-plotdata[plotdata$cell_type=="B memory",]
plotdata$CT<-as.vector(plotdata$CT)
shared_clone_CR<-(plotdata$clone_id %>% table)[(plotdata$clone_id %>% table)==2] %>% names()
circos.par("track.height" = 0.1)
circos.initialize(results$CT, x = results$clone_id)
circos.track(results$CT, y = results$v_score,bg.col=c("grey","white","white","white","grey","grey"),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           (CELL_META$cell.ylim[2] + CELL_META$cell.ylim[1])/2, 
                           CELL_META$sector.index)
             })
#circos.trackHist(results.locus.data$CT, results.locus.data$Freq,ylim = )
line.index=0
for(id in shared_clone_CR){
  subdata=plotdata[plotdata$clone_id==id,]
  for(i in 1:(nrow(subdata)-1)){
    for(j in (i+1):nrow(subdata)){
      data1=subdata[i,]
      data2=subdata[j,]
      diff=abs(which(CT_TIME_ORDER==data1$CT[1])-which(CT_TIME_ORDER==data2$CT[1]))
      if(!(diff %in% c(3))){
        next
      }
      line.color=case_when(diff==1~"red",
                           diff==2~"green",
                           diff==4~"green",
                           diff==5~"red",
                           .default = "blue")
      circos.link(data1$CT[1], as.numeric(data1$clone_id[1]),
                  data2$CT[1], as.numeric(data2$clone_id[1]),
                  h = 0.3,col = line.color)
      line.index=line.index+1
      message(line.index)
    }
  }
}

circos.clear()



#read TCR data
test<-readRDS('../analysis/all.samples.oddsremoved.tsneran.Azimuth.rds')
test$type<-test$predicted.celltype.l1.5
TCR_data<-NULL
for(patient in PATIENTS){
  for(i in 1:6){
    file=file.path(paste0(patient,"_",i),"outs/per_sample_outs",paste0(patient,"_",i),"vdj_t/vdj_results/filtered_contig_igblast_db-pass.tsv")
    if(file.exists(file)){
      TCR_sample_data=readChangeoDb(file)
      TCR_sample_data$individual=patient
      TCR_sample_data$CT=CT_TIME[i]
      TCR_sample_data$cell_id=paste0(patient,"_",CT_TIME[i],"_",TCR_sample_data$cell_id)
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


TCR_prim<-TCR_data %>% dplyr::filter(productive)
# remove cells with multiple alpha/beta chain
multi_alpha<-table(dplyr::filter(TCR_data, locus == "TRA")$cell_id)
multi_beta<-table(dplyr::filter(TCR_data, locus == "TRB")$cell_id)
multi_chain_cells<-c(names(multi_alpha)[multi_alpha > 1],names(multi_beta)[multi_beta > 1]) %>% unique()
TCR_prim <- dplyr::filter(TCR_prim,!cell_id %in% multi_chain_cells)

# split cells by heavy and light chains
alpha_cells<-dplyr::filter(TCR_prim, locus == "TRA")$cell_id
beta_cells<-dplyr::filter(TCR_prim, locus == "TRB")$cell_id
no_pair_cells<-alpha_cells[which(!alpha_cells %in% beta_cells)]

TCR_prim<-dplyr::filter(TCR_prim, !cell_id %in% no_pair_cells)
dist_nearest<-distToNearest(TCR_prim,nproc=1,locusValues = c("TRB","TRA"))

# find threshold for cloning automatically
threshold_output <- findThreshold(dist_nearest$dist_nearest,
                                  method = "gmm", model = "gamma-norm",
                                  cutoff = "user", spc = 0.995,progress = T)
#threshold=0.05070953
threshold<-threshold_output@threshold

# call clones using hierarchicalClones
results <- hierarchicalClones(TCR_prim, cell_id = 'cell_id',
                              threshold = threshold, only_heavy = FALSE,
                              split_light = TRUE, summarize_clones = FALSE)
results$sample<-getField(results$cell_id,"_",1:3)
results$sequence_id<-paste(results$sequence_id,results$sample,sep="_")
celltype.meta<-test$type %>% as.data.frame()
colnames(celltype.meta)[1]<-"cell_type"
celltype.meta<-rownames_to_column(celltype.meta,var="cell_id")
results<-left_join(results,celltype.meta)
results<-dplyr::filter(results,cell_type %in% c("CD4 T","CD8 T","dnT","MAIT"))
saveRDS(results,"../analysis/TCR_filtered.clone_type.rds")

results<-readRDS("../analysis/TCR_filtered.clone_type.rds")
results$CT<-factor(results$CT,levels = CT_TIME_ORDER)

plotdata<-results[,c("CT","clone_id","individual","cell_type")] %>% unique()
plotdata<-plotdata[plotdata$individual==PATIENTS[1] & plotdata$cell_type=="CD4 T",]
plotdata$CT<-as.vector(plotdata$CT)
shared_clone_TCR<-(plotdata$clone_id %>% table)[(plotdata$clone_id %>% table)==2] %>% names()
circos.par("track.height" = 0.1)
circos.initialize(results$CT, x = results$clone_id)
circos.track(results$CT, y = results$v_score,bg.col=c("grey","white","white","white","grey","grey"),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           (CELL_META$cell.ylim[2] + CELL_META$cell.ylim[1])/2, 
                           CELL_META$sector.index)
             })
line.index=0
for(id in shared_clone_TCR){
  subdata=plotdata[plotdata$clone_id==id,]
  for(i in 1:(nrow(subdata)-1)){
    for(j in (i+1):nrow(subdata)){
      data1=subdata[i,]
      data2=subdata[j,]
      diff=abs(which(CT_TIME_ORDER==data1$CT[1])-which(CT_TIME_ORDER==data2$CT[1]))
      if(!(diff %in% c(1,5))){
        next
      }
      line.color=case_when(diff==1~"red",
                           diff==2~"green",
                           diff==4~"green",
                           diff==5~"red",
                           .default = "blue")
      circos.link(data1$CT[1], as.numeric(data1$clone_id[1]),
                  data2$CT[1], as.numeric(data2$clone_id[1]),
                  h = 0.3,col = line.color)
      line.index=line.index+1
      message(line.index)
    }
  }
}


circos.clear()

plotSegmentUsageByCT(results,"TRA","v_call")
plotSegmentUsageByCT(results,"TRA","j_call")
plotSegmentUsageByCT(results,"TRB","v_call")
plotSegmentUsageByCT(results,"TRB","d_call")
plotSegmentUsageByCT(results,"TRB","j_call")

detectCircadianVDJ(results,out.dir = "../analysis/OS_VDJ/",immune.receptor = "tcr")
detectCircadianVDJ(results.BCR,out.dir = "../analysis/OS_VDJ/",immune.receptor = "bcr")

#compare expression by paired time point
#expression_mat_by_type_individual_CT<-srt2bulkMatrix(test,c("patient","type","CT"))
test<-readRDS('../analysis/all.samples.oddsremoved.tsneran.Azimuth.healthonly.rds')
#saveRDS(expression_mat_by_type_individual_CT,"../analysis/expression_mat_by_type_individual_CT.rds")
expression_mat_by_type_individual_CT<-readRDS("../analysis/expression_mat_by_type_individual_CT.rds")
mat.test<-expression_mat_by_type_individual_CT[,grepl("CT9$|CT17",colnames(expression_mat_by_type_individual_CT))&
                                                 grepl("CD14_Mono",colnames(expression_mat_by_type_individual_CT))]
mat.test<-mat.test[,colSums(mat.test)!=0]
mat.test<-mat.test[rowMin(mat.test %>% as.matrix())>=3,]
mat.test<-mat.test[!(grepl("^AC[0-9].*",rownames(mat.test))|grepl("^AL[0-9].*",rownames(mat.test))|grepl("^FP[0-9].*",rownames(mat.test))),]
condition<-factor(getField(colnames(mat.test),"-",3))
dornor<-factor(getField(colnames(mat.test),"-",1))
colData <- data.frame(row.names=colnames(mat.test), condition, dornor)

#dds<-DESeqDataSetFromMatrix(countData = mat.test, colData = colData, design = ~ condition+dornor)
dds<-DESeqDataSetFromMatrix(countData = mat.test, colData = colData, design = ~ condition)
dds.res <- DESeq(dds, fitType = 'mean')

plotdata<-results(dds.res) %>% as.data.frame()
plotdata$significant<-ifelse((-log10(plotdata$pvalue))>=2,"yes","no")
plotdata$gene<-rownames(plotdata)
ggplot(plotdata)+geom_point(aes(x=log2FoldChange,y=-log10(pvalue),color=significant))+
  geom_text_repel(data=plotdata[plotdata$significant=="yes",],aes(x=log2FoldChange,y=-log10(pvalue),label=gene))+
  geom_hline(yintercept=2,linetype=2)+
  geom_vline(xintercept=1,linetype=2)+
  geom_vline(xintercept=-1,linetype=2)+
  scale_color_manual(values = c("grey","red"))+theme_classic2()+NoLegend()
set3<-plotdata[plotdata$pvalue<0.01,"gene"]

#plotDispEsts(dds.res)
gene_ids=plotdata[plotdata$significant=="yes"&plotdata$log2FoldChange<0&!(plotdata$gene %in% c("HBB","HBA1","HBA2")),"gene"]
gene_ids=bitr(gene_ids,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID
gene_ids_universe=expression_mat_by_type_individual_CT[,grepl(gsub(" ","_","CD4 T"),colnames(expression_mat_by_type_individual_CT))]
gene_ids_universe=gene_ids_universe[rowMeans(gene_ids_universe)>1,] %>% rownames()
gene_ids_universe=bitr(gene_ids_universe,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID
#out=enrichGO(gene_ids,OrgDb=org.Hs.eg.db,universe=gene_ids_universe)
out=enrichPathway(gene_ids,universe=gene_ids_universe)
#out=enrichMKEGG(gene_ids)
dotplot(out)
out@result


#analyse individual-level JTK cycle results
JTK_CYCLEouts<-getAllJTK_CYCLEouts('../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5_by_individual2',celltypes=celltypes)
individual.cycIDs<-JTK_CYCLEouts[JTK_CYCLEouts$ADJ.P<0.01,"CycID"] %>% unique()

bulk_result<-read.delim('/lustre/home/acct-medll/medll/data/analysis/tauFisher_test/GSE113883_hs_full_adj_meta2d/meta2d_GSE113883_hs_full_adj.txt')
bulk_result_IDs<-bulk_result[bulk_result$JTK_pvalue<0.01,"CycID"]
JTK_CYCL<-NULL
for(block_data in list.files("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5_health/")){
  block_path=paste0("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5_health/",block_data,"/")
  block=getJTK_CYCLEouts(block_path,celltypes=celltypes,filter.data = F)
  if(is.null(JTK_CYCL)){
    JTK_CYCL=block
  }else{
    JTK_CYCL=rbind(JTK_CYCL,block)
  }
}
pseudo_bulk_IDs<-dplyr::filter(JTK_CYCL,ADJ.P<0.01)$CycID %>% unique()
JTK_CYCL_filter<-dplyr::filter(JTK_CYCL,ADJ.P<0.01)
#shared<-intersect(pseudo_bulk_IDs,bulk_result_IDs)
vennplot<-VennDiagram::venn.diagram(list("sampling_from_pool"=pseudo_bulk_IDs,"bulk"=bulk_result_IDs,"SC_individual"=individual.cycIDs),
                                    fill = generateColor(3),col = "black",margin=0.05,
                                    filename = NULL)
grid.draw(vennplot)
shared<-Reduce(intersect,x = list(pseudo_bulk_IDs,bulk_result_IDs,individual.cycIDs))
JTK_CYCLEouts[JTK_CYCLEouts$CycID%in%shared&JTK_CYCLEouts$ADJ.P<0.01,] %>% arrange(ADJ.P)


ggplot(JTK_CYCLEouts[-log10(JTK_CYCLEouts$ADJ.P)<=-log10(0.01),])+geom_point(aes(x=AMP,y=-log10(ADJ.P)))+
  geom_point(data=JTK_CYCLEouts[-log10(JTK_CYCLEouts$ADJ.P)>-log10(0.01),],aes(x=AMP,y=-log10(ADJ.P)),color="red")+
  #geom_text(data=JTK_CYCLEouts[-log10(JTK_CYCLEouts$ADJ.P)>-log10(0.01),],aes(x=AMP,y=-log10(ADJ.P),label=CycID),nudge_x = 0.5,nudge_y = 0.05)+
  geom_hline(yintercept=-log10(0.01),linetype=2)+facet_wrap(~celltype)

#reactome/enrichGO analysis
i=4
gene_ids=bitr(JTK_CYCLEouts[JTK_CYCLEouts$ADJ.P<0.01&JTK_CYCLEouts$celltype==celltypes[i],"CycID"],fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID
gene_ids_universe=expression_mat_by_type_individual_CT[,grepl(gsub(" ","_",celltypes[i]),colnames(expression_mat_by_type_individual_CT))]
gene_ids_universe=gene_ids_universe[rowMeans(gene_ids_universe)>1,] %>% rownames()
gene_ids_universe=bitr(gene_ids_universe,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID
#out=enrichGO(gene_ids,OrgDb=org.Hs.eg.db,universe=gene_ids_universe)
out=enrichPathway(gene_ids,universe=gene_ids_universe)
#out=enrichMKEGG(gene_ids)
dotplot(out)
out@result


#VDJ结果作图
JTK_CYCLEouts<-read.delim('../analysis/OS_VDJ/JTKresult_TRA_j_call.txt')
ggplot(JTK_CYCLEouts[log2(JTK_CYCLEouts$AMP)<=1|-log10(JTK_CYCLEouts$ADJ.P)<=-log10(0.05),])+geom_point(aes(x=log2(AMP+1),y=-log10(ADJ.P)))+
  geom_point(data=JTK_CYCLEouts[log2(JTK_CYCLEouts$AMP)>1&-log10(JTK_CYCLEouts$ADJ.P)>-log10(0.05),],aes(x=log2(AMP+1),y=-log10(ADJ.P)))+
  geom_text(data=JTK_CYCLEouts[log2(JTK_CYCLEouts$AMP)>1&-log10(JTK_CYCLEouts$ADJ.P)>-log10(0.05),],aes(x=log2(AMP+1),y=-log10(ADJ.P),label=CycID),nudge_x = 0.5,nudge_y = 0.05)+
  geom_hline(yintercept=-log10(0.05),linetype=2)+theme_classic2()+ggtitle("TRA_segment_J")

#melatonin data analysis
melatonin_data_test<-data.frame(concentrate=c(1000,500,250,125,62.5,31.25,15.63,0),OD=c(0.888,1.096,1.641,2.056,2.399,2.542,2.769,2.843))
ggplot(melatonin_data_test)+geom_point(aes(x=log2(concentrate),y=log10(OD)))+
  theme(panel.background=element_rect(fill = "white"),
        axis.line = element_line())
lm.melatonin<-lm(log2(concentrate)~log10(OD),melatonin_data_test[-8,])
plot(lm.melatonin)
new<-data.frame("OD"=c(1.745,1.708,1.623,1.455))
data<-2^predict(lm.melatonin,new) %>% as.data.frame()
colnames(data)<-"concentrate"
data$individual<-c(1,2,1,2)
data$CT<-c("CT10","CT10","CT16","CT16")
ggplot(data)+geom_boxplot(aes(x=CT,y=concentrate,group=CT))+
  #stat_compare_means(aes(x=CT,y=concentrate,group=CT),method = "t")+
  theme(panel.background=element_rect(fill = "white"),
        axis.line = element_line())+
  xlab("")+ylab("concentrate (pg/ml)")


#add NI aging labels
test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.Azimuth.rds")
test<-FindNeighbors(test, reduction = "integrated.cca", dims = 1:30)
test<-FindClusters(test, resolution = 0.8)
dimplot_publication(test,"RNA_snn_res.0.8",colors = generateColor(36,col.dist.min = 0.25),label = T)
#saveRDS(test,"../analysis/all.samples.oddsremoved.tsneran.Azimuth.rds")
all_markers<-FindAllMarkers(test,logfc.threshold = log(2),min.diff.pct = 0.3,min.pct = 0.5,only.pos = T)
#saveRDS(all_markers,"../analysis/all.samples.oddsremoved.tsneran.resolution0.5.markers.rds")
all_markers<-readRDS("../analysis/all.samples.oddsremoved.tsneran.resolution0.5.markers.rds")
dplyr::filter(all_markers,cluster==11) %>% arrange(desc(pct.1-pct.2))
FeaturePlot(test,"ANXA1")

test<-TypeCluster(test,c(9,25),type = "B_naive",col.name = "RNA_snn_res.0.5",new.meta = "manual_NI")
test<-TypeCluster(test,c(15),type = "Plasma cell",col.name = "RNA_snn_res.0.5",new.meta = "manual_NI")
test<-TypeCluster(test,c(18),type = "pDC",col.name = "RNA_snn_res.0.5",new.meta = "manual_NI")

test<-TypeCluster(test,c(1,8),type = "CD14⁺Monocytes",col.name = "RNA_snn_res.0.5",new.meta = "manual_NI")
test<-TypeCluster(test,c(11),type = "CD16⁺Monocytes",col.name = "RNA_snn_res.0.5",new.meta = "manual_NI")
test<-TypeCluster(test,c(13,23),type = "mDC",col.name = "RNA_snn_res.0.5",new.meta = "manual_NI")

test<-TypeCluster(test,c(21),type = "HSC_CD34",col.name = "RNA_snn_res.0.5",new.meta = "manual_NI")

dimplot_publication(test,"manual_NI",colors = generateColor(33,col.dist.min = 0.25),label = T)
#some of cell type in NI paper did not show up, subset TNK for further annotation
test.TNK<-readRDS("../analysis/all.samples.TNK.integrated.bottom.level.rds")
dimplot_publication(test.TNK,group.by = "seurat_clusters",colors = generateColor(35,rgb.max = 0.85,alpha = 0.8,col.dist.min = 0.25,seed = 2025),label = T)
all_markers<-FindAllMarkers(test.TNK,logfc.threshold = log(1.5),min.diff.pct = 0.2,min.pct = 0.5,only.pos = T)
saveRDS(all_markers,"../analysis/tmpTNK.markers.rds")
all_markers<-readRDS("../analysis/tmpTNK.markers.rds")
dplyr::filter(all_markers,cluster==13) %>% arrange(desc(pct.1-pct.2))
VlnPlot(test.TNK,rev(c(TNK_markers_NI,"CD14","MS4A1","JCHAIN")),pt.size = 0,stack = T,flip = T,fill.by = "ident")+
  scale_fill_manual(values=generateColor(34,rgb.max = 0.85,alpha = 0.8,col.dist.min = 0.25,seed = 2025))+
  theme(strip.text.y=element_text(size=8,face="plain",hjust =0),axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank())+scale_x_discrete(limits=as.character(c(20,19,0,2,6,18,33,27,5,10,8,29,15,21,25,
                                                                       1,31,30,32,3,14,16,9,11,28,13,12,
                                                                       4,22,7,24,
                                                                      17,23,
                                                                       26)))
#about 30 min
test.TNK<-RunAzimuth(test.TNK,reference = "pbmcref")

test.TNK<-TypeCluster(test.TNK,c(6,0,33),type = "CD4_naive_CCR7",col.name = "seurat_clusters",new.meta = "manual_NI")
test.TNK<-TypeCluster(test.TNK,c(19,2,18,27,25),type = "CD4_TCM_AQP3",col.name = "seurat_clusters",new.meta = "manual_NI")
test.TNK<-TypeCluster(test.TNK,c(5,10),type = "CD4_TEM_ANXA1",col.name = "seurat_clusters",new.meta = "manual_NI")
test.TNK<-TypeCluster(test.TNK,c(15),type = "CD4_TEM_GNLY",col.name = "seurat_clusters",new.meta = "manual_NI")
test.TNK<-TypeCluster(test.TNK,c(21),type = "CD4_Treg_FOXP3",col.name = "seurat_clusters",new.meta = "manual_NI")
test.TNK<-TypeCluster(test.TNK,c(8,29),type = "CD4_TEM_GZMK",col.name = "seurat_clusters",new.meta = "manual_NI")

test.TNK<-TypeCluster(test.TNK,c(1,31),type = "CD8_naive_LEF1",col.name = "seurat_clusters",new.meta = "manual_NI")
test.TNK<-TypeCluster(test.TNK,c(3,14),type = "CD8_TEM_CMC1",col.name = "seurat_clusters",new.meta = "manual_NI")
test.TNK<-TypeCluster(test.TNK,c(16,9,11,28,23),type = "CD8_TEM_GNLY",col.name = "seurat_clusters",new.meta = "manual_NI")
test.TNK<-TypeCluster(test.TNK,c(30),type = "dnT_LYST",col.name = "seurat_clusters",new.meta = "manual_NI")
test.TNK<-TypeCluster(test.TNK,c(32),type = "CD8_TEM_ZNF683",col.name = "seurat_clusters",new.meta = "manual_NI")
test.TNK<-TypeCluster(test.TNK,c(13),type = "CD8_MAIT_SLC4A10",col.name = "seurat_clusters",new.meta = "manual_NI")
test.TNK<-TypeCluster(test.TNK,c(12),type = "γδT",col.name = "seurat_clusters",new.meta = "manual_NI")
test.TNK<-TypeCluster(test.TNK,c(7,24),type = "NK_CD56ᵈⁱᵐ",col.name = "seurat_clusters",new.meta = "manual_NI")
test.TNK<-TypeCluster(test.TNK,c(4,22),type = "NK_CD56ᵇʳⁱᵍʰᵗ",col.name = "seurat_clusters",new.meta = "manual_NI")

test.TNK<-TypeCluster(test.TNK,c(26),type = "TNK_proliferatig_MKI67",col.name = "seurat_clusters",new.meta = "manual_NI")
test.TNK<-TypeCluster(test.TNK,c(17,20,23),type = "debris",col.name = "seurat_clusters",new.meta = "manual_NI")

dimplot_publication(test.TNK,group.by = "manual_NI")

saveRDS(test.TNK,"../analysis/all.samples.TNK.integrated.bottom.level.rds")

tnk.meta<-test.TNK$manual_NI

test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.Azimuth.rds")
test<-AddMetaData(test,metadata = tnk.meta,col.name = "manual_NI")
dimplot_publication(test,group.by = "manual_NI")
VlnPlot(test,rev(c(TNK_markers_NI,"CD14","MS4A1","JCHAIN")),group.by = "manual_NI",pt.size = 0,stack = T,flip = T,fill.by = "ident")+
  scale_fill_manual(values=generateColor(34,rgb.max = 0.85,alpha = 0.8,col.dist.min = 0.25,seed = 2025))+
  theme(strip.text.y=element_text(size=8,face="plain",hjust =0),axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),axis.text.x = element_text(angle=90,size=8))
saveRDS(test,"../analysis/all.samples.oddsremoved.tsneran.Azimuth.rds")

VlnPlot(test,rev(c(TNK_markers_NI,"CD14","MS4A1","JCHAIN")),pt.size = 0,stack = T,flip = T,fill.by = "ident")+
  scale_fill_manual(values=generateColor(25,rgb.max = 0.85,alpha = 0.8,col.dist.min = 0.25,seed = 2025))+
  theme(strip.text.y=element_text(size=8,face="plain",hjust =0),axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank())

#network showing gene sharing between cell types
JTK_CYCL<-getAllJTK_CYCLEouts("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5_by_individual2",celltypes = celltypes)
JTK.result.filtered<-JTK_CYCL[JTK_CYCL$ADJ.P<=0.05,]
test$type<-test$predicted.celltype.l1.5
plotNetWorkCircadian(test,JTK.result.filtered)

#shift phase of individual that have a different pattern among population




#sychronization plot
test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.Azimuth.healthonly.rds")
JTK_CYCL<-getAllJTK_CYCLEouts("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5_by_individual2",celltypes = celltypes)
JTK.result.filtered<-JTK_CYCL[JTK_CYCL$ADJ.P<=0.01,]
test$type<-test$predicted.celltype.l1.5
JTK.result.filtered$CycID %>% table() %>% sort %>% tail(.,20)
JTK.result.filtered[JTK.result.filtered$CycID=="PNPLA6",]
plotRadarCellType(test,JTK.result.filtered[JTK.result.filtered$CycID=="CXCR4","celltype"],"CXCR4")
plotRadarCellType(test,JTK.result.filtered[JTK.result.filtered$CycID=="PNPLA6","celltype"],"PNPLA6")

plotRadarCellType(test,c("ASDC","B memory","B naive","CD16 Mono","CD4 T","CD8 T","cDC2","gdT","NK"),"CXCR4")
plotRadarCellType(test,c("CD14 Mono","CD4 T","CD8 T","gdT","NK"),"CCR7")

#gene enrichment use topGO
test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.Azimuth.healthonly.rds")
test$type<-test$predicted.celltype.l1.5
JTK_CYCL<-getAllJTK_CYCLEouts("../analysis/OSgenePredict_PseudoBulk_20k_all_genes_level1.5_by_individual2",celltypes = celltypes)
JTK.result.filtered<-JTK_CYCL[JTK_CYCL$ADJ.P<=0.05,]
specific.genes<-getOnceOccured(JTK.result.filtered$CycID)
JTK.result.filtered<-JTK.result.filtered[JTK.result.filtered$CycID %in% specific.genes,]
result_list<-list()
for(celltype in celltypes){
  celltype.specific.genes=JTK.result.filtered[JTK.result.filtered$celltype==celltype&JTK.result.filtered$ADJ.P<=0.01,"CycID"]
  celltype.genes=FetchExpressedGenesCellType(test,celltype)
  results=enrichGObyHGNC(celltype.specific.genes,celltype.genes)
  result_list[[celltype]]=results
}
saveRDS(result_list,"../analysis/enrichGO.rds")

for(i in 1:length(celltypes)){
  message(celltypes[i])
  try(goplot(result_list[[i]],showCategory = 8))
  Sys.sleep(5)
}

#run viper to show tf activities
small.test<-subset(test,patient=="TFSH190500A_HZD")
expression_data <- as.matrix(GetAssayData(small.test, slot = "data"))
# Load human regulons (for example)
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs
# Run VIPER
tf_activities <- run_viper(expression_data, regulons, 
                           options = list(method = "scale", minsize = 4, 
                           eset.filter = FALSE, cores = 39, verbose = FALSE))
#read obj that contain TF activity info
test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.Azimuth.healthonly.rds")



#cell proliferation time preferance in a day
test<-readRDS("../analysis/all.samples.oddsremoved.tsneran.Azimuth.healthonly.rds")
proliferation_genes <- c("MKI67", "PCNA", "TOP2A", "CCNB1", "CCNB2", "CCND1", "CCND2", "CCNE1", "CCNE2", "CDK1", "CDK2", "CDK4", "CDK6")
test <- AddModuleScore(
  test,
  features = list(proliferation_genes),
  name = "Proliferation_Score"
)
plotPseudobulkByCTByIndividual(test,"CD14 Mono","MKI67")

#use seacell(metacell) to analyse data
test<-readAllMetaCell()
samples=Cells(test) %>% getField(.,"_",1:3)
names(samples)=Cells(test)
test<-AddMetaData(test,samples,col.name = "sample")
test[["RNA"]] <- split(test[["RNA"]], f = test$sample)
test<-NormalizeData(test)
test<-FindVariableFeatures(test,nfeatures = 5000)
test<-ScaleData(test)
test<-RunPCA(test)
test<-FindNeighbors(test)
test<-FindClusters(test,resolution = 1)
test<-RunUMAP(test,dims = 1:30,reduction.name = "umap.unintegrated")
DimPlot(test, reduction = "umap.unintegrated", group.by = c("sample", "type"))

test <- IntegrateLayers(
  object = test, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
test <- FindNeighbors(test, reduction = "harmony", dims = 1:30)
test <- FindClusters(test, resolution = 2, cluster.name = "harmony_clusters")
test <- RunUMAP(test, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(test, reduction = "umap.harmony", group.by = c("sample", "type"))

## demux 4 sample mix
source('/tmpdata/LyuLin/script/demuxer/src/demux_core_functions.R')
sample_path="/tmpdata/LyuLin/analysis/circadian/cellranger/X5_250424MIX07/outs/per_sample_outs/X5_250424MIX07/count/sample_filtered_feature_bc_matrix/"
srt<-seuratWrap1(sample_path,min.features = 500)
srt<-seuratWrap2(srt)
srt<-seuratWrap3(srt,res = 0.01)
DimPlot(srt)

#sc_merged0<-readAFFileWithFiltering('/tmpdata/LyuLin/analysis/circadian/cellranger/X5_250424MIX02/outs/per_sample_outs/X5_250424MIX02/count',sc.mode = T)
#sc_merged<-sc_merged0
#std<-readAFStandard('/tmpdata/LyuLin/analysis/circadian/RNASeq')
#sc_merged<-filterInformativePositionsBySTD(sc_merged,std[c(1,2,3,4,6)])
#sc_merged<-filterDoublet(sc_merged,threshould = 1)
#sc_merged<-callExactBase(sc_merged)
#demux_out<-demuxBySTD(sc_merged,std[c(1,2,3,4,6)])
demux_out<-demux('/tmpdata/LyuLin/analysis/circadian/cellranger/X5_250424MIX07/outs/per_sample_outs/X5_250424MIX07/count','/tmpdata/LyuLin/analysis/circadian/RNASeq')

srt<-AddMetaData(srt,demux_out,col.name = "individual")
srt@meta.data[is.na(srt@meta.data$individual),"individual"]<-"doublet"
srt@meta.data$CT<-case_when(srt@meta.data$individual == "KD" ~ "CT19",
                            srt@meta.data$individual == "ZYR" ~ "CT19",
                            srt@meta.data$individual == "JJC" ~ "CT11",
                            srt@meta.data$individual == "LYH" ~ "CT19",
                            TRUE ~ "unknown")
new_cell_ids <- paste0("TF_",srt@meta.data$individual,"_",srt@meta.data$CT,"_", colnames(srt))
srt@meta.data$orig.ident=rownames(srt@meta.data)
srt<-RenameCells(srt,new.names = new_cell_ids)
DimPlot(srt,group.by = "individual")
saveRDS(srt,'/tmpdata/LyuLin/analysis/circadian/R/X5_250424MIX07.demuxed.rds')

srt$individual %>% table %>% as.data.frame() %>% ggplot() + 
  geom_bar(aes(x=get("."),y=Freq),stat="identity",fill="blue",color="black",width=0.75,)+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  xlab("")+ylab("# of cell")+scale_x_discrete(limit=c("JJC","KD","LYH","ZYR","unknown","doublet"))+
  geom_text(aes(x=get("."),y=Freq+1000,label=Freq))+scale_y_continuous(expand = c(0,0),limits = c(0,20000))+
  ggtitle("sample: X5_250424MIX06\ndemultiplex with SNP on chrM\nedit dist = 0")+
  theme(title = element_text(size = 8,face = "plain",family="Arial"),axis.title = element_text(size=12))
DimPlot(srt,group.by = "individual")
srt@meta.data %>% ggplot(.)+geom_violin(aes(x=individual,y=nCount_RNA,fill=individual))+ylim(c(0,30000))+
  scale_x_discrete(limit=c("doublet","JJC","KD","LYH","ZYR","unknown"))+theme_linedraw()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))

srt@meta.data %>% ggplot(.)+geom_violin(aes(x=individual,y=nFeature_RNA,fill=individual))+
  scale_x_discrete(limit=c("doublet","JJC","KD","LYH","ZYR","unknown"))+theme_linedraw()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))

# use single sample to see performance of filterDoublet
sc_merged0<-readAFFileWithFiltering('~/Downloads/ZXL.allele.freq.cell.tsv',sc.mode = T)
sc_merged<-sc_merged0
std<-readAFStandard('/tmpdata/LyuLin/analysis/circadian/RNASeq')
sc_merged<-filterInformativePositionsBySTD(sc_merged,std[c(1,2,3,4,6)])
sc_merged<-filterDoublet(sc_merged,plt.db.score = c(-0.1,20))

# manual annotation level1, multiplexed sample (batch3)
setwd('/tmpdata/LyuLin/analysis/circadian/R')
srt<-mergeMultiplexedSamplesByID(
  dir.path='/tmpdata/LyuLin/analysis/circadian/R',
  sample.ids=c("X5_250409MIX01","X5_250424MIX02","X5_250424MIX03",
               "X5_250424MIX04","X5_250424MIX05","X5_250424MIX06",
               "X5_250424MIX07")
  )
srt<-dimentionalReductionMergedSrt(srt)
DimPlot(srt,label = T)
# forget to add sample column in previous step, manually add it, and previous pipeline has been modified
#srt$sample<-case_when(paste0(srt$individual,srt$CT) %in% c("KDCT15","ZYRCT15","JJCCT19","LYHCT15") ~ "X5_250409MIX01",
#                      paste0(srt$individual,srt$CT) %in% c("KDCT35","ZYRCT27","JJCCT15","LYHCT27") ~ "X5_250424MIX02",
#                      paste0(srt$individual,srt$CT) %in% c("KDCT23","ZYRCT23","JJCCT35","LYHCT31") ~ "X5_250424MIX03",
#                      paste0(srt$individual,srt$CT) %in% c("KDCT31","ZYRCT11","JJCCT31","LYHCT23") ~ "X5_250424MIX04",
#                      paste0(srt$individual,srt$CT) %in% c("KDCT27","ZYRCT31","JJCCT23","LYHCT11") ~ "X5_250424MIX05",
#                      paste0(srt$individual,srt$CT) %in% c("KDCT11","ZYRCT35","JJCCT27","LYHCT35") ~ "X5_250424MIX06",
#                      paste0(srt$individual,srt$CT) %in% c("KDCT19","ZYRCT19","JJCCT11","LYHCT19") ~ "X5_250424MIX07")
srt<-integrateMergedSrt(srt)
saveRDS(srt,"/tmpdata/LyuLin/analysis/circadian/R/7mixed.integrated.rds")

DimPlot(srt,label = T)
FeaturePlot(srt,MAIN_GROUP_MARKERS)
srt<-TypeCluster(srt,cluster = c(0,2,3,4,6,7,9,10,11),type = "TNK",new.meta = "manual.level1")
srt<-TypeCluster(srt,cluster = c(5,12),type = "B",new.meta = "manual.level1")
srt<-TypeCluster(srt,cluster = c(1,8,13,15),type = "Myeloid",new.meta = "manual.level1")
srt<-TypeCluster(srt,cluster = c(14),type = "Stem",new.meta = "manual.level1")
DimPlot(srt,label = T,group.by = "manual.level1")

srt.B<-subset(srt,manual.level1=="B")
srt.B<-dimentionalReductionSubsetSrt(srt.B)
srt.B<-integrateSubset(srt.B)
DimPlot(srt.B,label = T)
srt.B<-TypeCluster(srt.B,cluster = c(4,7,9),type = "doublet",new.meta = "manual.level2")
srt.B<-TypeCluster(srt.B,cluster = c(4,7,9),type = "doublet",new.meta = "manual.level1")
srt.B<-TypeCluster(srt.B,cluster = c(1),type = "cB01_TCL1A+B_naive",new.meta = "manual.level2")
srt.B<-TypeCluster(srt.B,cluster = c(0,2,3,8),type = "cB03_IgG+IgA+B_memory",new.meta = "manual.level2")
srt.B<-TypeCluster(srt.B,cluster = c(5,10,11),type = "cB04_Plasma",new.meta = "manual.level2")
srt.B<-TypeCluster(srt.B,cluster = c(6),type = "debris",new.meta = "manual.level2")
srt.B<-TypeCluster(srt.B,cluster = c(6),type = "debris",new.meta = "manual.level1")
DimPlot(srt.B,label = T,group.by = "manual.level2")

srt.B<-TypeCluster(srt.B,cluster = c(4,7,9),type = "doublet",new.meta = "manual_NI")
srt.B<-TypeCluster(srt.B,cluster = c(1),type = "B_naive",new.meta = "manual_NI")
srt.B<-TypeCluster(srt.B,cluster = c(0,2,3,8),type = "B_memory",new.meta = "manual_NI")
srt.B<-TypeCluster(srt.B,cluster = c(5,10,11),type = "Plasma cell",new.meta = "manual_NI")
srt.B<-TypeCluster(srt.B,cluster = c(6),type = "debris",new.meta = "manual_NI")
DimPlot(srt.B,label = T,group.by = "manual_NI")
saveRDS(srt.B,"/tmpdata/LyuLin/analysis/circadian/R/7mixed.integrated.B.rds")


srt.T<-subset(srt,manual.level1=="TNK")
srt.T<-dimentionalReductionSubsetSrt(srt.T)
srt.T<-integrateSubset(srt.T)

srt.T<-FindClusters(srt.T,resolution =1 )
DimPlot(srt.T,label = T)

srt.T<-TypeCluster(srt.T,cluster = c(24),type = "doublet",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(24),type = "doublet",new.meta = "manual.level1")
srt.T<-TypeCluster(srt.T,cluster = c(0,1,12,25,33),type = "cTNK01_CD4+CCR7+Tn",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(3,10),type = "cTNK02_CD8+CCR7+Tn",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(4,6,8,39),type = "cTNK03_CD4+PTGER2+Tcm",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(11,16,37),type = "cTNK04_CD8+GZMK+Tem",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(2,9,13,17,19,20,27),type = "cTNK05_GNLY+NK",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(7,38),type = "cTNK06_gdT",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(5,32,36),type = "cTNK07_CD8+SLC4A10+MAIT",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(30),type = "cTNK08_CD4+NKG7+CTL",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(21),type = "cTNK09_CD4+FOXP3+Treg",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(26),type = "cTNK10_MKI67+proliferating",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(28),type = "cTNK11_LYST+dnT",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(29),type = "cTNK12_CD8+cTRBV+Tcm-like",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(22),type = "cTNK13_CD8+IKZF2+IEL",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(14,18,23,31,34),type = "cTNK14_CD8+GZMH+TEMRA/TEFF",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(15,35),type = "debris",new.meta = "manual.level2")
srt.T<-TypeCluster(srt.T,cluster = c(15,35),type = "debris",new.meta = "manual.level1")
DimPlot(srt.T,label = T,group.by = "manual.level2")

srt.T<-TypeCluster(srt.T,cluster = c(24),type = "doublet",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(1,33),type = "CD4_TCM_AQP3",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(6,8,39),type = "CD4_TEM_ANXA1",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(30),type = "CD4_TEM_GNLY",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(4),type = "CD4_TEM_GZMK",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(21),type = "CD4_Treg_FOXP3",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(0,12,25),type = "CD4_naive_CCR7",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(5,32,36),type = "CD8_MAIT_SLC4A10",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(11,16,29,37),type = "CD8_TEM_CMC1",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(14,18,23,31,34),type = "CD8_TEM_GNLY",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(22),type = "CD8_TEM_ZNF683",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(3,10),type = "CD8_naive_LEF1",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(17,20),type = "NK_CD56ᵇʳⁱᵍʰᵗ",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(2,9,13,19,27),type = "NK_CD56ᵈⁱᵐ",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(26),type = "TNK_proliferatig_MKI67",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(28),type = "dnT_LYST",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(7,38),type = "γδT",new.meta = "manual_NI")
srt.T<-TypeCluster(srt.T,cluster = c(15,35),type = "debris",new.meta = "manual_NI")
DimPlot(srt.T,label = T,group.by = "manual_NI")
saveRDS(srt.T,"/tmpdata/LyuLin/analysis/circadian/R/7mixed.integrated.TNK.rds")

srt.M<-subset(srt,manual.level1=="Myeloid")
srt.M<-dimentionalReductionSubsetSrt(srt.M)
srt.M<-integrateSubset(srt.M)
DimPlot(srt.M,label = T)

srt.M<-TypeCluster(srt.M,cluster = c(5,6,7),type = "doublet",new.meta = "manual.level2")
srt.M<-TypeCluster(srt.M,cluster = c(5,6,7),type = "doublet",new.meta = "manual.level1")
srt.M<-TypeCluster(srt.M,cluster = c(0,1,2,8),type = "cM01_CD14+FCGR3A-Monocyte",new.meta = "manual.level2")
srt.M<-TypeCluster(srt.M,cluster = c(3),type = "cM02_CD14-FCGR3A+Monocyte",new.meta = "manual.level2")
srt.M<-TypeCluster(srt.M,cluster = c(11),type = "cM03_CLEC9A+cDC1",new.meta = "manual.level2")
srt.M<-TypeCluster(srt.M,cluster = c(4),type = "cM04_FCER1A+cDC2",new.meta = "manual.level2")
srt.M<-TypeCluster(srt.M,cluster = c(10),type = "cM05_LILRA4+pDC",new.meta = "manual.level2")
srt.M<-TypeCluster(srt.M,cluster = c(9),type = "debris",new.meta = "manual.level2")
srt.M<-TypeCluster(srt.M,cluster = c(9),type = "debris",new.meta = "manual.level1")
DimPlot(srt.M,label = T,group.by = "manual.level2")

srt.M<-TypeCluster(srt.M,cluster = c(5,6,7),type = "doublet",new.meta = "manual_NI")
srt.M<-TypeCluster(srt.M,cluster = c(5,6,7),type = "doublet",new.meta = "manual_NI")
srt.M<-TypeCluster(srt.M,cluster = c(0,1,2,8),type = "CD14⁺Monocytes",new.meta = "manual_NI")
srt.M<-TypeCluster(srt.M,cluster = c(3),type = "CD16⁺Monocytes",new.meta = "manual_NI")
srt.M<-TypeCluster(srt.M,cluster = c(4,11),type = "mDC",new.meta = "manual_NI")
srt.M<-TypeCluster(srt.M,cluster = c(10),type = "pDC",new.meta = "manual_NI")
srt.M<-TypeCluster(srt.M,cluster = c(9),type = "debris",new.meta = "manual_NI")
srt.M<-TypeCluster(srt.M,cluster = c(9),type = "debris",new.meta = "manual_NI")
DimPlot(srt.M,label = T,group.by = "manual_NI")

saveRDS(srt.M,"/tmpdata/LyuLin/analysis/circadian/R/7mixed.integrated.M.rds")

srt.S<-subset(srt,manual.level1=="Stem")
srt.S<-dimentionalReductionSubsetSrt(srt.S)
srt.S<-integrateSubset(srt.S)
DimPlot(srt.S,label = T)

srt.S<-TypeCluster(srt.S,cluster = c(0,4,5,6),type = "doublet",new.meta = "manual.level2")
srt.S<-TypeCluster(srt.S,cluster = c(0,4,5,6),type = "doublet",new.meta = "manual.level1")
srt.S<-TypeCluster(srt.S,cluster = c(1,2,3),type = "cS01_CD34+HSPC",new.meta = "manual.level2")
srt.S<-TypeCluster(srt.S,cluster = c(7),type = "cM05_MS4A2+Granulocyte",new.meta = "manual.level2")
DimPlot(srt.S,label = T,group.by = "manual.level2")

srt.S<-TypeCluster(srt.S,cluster = c(0,4,5,6),type = "doublet",new.meta = "manual_NI")
srt.S<-TypeCluster(srt.S,cluster = c(1,2,3),type = "HSC_CD34",new.meta = "manual_NI")
srt.S<-TypeCluster(srt.S,cluster = c(7),type = "MS4A2_Granulocyte",new.meta = "manual_NI")
DimPlot(srt.S,label = T,group.by = "manual_NI")

saveRDS(srt.S,"/tmpdata/LyuLin/analysis/circadian/R/7mixed.integrated.S.rds")

setwd('/tmpdata/LyuLin/analysis/circadian/R/')
generateAnnotationFile(TNK.file = '7mixed.integrated.TNK.rds',B.file = '7mixed.integrated.B.rds',
                       Myeloid.file = '7mixed.integrated.M.rds',Stem.file = '7mixed.integrated.S.rds',
                       cols = c("manual.level1","manual.level2","manual_NI"),
                       out.path = "cell.annotation.7mixed.integrated.tsv")
 
annotation.meta<-read.delim('cell.annotation.7mixed.integrated.tsv',header = T)
rownames(annotation.meta)<-annotation.meta$cell_id
srt<-AddMetaData(srt,annotation.meta)
srt$type<-srt$manual.level2

#found that LYH CT31 was wrongly assingned as CT21
srt@meta.data[srt@meta.data$CT=="CT21","CT"]<-"CT31"
old.cell.ids<-Cells(srt)
new.cell.ids<-gsub("CT21","CT31",old.cell.ids)
srt<-RenameCells(srt,new.names = new.cell.ids)
srt@meta.data[c("CT","individual")] %>% table()

saveRDS(srt,"7mixed.integrated.annotated.raw.rds")

srt<-subset(srt,doublet_score<=1&manual.level1!="doublet"&manual.level1!="debris")
saveRDS(srt,"7mixed.integrated.annotated.clean.rds")

srt<-SCTransform(srt)
saveRDS(srt,"7mixed.integrated.annotated.clean.sct.rds")

srt<-readRDS('7mixed.integrated.annotated.clean.sct.rds')
srt<-RunAzimuth(srt, reference = "pbmcref")
srt<-modAzimuthAnnotation(srt)
saveRDS(srt,"7mixed.integrated.annotated.clean.sct.Azimuth.rds")

srt$type<-srt$predicted.celltype.l1.5
DefaultAssay(srt)<-"RNA"
plotPseudobulkByCTByIndividual(srt,"CD16 Mono","NR1D2",normalize.data = T,
                               time.points = c("CT11","CT15","CT19","CT23","CT27","CT31","CT35"))

# generate cell annotation file for seacell
generateAnnotationFile(full.file = srt,cols=c("manual.level1","manual.level2","predicted.celltype.l2","manual_NI"),
                       out.path = "/tmpdata/LyuLin/analysis/circadian/R/cell.annotation.clean.sct.Azimuth.tsv")

# read metacell 
srt.metacell<-readAllMetaCell('/tmpdata/LyuLin/analysis/circadian/R/preparation_seacell.0.04/',std.cellranger.out = F)
#srt.metacell$sample<-case_when(paste0(srt.metacell$individual,srt.metacell$CT) %in% c("KDCT15","ZYRCT15","JJCCT19","LYHCT15") ~ "X5_250409MIX01",
#                      paste0(srt.metacell$individual,srt.metacell$CT) %in% c("KDCT35","ZYRCT27","JJCCT15","LYHCT27") ~ "X5_250424MIX02",
#                      paste0(srt.metacell$individual,srt.metacell$CT) %in% c("KDCT23","ZYRCT23","JJCCT35","LYHCT31") ~ "X5_250424MIX03",
#                      paste0(srt.metacell$individual,srt.metacell$CT) %in% c("KDCT31","ZYRCT11","JJCCT31","LYHCT23") ~ "X5_250424MIX04",
#                      paste0(srt.metacell$individual,srt.metacell$CT) %in% c("KDCT27","ZYRCT31","JJCCT23","LYHCT11") ~ "X5_250424MIX05",
#                      paste0(srt.metacell$individual,srt.metacell$CT) %in% c("KDCT11","ZYRCT35","JJCCT27","LYHCT35") ~ "X5_250424MIX06",
#                      paste0(srt.metacell$individual,srt.metacell$CT) %in% c("KDCT19","ZYRCT19","JJCCT11","LYHCT19") ~ "X5_250424MIX07")
srt.metacell$sample<-paste0(srt.metacell$individual,"_",srt.metacell$CT)
srt.metacell<-subset(srt.metacell,predicted.celltype.l2.purity>=0.90)
srt.metacell<-integrateSubset(srt.metacell,k.weight = 50)
DimPlot(srt.metacell,group.by = "type")
plotMetaCellByIndividual(srt.metacell,cell.type = "CD14 Mono",feature = "NR1D2",norm.dist = T)

plotPseudobulkByCTByIndividual(srt,"CD14 Mono","NR1D2",
                               time.points = c("CT11","CT15","CT19","CT23","CT27","CT31","CT35"),
                               normalize.data = T,plot = "errorbar",proportion = 0.5)
saveRDS(srt.metacell,"seacell.0.04.rds")

#use drc package to analyse ELISA data
## melatonin data
std_data <- data.frame(
  #Concentration = c(1000, 500, 250, 125, 62.5, 31.25, 15.63),  # Standard concentrations
  #OD = c(0.701, 1.078, 1.450, 1.735, 1.964, 2.200, 2.694)      # Corresponding OD values
  Concentration = c(1000, 500, 250, 125, 62.5, 31.25, 15.63, 0),
  OD = c(0.738, 1.151, 1.686, 2.411, 2.968, 3.101, 3.264, 3.259) 
)
fit <- drm(OD ~ Concentration, data = std_data, fct = LL.4())
summary(fit)

all_melatonin_data<-NULL
# Example unknown ODs
unknown_ODs<-c(1.372, 1.295, 1.133, 1.208, 1.316, 1.338, 1.300) # KD
unknown_ODs<-c(1.333, 1.219, 1.405, 1.335, 1.329, 1.382, 1.354) # WMY
unknown_ODs<-c(1.255, 1.149, 1.270, 1.092, 1.228, 1.297, 1.221) # ZYR
unknown_ODs1<-c(2.354, 2.138, 2.112, 2.107, 2.116, 2.224, 2.161) # LYH1
unknown_ODs2<-c(2.482, 2.375, 2.431, 2.141, 2.404, 2.185, 2.088) # LYH2
unknown_ODs<-(unknown_ODs1+unknown_ODs2)/2
unknown_ODs1<-c(2.077, 2.207, 1.795, 1.979, 2.067, 2.589, 2.313) # JJC1
unknown_ODs2<-c(2.583, 1.989, 2.027, 2.239, 2.097, 1.963, 2.234) # JJC2
unknown_ODs<-(unknown_ODs1+unknown_ODs2)/2
# Predict concentrations (inverse prediction)
unknown_conc <- ED(fit, unknown_ODs, type = "absolute", display = FALSE)
unknown_conc<-as.data.frame(unknown_conc)
unknown_conc$CT<-c("CT11","CT15","CT19","CT23","CT27","CT31","CT35")
unknown_conc$individual<-"KD"
unknown_conc$individual<-"WMY"
unknown_conc$individual<-"ZYR"
unknown_conc$individual<-"LYH"
unknown_conc$individual<-"JJC"
if(is.null(all_melatonin_data)){
  all_melatonin_data=unknown_conc
}else{
  all_melatonin_data=rbind(all_melatonin_data,unknown_conc)
}
ggplot(all_melatonin_data)+geom_point(aes(x=CT,y=Estimate))+
  geom_line(aes(x=CT,y=Estimate,group=individual,color=individual))+
  facet_wrap(~individual,scale="free_y")+
  ggtitle("melatonin of individuals")+theme(axis.text.x = element_text(angle=60,hjust=1))
all_melatonin_data<-all_melatonin_data %>% group_by(individual) %>% mutate(relative_concentration=(Estimate-mean(Estimate))/sd(Estimate))
ggplot(all_melatonin_data)+geom_boxplot(aes(x=CT,y=relative_concentration))+
  geom_point(aes(x=CT,y=relative_concentration))+ylab("z-score")+
  ggtitle("z-score of melatonin")+theme(axis.text.x = element_text(angle=60,hjust=1))
saveRDS(all_melatonin_data,"ELISA_melatonin.rds")
## cortisol data
std_data <- data.frame(
  Concentration = c(400, 200, 100, 50, 25, 12.5, 6.25, 0),  # Standard concentrations
  OD = c(0.079, 0.094, 0.107, 0.239, 0.368, 0.574, 0.736, 0.782)      # Corresponding OD values
)
fit <- drm(OD ~ Concentration, data = std_data, fct = LL.4())
summary(fit)

all_cortisol_data<-NULL
unknown_ODs<-c(0.198, 0.229, 0.317, 0.357, 0.284, 0.302, 0.272) # KD
unknown_ODs<-c(0.157, 0.253, 0.258, 0.227, 0.231, 0.204, 0.196) # WMY
unknown_ODs<-c(0.200, 0.199, 0.214, 0.269, 0.183, 0.188, 0.278) # ZYR
unknown_ODs<-c(0.440, 0.400, 0.531, 0.675, 0.351, 0.428, 0.344) # LYH
unknown_ODs<-c(0.088, 0.097, 0.132, 0.105, 0.130, 0.084, 0.087) # JJC

unknown_conc <- ED(fit, unknown_ODs, type = "absolute", display = FALSE)
unknown_conc<-as.data.frame(unknown_conc)
unknown_conc$CT<-c("CT11","CT15","CT19","CT23","CT27","CT31","CT35")
unknown_conc$individual<-"KD"
unknown_conc$individual<-"WMY"
unknown_conc$individual<-"ZYR"
unknown_conc$individual<-"LYH"
unknown_conc$individual<-"JJC"
if(is.null(all_cortisol_data)){
  all_cortisol_data=unknown_conc
}else{
  all_cortisol_data=rbind(all_cortisol_data,unknown_conc)
}
ggplot(all_cortisol_data)+geom_point(aes(x=CT,y=Estimate))+
  geom_line(aes(x=CT,y=Estimate,group=individual,color=individual))+
  facet_wrap(~individual,scale="free_y")+
  ggtitle("cortisol of individuals")+theme(axis.text.x = element_text(angle=60,hjust=1))
all_cortisol_data<-all_cortisol_data %>% group_by(individual) %>% mutate(relative_concentration=(Estimate-mean(Estimate))/sd(Estimate))
ggplot(all_cortisol_data)+geom_boxplot(aes(x=CT,y=relative_concentration))+
  geom_point(aes(x=CT,y=relative_concentration))+ylab("z-score")+
  ggtitle("z-score of cortisol")+theme(axis.text.x = element_text(angle=60,hjust=1))
saveRDS(all_cortisol_data,"ELISA_cortisol.rds")

# analyse metacell result
AllJTKresult<-readJTKFromMetaCells('/tmpdata/LyuLin/analysis/circadian/R/seacell.meta2d.0.04')
AllJTKresult.filtered<-filter(AllJTKresult,ADJ.P<0.05,PER>=20,PER<=28)
#filter for genes oscillate in at least 2 people
filtered.genes<-AllJTKresult.filtered[c("CycID","celltype")] %>% table() %>% as.data.frame() %>% dplyr::filter(.,Freq>=2)
filtered.genes<-filtered.genes$CycID %>% unique()
AllJTKresult.filtered<-AllJTKresult.filtered[AllJTKresult.filtered$CycID %in% filtered.genes,]
(AllJTKresult.filtered[c("CycID","celltype")] %>% unique())$celltype %>% table() %>% as.data.frame()

plotMetaCellByIndividual(srt.metacell,cell.type = "CD14 Mono",feature = "RPS27",norm.dist = T)
srt<-readRDS('7mixed.integrated.annotated.clean.sct.Azimuth.rds')

# scrutinize batch 
srt@meta.data[is.na(srt@meta.data$sample),"sample"]="X5_250424MIX03"
saveRDS(srt,'7mixed.integrated.annotated.clean.sct.Azimuth.rds')

srt$type<-srt$predicted.celltype.l2
matBySample<-srt2bulkMatrix(srt.metacell,c("individual","CT","type"),layer = "data",normalize = T)
matBySample<-rownames_to_column(matBySample,"feature")
matBySample.nr<-matBySample %>% gather(key="key",value="count",-feature)
matBySample.nr$individual<-getField(matBySample.nr$key,"-",1)
matBySample.nr$CT<-getField(matBySample.nr$key,"-",2)
matBySample.nr$type<-getField(matBySample.nr$key,"-",3)
matBySample.nr$key<-getField(matBySample.nr$key,"-",c(1,2))
matBySample.nr$batch<-BATCH[matBySample.nr$key]
matBySample.nr[is.na(matBySample.nr$count),"count"]<-0
saveRDS(matBySample.nr,'matBySample.rds')

calculateCorBatch("CD4_TCM","LYZ",matBySample.nr)
batch.by.individual.type=NULL
calculateCorBatch<-function(arg.type,arg.feature,mat){
  this.data=dplyr::filter(mat,feature==arg.feature,type==arg.type)[c("batch","individual","count")] %>% spread(key="individual",value="count")
  this.data=column_to_rownames(this.data,"batch")
  print(cor(this.data))
  this.cor=cor(this.data)[upper.tri(cor(this.data), diag = FALSE)]
  this.batch.data=data.frame("feature"=rep(arg.feature,length(this.cor)),"type"=rep(arg.type,length(this.cor)),"cor"=this.cor)
  return(this.batch.data)
}

srt.metacell$type<-srt.metacell$predicted.celltype.l2.main
classic_oscillating<-AllJTKresult.filtered[AllJTKresult.filtered$CycID %in% CIRCADIAN_GENES_MAIN,]
ggplot(classic_oscillating)+geom_point(aes(x=CycID,y=LAG,color=celltype,shape=individual))

classic_oscillating$cor_to_batch<-NA
for(i in 1:nrow(classic_oscillating)){
  classic_oscillating$cor_to_batch[i]=median(calculateCorBatch(gsub(" ","_",classic_oscillating$celltype[i]),
                                                        classic_oscillating$CycID[i],matBySample.nr)$cor)
}

srt.drop.list<-list()

srt.drop<-Read10X('/tmpdata/LyuLin/tmp/tmp/X5_250424MIX07/raw_feature_bc_matrix/')
srt.drop <- CreateSeuratObject(srt.drop, min.cells = 3, min.features = 0)
ggplot(srt.drop@meta.data)+geom_density(aes(x=nFeature_RNA))+xlim(c(6,500))+
  geom_vline(xintercept = 50,linetype=2)+geom_vline(xintercept = 175,linetype=2)+
  annotate("text", x = 50, y = Inf, label = "50", vjust = 2, size = 3, color = "red",hjust=-0.5) +
  annotate("text", x = 150, y = Inf, label = "175", vjust = 2, size = 3, color = "red",hjust=-0.5)+
  ylab("density")+ggtitle("sample: X5_250424MIX07\nnFeature of empty drop")
srt.drop.list<-c(srt.drop.list,subset(srt.drop,nFeature_RNA>50&nFeature_RNA<=175))
rm(srt.drop)
gc()

saveRDS(srt.drop.list,"/tmpdata/LyuLin/analysis/circadian/R/empty.drop.rds")
