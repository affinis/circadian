# This is a script for generating seacell (a python tool for calculating metacells) inputs with cellranger multi output with a pre-generated cellannotation file

# upstream: cellranger, circadian_core.R:generateAnnotationFile
# downstream: <bash>runSeaCells.sh, <slurm>runSeaCells.slurm
# dependency: DropletUtils (install.packages("DropletUtils"))
# caller: NSF
##
# speed (n: number of smallest loop): n=340, t~=60 min; n=520, t~=80 min

#source('~/script/circadian/circadian_core.R')
library('DropletUtils')
library('Seurat')
library('tidyverse')

main.dir='/dssg/home/acct-medll/medll/analysis/cellranger'
annotation.file='/dssg/home/acct-medll/medll/analysis/circadian/R/cell.annotations.manual1_2_NI.predicted2.modforSeacells.tsv'
setwd(main.dir)

time_points<-c("CT09","CT13","CT17","CT21","CT25","CT29")
individuals<-c("HZD","LJQ","SLY","XSP","YXQ","ZF","ZXL","ZYJ","ZYX")

# specify input
#srt<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/7mixed.integrated.annotated.clean.sct.Azimuth.rds')
annotation<-read.delim(annotation.file)
#types<-unique(annotation$predicted.celltype.l2)

for(time_point in time_points){
  for(this.individual in individuals){
    out.dir=file.path(main.dir,paste0(this.individual,"_",time_point),
                      "outs/per_sample_outs",
                      paste0(this.individual,"_",time_point),"count/preparation_seacells")
    CR.out.dir=file.path(main.dir,paste0(this.individual,"_",time_point),
                         "outs/per_sample_outs",
                         paste0(this.individual,"_",time_point),"count/sample_filtered_feature_bc_matrix")
    if(!dir.exists(CR.out.dir)){
      next
    }
    if(!dir.exists(out.dir)){
      dir.create(out.dir)
    }
    if(!file.exists(file.path(CR.out.dir,"barcodes.tsv.gz"))){
      system(paste0('gzip ',file.path(CR.out.dir,"barcodes.tsv")))
      system(paste0('cp ',file.path(CR.out.dir,"barcodes.tsv.gz"),' ',file.path(CR.out.dir,"barcodes.tsv.gz.backup")))
      barcode.1=(gzfile(file.path(CR.out.dir,"barcodes.tsv.gz")) %>% scan(.,what="char"))[1]
      if(!startsWith(barcode.1,"TF_")){
        system(paste0('zcat ',file.path(CR.out.dir,"barcodes.tsv.gz")," | sed 's/^/TF_",this.individual,"_",time_point,"_/g' > ",file.path(CR.out.dir,"barcodes.tsv")))
        system(paste0('gzip -f ',file.path(CR.out.dir,"barcodes.tsv")))
      }
    }
    this.all.srt=CreateSeuratObject(counts=Read10X(CR.out.dir))
    if(!startsWith(Cells(this.all.srt)[1],"TF_")){
      this.all.srt=RenameCells(this.all.srt,add.cell.id=paste0("TF_",this.individual,"_",time_point))
    }
    print(head(Cells(this.all.srt)))
    valid_cell_id=intersect(Cells(this.all.srt),annotation$cell_id)
    if(length(valid_cell_id)==0){
      next
    }
    this.all.srt=subset(this.all.srt,cells=valid_cell_id)
    this.annotation=annotation[annotation$cell_id %in% valid_cell_id,]
    print(head(this.annotation))
    rownames(this.annotation)=this.annotation$cell_id
    #this.annotation=column_to_rownames(this.annotation,var="cell_id")
    this.all.srt=AddMetaData(this.all.srt,this.annotation)
    for (this.type in unique(this.annotation$predicted.celltype.l2)) {
      message(paste0(time_point,"_",this.individual,"_",this.type))
      this.srt=subset(this.all.srt,subset=predicted.celltype.l2==this.type)
      print(head(Cells(this.srt)))
      if(length(Cells(this.srt))==1){
        next
      }
      # Assuming your Seurat object is named 'seurat_obj'
      # Extract the raw count matrix (use counts, not normalized data)
      count_matrix=GetAssayData(object = this.srt, assay = "RNA", layer = "counts")
      
      # Get cell barcodes
      barcodes=colnames(this.srt)
      
      # Get gene features
      features=data.frame(
        gene_id = rownames(this.srt),
        gene_name = rownames(this.srt),
        feature_type = "Gene Expression"
      )
      
      # Create output directory
      output_dir=paste0(this.individual,"_",time_point,"_",gsub(" ","-",this.type),"_filtered_bc_matrix")
      output_dir=file.path(out.dir,output_dir)
      dir.create(output_dir, showWarnings = FALSE)
      
      # Write files in Cell Ranger format
      write10xCounts(
        path = output_dir,
        x = count_matrix,
        barcodes = barcodes,
        gene.id = features$gene_id,
        gene.symbol = features$gene_name,
        gene.type = features$feature_type,
        overwrite = TRUE,version = "3",
        chemistry = "Single Cell 5' v3"
      )
    }
    # Gzip the files (optional if write10xCounts didn't do it automatically)
    #system(paste("gzip", file.path(output_dir, "barcodes.tsv")))
    #system(paste("gzip", file.path(output_dir, "features.tsv")))
    #system(paste("gzip", file.path(output_dir, "matrix.mtx")))
  }
}

