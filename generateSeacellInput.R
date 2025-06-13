# This is a script for generating seacell (a python tool for calculating metacells) inputs with seurat objects
# upstream: <Azimuth>RunAzimuth
# downstream: <bash>runSeaCells.sh
# dependency: DropletUtils (install.packages("DropletUtils"))
# caller: NSF

source('/tmpdata/LyuLin/script/circadian/circadian_core.R')
library('DropletUtils')
setwd('/tmpdata/LyuLin/analysis/circadian/R/preparation_seacell')

time_points<-c("CT11","CT15","CT19","CT23","CT27","CT31","CT35")
individuals<-c("ZYR","JJC","KD","LYH")

# specify input
srt<-readRDS('../7mixed.integrated.annotated.clean.sct.Azimuth.rds')

for(time_point in time_points){
  for(this.individual in individuals){
    this.srt=subset(srt,CT==time_point&individual==this.individual)
    
    # Assuming your Seurat object is named 'seurat_obj'
    # Extract the raw count matrix (use counts, not normalized data)
    count_matrix=GetAssayData(object = this.srt, assay = "RNA", slot = "counts")
    
    # Get cell barcodes
    barcodes=colnames(this.srt)
    
    # Get gene features
    features=data.frame(
      gene_id = rownames(this.srt),
      gene_name = rownames(this.srt),
      feature_type = "Gene Expression"
    )
    
    # Create output directory
    output_dir=paste0(this.individual,"_",time_point,"_filtered_bc_matrix")
    dir.create(output_dir, showWarnings = FALSE)
    
    # Write files in Cell Ranger format
    write10xCounts(
      path = output_dir,
      x = count_matrix,
      barcodes = barcodes,
      gene.id = features$gene_id,
      gene.symbol = features$gene_name,
      gene.type = features$feature_type,
      overwrite = TRUE
    )
    
    # Gzip the files (optional if write10xCounts didn't do it automatically)
    system(paste("gzip", file.path(output_dir, "barcodes.tsv")))
    system(paste("gzip", file.path(output_dir, "features.tsv")))
    system(paste("gzip", file.path(output_dir, "matrix.mtx")))
  }
}

