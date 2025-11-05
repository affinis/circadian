# This is a script for generating seacell (a python tool for calculating metacells) inputs with seurat objects

# upstream: <Azimuth>RunAzimuth
# downstream: <bash>runSeaCells.sh, <slurm>runSeaCells.slurm
# dependency: DropletUtils (install.packages("DropletUtils"))
# caller: NSF
##
# speed (n: number of smallest loop): n=340, t~=60 min; n=520, t~=80 min

source('~/script/circadian/circadian_core.R')
library('DropletUtils')

out.dir="~/analysis/circadian/R/preparation_seacell_PseudoCellranger"
if(!dir.exists(out.dir)){
  dir.create(out.dir)
}
setwd(out.dir)

time_points<-c("CT11","CT15","CT19","CT23","CT27","CT31","CT35")
individuals<-c("ZYR","JJC","KD","LYH")

# specify input
srt<-readRDS('~/analysis/circadian/R/16individual.srt.annotated.rds')
srt$type<-srt$predicted.celltype.l2
types<-unique(srt$type)

for(time_point in time_points){
  for(this.individual in individuals){
    for (this.type in types) {
      message(paste0(time_point,this.individual,this.type,collapse = "_"))
      n.cell=nrow(srt@meta.data[srt@meta.data$CT==time_point&srt@meta.data$individual==this.individual&srt@meta.data$type==this.type,])
      if(n.cell<25){
        next
      }
      this.srt=subset(srt,CT==time_point&individual==this.individual&type==this.type)
      
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
      pseudo_cellranger_outdir=paste0(this.individual,"_",time_point,"/outs/per_sample_outs/",this.individual,"_",time_point,"/count/preparation_seacells/")
      output_dir=paste0(this.individual,"_",time_point,"_",gsub(" ","-",this.type),"_filtered_bc_matrix")
      output_dir=paste0(pseudo_cellranger_outdir,output_dir)
      dir.create(output_dir, showWarnings = FALSE,recursive = T)
      
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

