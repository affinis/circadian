# This is a script for generating seacell (a python tool for calculating metacells) inputs with seurat objects
# upstream: 
# downstream: <bash>runSeaCells.sh
# dependency: DropletUtils (install.packages("DropletUtils"))
# caller: NSF

source('/tmpdata/LyuLin/script/circadian/circadian_core.R')
library('DropletUtils')
setwd('/tmpdata/LyuLin/analysis/circadian/R/preparation_seacell_EmptyDrop')

srt.list<-readRDS("/tmpdata/LyuLin/analysis/circadian/R/empty.drop.rds")
samples<-c("X5_250409MIX01","X5_250424MIX02","X5_250424MIX03","X5_250424MIX04","X5_250424MIX05","X5_250424MIX06","X5_250424MIX07")
all_cell_ids=NULL
for (i in 1:length(samples)) {
  # Assuming your Seurat object is named 'seurat_obj'
  # Extract the raw count matrix (use counts, not normalized data)
  this.srt=srt.list[[i]]
  new_cell_ids <- paste0(samples[i],"_",colnames(this.srt))
  this.srt=RenameCells(this.srt,new.names = new_cell_ids)
  sampled_ids=sample(new_cell_ids,10000,F)
  this.srt=subset(this.srt,cells=sampled_ids)
  all_cell_ids=c(all_cell_ids,sampled_ids)
  
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
  output_dir=paste0(samples[i],"_filtered_bc_matrix")
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
  )
  
  # Gzip the files (optional if write10xCounts didn't do it automatically)
  system(paste("gzip", file.path(output_dir, "barcodes.tsv")))
  system(paste("gzip", file.path(output_dir, "features.tsv")))
  system(paste("gzip", file.path(output_dir, "matrix.mtx")))
}

meta.data=data.frame('cell_id'=all_cell_ids,
                     'manual.level2'=rep("empty_droplet",length(all_cell_ids)),
                     'predicted.celltype.l2'=rep("empty_droplet",length(all_cell_ids)))
write_delim(meta.data,'/tmpdata/LyuLin/analysis/circadian/R/droplet.annotation.tsv',delim = '\t')
