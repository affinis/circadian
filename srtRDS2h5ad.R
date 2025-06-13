# devtools::install_github("cellgeni/sceasy")
# BiocManager::install(c("LoomExperiment", "SingleCellExperiment"))
# conda install anndata -c bioconda

sceasy::convertFormat(seurat_object, from="seurat", to="anndata",
                      outFile='filename.h5ad')
