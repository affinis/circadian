library(Seurat)

srt<-readRDS('/lustre/home/acct-medll/medll/data/analysis/all.samples.oddsremoved.tsneran.Azimuth.healthonly.rds')

srt<-SCTransform(srt,return.only.var.genes=FALSE,ncells=20000)

srt@assays$SCT@scale.data<-matrix()

saveRDS(srt,'/lustre/home/acct-medll/medll/data/analysis/all.samples.oddsremoved.tsneran.Azimuth.healthonly.SCTransform.rds')
