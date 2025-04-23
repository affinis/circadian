library(Seurat)
library(dorothea)
library(viper)


test<-readRDS("~/data/analysis/all.samples.oddsremoved.tsneran.Azimuth.healthonly.rds")
#test<-readRDS('~/data/analysis/TF_SLY_1.levelTop.rds')
#run viper to show tf activities
small.test<-test
#small.test<-subset(test,patient=="TFSH190500A_HZD")
expression_data <- as.matrix(LayerData(small.test,layer = "data"))
# Load human regulons (for example)
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs
# Run VIPER
tf_activities <- run_viper(expression_data, regulons, 
                           options = list(method = "scale", minsize = 4, 
                           eset.filter = FALSE, cores = 95, verbose = FALSE))
# Add TF activities to Seurat object
small.test[["TF"]] <- CreateAssayObject(counts = tf_activities)
saveRDS(small.test,"~/data/analysis/all.samples.oddsremoved.tsneran.Azimuth.healthonly.TF.rds")
