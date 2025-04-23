library(tauFisher)
library(tidyverse)

getField<-function(vec,sep,field){
  if(length(field)>1){
    strsplit(vec,sep) %>% lapply(.,`[`,field) %>% lapply(.,paste0,collapse=sep) %>% unlist()
  }else{
    strsplit(vec,sep) %>% lapply(.,`[`,field) %>% unlist()
  }
}

#using mm data, load bulk data
bulk_file <- system.file("extdata", "GSE157077_mouse_scn_control.tsv", 
                         package = "tauFisher", mustWork = TRUE)
bulk_data <- utils::read.delim(file = bulk_file, stringsAsFactors = FALSE)

#using hs data
bulk_data_hs<-read.delim('/lustre/home/acct-medll/medll/data/bulk_RNA-seq/GSE113883_TPM_processed_data.txt')
bulk_data_hs$GENE<-getField(bulk_data_hs$GENE,"__",1)
colnames(bulk_data_hs)[1]<-"ID"
colnames(bulk_data_hs)[2:length(colnames(bulk_data_hs))]<-getField(colnames(bulk_data_hs)[2:length(colnames(bulk_data_hs))],"_",c(6,1,2))
reps<-colnames(bulk_data_hs)[2:length(colnames(bulk_data_hs))] %>% getField(.,"_",2:3) %>% unique() 
rep_name<-paste0("REP_",1:length(reps))
i=1
for(i in 1:length(reps)){
  if(i %in% c(4,6,11)){
    next
  }
  colnames(bulk_data_hs)=gsub(reps[i],rep_name[i],colnames(bulk_data_hs))
}
colnames(bulk_data_hs)<-gsub("^X","CT_",colnames(bulk_data_hs))
colnames(bulk_data_hs)<-gsub("hr","",colnames(bulk_data_hs))
colnames(bulk_data_hs)<-gsub("_0","_",colnames(bulk_data_hs))
bulk_data_hs<-bulk_data_hs[,!(grepl("CT_27|CT_29|CT_25",colnames(bulk_data_hs)))]
bulk_data_hs<-bulk_data_hs[,grepl("REP|ID",colnames(bulk_data_hs))]
colnames(bulk_data_hs)<-gsub("REP_10","REP_4",colnames(bulk_data_hs))
colnames(bulk_data_hs)<-gsub("REP_9","REP_6",colnames(bulk_data_hs))
bulk_data<-bulk_data_hs


bulk_data_hs_mod<-bulk_data_hs
bulk_data_hs_mod<-gather(bulk_data_hs_mod,key=sample,value=expression,-GENE)
bulk_data_hs_mod$CT<-getField(bulk_data_hs_mod$sample,"_",6)
bulk_data_hs_mod[grepl("^PER1",bulk_data_hs_mod$GENE),] %>% ggplot(.)+geom_point(aes(x=CT,y=expression))


# Take the average expression if there are non-unique genes
bulk<-bulk_data[stats::complete.cases(bulk_data[,-1]), ] %>%
  dplyr::mutate(ID = toupper(ID)) %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise_all(mean) %>%
  as.data.frame()

#Set the rownames to be the genes and remove the ID column 
rownames(bulk)<-bulk$ID
bulk<-bulk[,-1]

#parse the column names to obtain a vector of time points
time<-as.numeric(vapply(stringr::str_split(colnames(bulk), "_"),
                          '[', 2, FUN.VALUE = character(1)) )
#replicate<-as.numeric(vapply(stringr::str_split(colnames(bulk), "_"),
#                              '[', 4, FUN.VALUE = character(1)) )
cycle<-0
time_adj<-NULL
for(i in 1:length(time)){
  if(is.null(time_adj)){
    time_adj=time[i]
    next
  }
  if(time[i]<time[i-1]){
    cycle=cycle+1
  }
  time_adj=c(time_adj,time[i]+cycle*24)
}

# adjust time - each replicate is now the next 'period' so there's 3 sets
#time_adj <- time + 24*(as.numeric(replicate) - 1)

bulk_adj = bulk
colnames(bulk_adj)<-time_adj
bulk_adj = bulk_adj[, order(as.numeric(colnames(bulk_adj)))]
time_adj <- time_adj[order(as.numeric(time_adj))]

#using mm pseudobulk, load pseudo bulk data from single cell transcriptome
#pseudobulk_file <- system.file("extdata", "GSE132608_scn_pseudobulk.RData", 
#                               package = "tauFisher", mustWork = TRUE)
#load(pseudobulk_file)


#using human pseudobulk data
srts<-getSamplesOneDayToList(HEALTH[5])
pseudobulk_data_hs<-NULL
for(i in 1:length(srts)){
  this.data=LayerData(srts[[i]], assay = "RNA", layer = "counts") %>% rowSums() %>% as.data.frame()
  colnames(this.data)=CT_TIME[i]
  this.data=rownames_to_column(this.data,var="ID")
  if(is.null(pseudobulk_data_hs)){
    pseudobulk_data_hs=this.data
    next
  }else{
    pseudobulk_data_hs=left_join(pseudobulk_data_hs,this.data,by="ID")
  }
}
pseudobulk_data_hs[is.na(pseudobulk_data_hs)]<-0
rownames(pseudobulk_data_hs)<-pseudobulk_data_hs$ID
pseudobulk_data_hs$ID<-NULL
pseudobulk_time<-as.numeric(vapply(stringr::str_split(colnames(pseudobulk_data_hs), "CT"),
                                   '[', 2, FUN.VALUE = character(1)))

#
#
#run_metacycle(df=bulk_adj, 
#              timepoints=time_adj, 
#              out_file="/lustre/home/acct-medll/medll/data/analysis/tauFisher_test/GSE113883_hs_full_adj.txt", 
#              out_dir="/lustre/home/acct-medll/medll/data/analysis/tauFisher_test/GSE113883_hs_full_adj_meta2d", 
#              method=c("JTK","LS"))

meta2d_file = "/lustre/home/acct-medll/medll/data/analysis/tauFisher_test/GSE113883_hs_full_adj_meta2d/meta2d_GSE113883_hs_full_adj.txt"
genes = find_periodic_genes(input_file=meta2d_file,
                            test_genes=rownames(pseudobulk_data_hs),
                            per_method=c("LS","JTK"),
                            thres=1)

print(genes)

# Capitalize all gene names to standardize format
rownames(bulk_adj) <- toupper(rownames(bulk_adj))
# Subset  
chosen_genes <- genes$JTK[genes$JTK %in% rownames(bulk_adj)]
bulk_subset <- t(data.frame(bulk_adj[chosen_genes, ]))
rownames(bulk_subset) <- time_adj

# Order the genes 
chosen_genes <- chosen_genes[order(chosen_genes)]
bulk_subset <- bulk_subset[, order(colnames(bulk_subset))]

# log2 transform
bulk_log <- log2(bulk_subset+1)

#run FDA
nrep = 8 # number of replicates
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
fda_diff <- create_DiffMatrix(genes=chosen_genes, dat=fda_expression2)
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

#Run multinomial regression
ndims = 2
pc_data<-data.frame(pc$x[,1:ndims])

# Get the times and relevel them so the smallest CT is the reference level
pc_data$CT24<-as.numeric(vapply(stringr::str_split(rownames(pc_data), pattern = '_'),
                                  '[', 1, FUN.VALUE = character(1) ))
pc_data$CT24_relevel<-stats::relevel(factor(pc_data$CT24),
                      ref = as.character(min(train_time)))

mod <- nnet::multinom(CT24_relevel ~ PC1 + PC2, data = pc_data, trace=T)

# Subset the data on the chosen genes
rownames(pseudobulk_data_hs) <- toupper(rownames(pseudobulk_data_hs))
pseudobulk_subset <- t(data.frame(pseudobulk_data_hs[chosen_genes, ]))

# Order the genes
pseudobulk_subset <- pseudobulk_subset[, order(colnames(pseudobulk_subset))]

# log2 transform
pseudobulk_log <- log2(pseudobulk_subset+1)

#Calculate differences for each gene pair
pseudobulk_diff <- create_DiffMatrix(genes=chosen_genes, dat=pseudobulk_log)
pseudobulk_diff_scaled = scale_DiffMatrix(diffs=pseudobulk_diff)

# Project data onto PCA space
pc_pred <- stats::predict(pc, newdata = pseudobulk_diff_scaled)

# Predict
pred_vals = stats::predict(mod, newdata = pc_pred[,1:ndims])
pred_vals = as.numeric(as.character(pred_vals))

calc_error(truth=pseudobulk_time, pred=pred_vals)
