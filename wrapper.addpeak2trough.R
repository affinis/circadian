#source('/tmpdata/LyuLin/script/circadian/circadian_core.R')
source('~/script/circadian/circadian_core.R')

# upstream: wrapper.compareExpressionWithDroplet.R
# downstream: 
# dependency: 
# caller: NSF
##
# speed: ~30 min

### need modify every run
#srt.metacell<-readRDS('~/analysis/circadian/R/seacell.16individual.rds')
#srt.metacell<-readRDS('~/analysis/circadian/R/CYTOF.300k.reannotation.srt.rds')
srt.metacell<-readRDS('~/analysis/circadian/R/16individual.unspliced.srt.rds')

#srt.metacell$type<-srt.metacell$manual.new

### need modify every run
#Allresult.filtered<-rbind(readRDS('~/analysis/circadian/R/JTK.12individual.raw.rds'),readRDS('~/analysis/circadian/R/JTK.4individual.raw.rds'))
#Allresult.filtered<-dplyr::filter(Allresult.filtered,ADJ.P<0.05)
Allresult.filtered<-readRDS('~/analysis/circadian/R/JTK.unspliced.filtered.16individual.addmedianexp.rds')
Genes<-Allresult.filtered$CycID %>% unique()
srt.metacell<-NormalizeData(srt.metacell,normalization.method = "RC",scale.factor=1000000)
exp_mat<-LayerData(srt.metacell,features=Genes,layer = "data") %>% as.data.frame()
exp_mat<-rownames_to_column(exp_mat,var = "feature")
exp_mat<-gather(exp_mat,key="observation",value="expression",-feature)
exp_mat$individual<-getField(exp_mat$observation,"_",1)
exp_mat$CT<-getField(exp_mat$observation,"_",2)
exp_mat$celltype<-getField(exp_mat$observation,"_",3)

#following 4 lines only work for cytof data to replace 'exp_mat$celltype<-getField(exp_mat$observation,"_",4) %>% gsub("-"," ",.)'
#annotation<-srt.metacell$type %>% as.data.frame()
#colnames(annotation)<-"celltype"
#annotation<-rownames_to_column(annotation,"observation")
#exp_mat<-left_join(exp_mat,annotation,by="observation")

mean_exp_mat<-exp_mat %>% group_by(feature,individual,CT,celltype) %>% summarise(mean.expression=mean(expression))
mean_exp_mat<-mean_exp_mat %>% group_by(feature,individual,celltype) %>% mutate(max.expression=max(mean.expression),min.expression=min(mean.expression))
mean_exp_mat$fold_peak2trough<-mean_exp_mat$max.expression/(1+mean_exp_mat$min.expression)
# following 1 line for cytof
#mean_exp_mat$fold_peak2trough<-mean_exp_mat$max.expression/(mean_exp_mat$min.expression)

mean_exp_mat<-mean_exp_mat[c("feature","individual","celltype","fold_peak2trough")] %>% unique()
colnames(mean_exp_mat)[1]<-"CycID"
mean_exp_mat$celltype<-mean_exp_mat$celltype %>% gsub("_"," ",.)
print(head(mean_exp_mat))
print(head(Allresult.filtered))

Allresult.filtered<-left_join(Allresult.filtered,mean_exp_mat,by=c("CycID","individual","celltype"))

### need modify every run
saveRDS(Allresult.filtered,"~/analysis/circadian/R/JTK.unspliced.filtered.16individual.addmedianexp.addp2t.rds")
