source('/tmpdata/LyuLin/script/circadian/circadian_core.R')

# upstream: wrapper.compareExpressionWithDroplet.R
# downstream: 
# dependency: 
# caller: NSF
##
# speed: ~30 min

srt.metacell<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/seacell.0.08.bytype.rds')
Allresult.filtered<-readRDS('/tmpdata/LyuLin/analysis/circadian/R/LS.result.filtered.addp2bkg.addcor2batch.bytype.0.08.rds')
Genes<-Allresult.filtered$CycID %>% unique()
srt.metacell<-NormalizeData(srt.metacell,normalization.method = "RC",scale.factor=1000000)
exp_mat<-LayerData(srt.metacell,features=Genes,layer = "data") %>% as.data.frame()
exp_mat<-rownames_to_column(exp_mat,var = "feature")
exp_mat<-gather(exp_mat,key="observation",value="expression",-feature)
exp_mat$individual<-getField(exp_mat$observation,"_",2)
exp_mat$CT<-getField(exp_mat$observation,"_",3)
exp_mat$celltype<-getField(exp_mat$observation,"_",4) %>% gsub("-"," ",.)
mean_exp_mat<-exp_mat %>% group_by(feature,individual,CT,celltype) %>% summarise(mean.expression=mean(expression))
mean_exp_mat<-mean_exp_mat %>% group_by(feature,individual,celltype) %>% mutate(max.expression=max(mean.expression),min.expression=min(mean.expression))
mean_exp_mat$fold_peak2trough<-mean_exp_mat$max.expression/(1+mean_exp_mat$min.expression)
mean_exp_mat<-mean_exp_mat[c("feature","individual","celltype","fold_peak2trough")] %>% unique()
colnames(mean_exp_mat)[1]<-"CycID"

Allresult.filtered<-left_join(Allresult.filtered,mean_exp_mat,by=c("CycID","individual","celltype"))

saveRDS(Allresult.filtered,"/tmpdata/LyuLin/analysis/circadian/R/LS.result.filtered.addp2bkg.addcor2batch.addp2t.bytype.0.08.rds")
