source('/tmpdata/LyuLin/script/circadian/circadian_core.R')
setwd('/tmpdata/LyuLin/analysis/circadian/R/preparation_seacell')

time_points<-c("CT11","CT15","CT19","CT23","CT27","CT31","CT35")
individuals<-c("ZYR","JJC","KD","LYH")

srt<-readRDS('../7mixed.integrated.annotated.clean.sct.Azimuth.rds')

for(time_point in time_points){
  for(this.individual in individuals){
    this.srt=subset(srt,CT==time_point&individual==this.individual)
    saveRDS(this.srt,paste0(this.individual,"_",time_point,".rds"))
  }
}

