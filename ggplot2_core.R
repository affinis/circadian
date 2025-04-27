library("ggplot2")
#library('ggpubr')
library("ggsci")
library("forcats")
library("scales")
library("gggenes")
library("tidyr")

TransformDiscreteDf<-function(data){
  cns=colnames(data)
  rns=rownames(data)
  ret=data.frame("observations"=character(0),"features"=character(0),"groups"=character(0))
  for (i in cns){
#    print(i)
    for (j in rns){
      value=data[j,i]
      subdf=data.frame(observations=j,features=i,groups=value)
      ret=rbind(ret,subdf)
    }
  }
  return(ret)
}

#transform a numeric matrix into 'dds' form, rows were defined to observations, cols were defined to features.
TransformContinuousDf<-function(df,minors.threshold=0){
  data=t(df)
  cns=colnames(data)
  rns=rownames(data)
  ret=data.frame("observations"=character(0),"features"=character(0),"values"=integer(0))
  if(is.null(minors.threshold)){
    for(i in cns){
      print(i)
      subdf=data.frame("observations"=rep(i,length(rns)),"features"=rns,"values"=data[,i] %>% as.vector())
      ret=rbind(ret,subdf)
    }
    return(ret)
  }
  for (i in cns){
    print(i)
    if(sum(data[,i])==0){
      next
    }
    features=(rownames(data))[data[,i]/sum(data[,i]) >= minors.threshold]
    if (length(features)==0){
      #      features=genus_count[,i] %>% which.max() %>% names()
      subdf=data.frame(observations=i,features="others",values=sum(data[,i]))
      ret=rbind(ret,subdf)
      next
    }
    values=data[features,i] %>% as.integer()
    observations=rep(i,length(features))
    subdf=data.frame(observations=observations,features=features,values=values)
    ret=rbind(ret,subdf)
    if (sum(data[,i])-sum(subdf$values)>0){
      otherline=data.frame(observations=i,features="others",values=sum(data[,i])-sum(subdf$values))
      ret=rbind(ret,otherline)
    }
  }
  return(ret)
}

#input: data frame with 1 column filled with discrete data
dplot<-function(data){
  ggplot(data,aes(x=get(colnames(data))))+
    geom_bar(fill="blue",width=0.5)+
    xlab(colnames(data))+
    theme(axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.background = element_rect(fill="white"))
}

#input: data frame with 1 column filled with continuous data
cplot<-function(data,type="density"){
  p=ggplot(data,aes(x=get(colnames(data))))
  if(type=="density"){
    p=p+geom_density()
  }else if(type=="hist"){
    p=p+geom_histogram(binwidth=1)
  }else if(type=="qq"){
    p=p+geom_qq()
  }else{
    stop("Unsupported type")
  }
  p+xlab(colnames(data))+
    theme(axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.background = element_rect(fill="white"))
}

#input: data frame with 2 columns filled with continuous data
ccplot<-function(data,type="dot"){
  coln=colnames(data)
  p=ggplot(data,aes(x=get(coln[1]),y=get(coln[2])))
  if(type=="dot"){
    p=p+geom_point()
  }else if(type=="line"){
    p=p+geom_line()
  }else{
    stop("Unsupported type")
  }
  p+xlab(coln[1])+ylab(coln[2])+
    theme(axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.background = element_rect(fill="white"))
}

dcplot<-function(data,type="bar"){
  coln=colnames(data)
  p=ggplot(data,aes(x=get(coln[1]),y=coln[2]))
  if(type=="bar"){
    p=p+geom_col()
  }else if(type=="box"){
    p=p+geom_boxplot()
  }else if(type=="violin"){
    p=p+geom_violin()
  }else{
    stop("Unsupported type")
  }
  p+xlab(coln[1])+ylab(coln[2])+
    theme(axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.background = element_rect(fill="white"))
}

ddplot<-function(data,type="bar",reverse.col=F,x.order=NULL,font.size=13){
  coln=colnames(data)
  data=data %>% group_by(.,get(coln[1]),get(coln[2])) %>% summarise(.,count=n())
  colnames(data)[1:2]=coln
  if(reverse.col){
    coln=rev(coln)
  }
  data=as.data.frame(data)
  if(!is.null(x.order)){
    data[,coln[1]]<-factor(data[,coln[1]],levels=x.order)
  }
  print(head(data))
  if(type=="bar"){
    p=ggplot(data,aes(x=get(coln[1]),y=count,fill=get(coln[2])))+
      geom_bar(stat="identity",position="fill",color="black",width=0.75)+theme_publication(base.size=font.size)+
      ylab("ratio")+xlab("")+scale_fill_jama()+scale_y_continuous(expand=c(0,0))+labs(fill="")
    p
  }else{
    stop("unsupported type")
  }
}

#input: data frame with 3 columns, first two are continuous data and the last column is discrete data
ccdplot<-function(data,width=0.5,label=F,label.size=4.5,hjust=0.05,position.dodge=1){
  if(label){
    y.expand.upper=hjust+0.05
  }else{
    y.expand.upper=0.05
  }
  coln=colnames(data)
  n.row=nrow(data)
  data=gather(data,feature,value,-coln[3])
  coln=colnames(data)
  p=ggplot(data)+geom_bar(aes(x=fct_inorder(get(coln[1])),y=get(coln[3]),group=get(coln[2]),fill=get(coln[2])),
                        stat="identity",position="dodge",color="black",width=width)+
    xlab(coln[1])+ylab(coln[3])+guides(fill=guide_legend(title=coln[2]))+scale_fill_manual(values=c("red","blue"))+
    scale_y_continuous(expand = c(0,0),limits = c(0,y.expand.upper+ceiling(max(data[3]))))
  if(label){
    p=p+geom_text(aes(x=fct_inorder(get(coln[1])),y=get(coln[3])+hjust,group=get(coln[2]),label=get(coln[3])),size=label.size,position=position_dodge(position.dodge))
  }
  p
}


#input: data frame with 3 columns, all of the column store discrete data
dddplot<-function(data,flip_coord=F){
  coln=colnames(data)
  n.row=nrow(data)
  if(flip_coord){
    ggplot(data)+geom_tile(aes(x=get(coln[1]),y=get(coln[2]),fill=get(coln[3])))+coord_flip()+
      xlab(coln[1])+ylab(coln[2])+guides(fill=guide_legend(title=coln[3]))+
      theme(axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.background = element_rect(fill="white"))
  }else{
    ggplot(data)+geom_tile(aes(x=get(coln[1]),y=get(coln[2]),fill=get(coln[3])))+
      xlab(coln[1])+ylab(coln[2])+guides(fill=guide_legend(title=coln[3]))+
      theme(axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.background = element_rect(fill="white"))
  }
}

#input: data frame with 3 columns, discrete, discrete, continuous data respectively
ddcplot<-function(data,colors=c("blue","red"),ob.order=NULL,feat.order=NULL,type="heatmap"){
  coln=colnames(data)
  n.row=nrow(data)
  if(is.null(ob.order)){
    data[,coln[1]]=factor(data[,coln[1]],levels=levels(fct_inorder(data[,coln[1]])))
  }else{
    data[,coln[1]]=factor(data[,coln[1]],levels=levels(fct_inorder(factor(ob.order))))
  }
  if(is.null(feat.order)){
    data[,coln[2]]=factor(data[,coln[2]],levels=levels(fct_inorder(data[,coln[2]])))
  }else{
    data[,coln[2]]=factor(data[,coln[2]],levels=levels(fct_inorder(factor(feat.order))))
  }
  if(type=="heatmap"){
    ggplot(data)+geom_tile(aes(x=get(coln[1]),y=get(coln[2]),fill=get(coln[3])))+coord_flip()+
      xlab(coln[1])+ylab(coln[2])+scale_fill_gradientn(colours=colors,guide="colourbar")+
      guides(fill=guide_colorbar(title=coln[3],keyheight=0.1,label=T,ncol=1))+
      theme(axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.text.x=element_text(angle=45,hjust=1),
            panel.background = element_rect(fill="white"))
  }else if(type=="bar"){
    ggplot(data)+geom_bar(aes(x=get(coln[1]),y=get(coln[3]),fill=get(coln[2])),stat="identity",position="fill")+
      xlab(coln[1])+ylab("ratio")+guides(fill=guide_legend(title=coln[2]))+scale_fill_simpsons("springfield")+
      theme(text=element_text(size=15),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.text.x=element_text(angle=45,hjust=1),
            panel.background = element_rect(fill="white"),
            legend.position = "bottom")
  }else{
    stop("Unsupported type")
  }
}

#input: data frame with n columns, all of the columns store discrete data, very suitable for metadata
#rows are observations, cols are features
ndplot<-function(data,ob.level=NULL,flip.coord=T){
  data=TransformDiscreteDf(data)
  if(!is.null(ob.level)){
    data[,1]=factor(data[,1],levels=ob.level)
  }else{
    data[,1]=factor(data[,1],levels=levels(fct_inorder(data[,1])))
  }
  if(flip.coord){
    p=dddplot(data,flip_coord=T)
  }else{
    p=dddplot(data)
  }
  if(length(data[,3]%>%unique())>10){
    p=p+scale_fill_gradientn(colours = c("black","white"))
  }else{
    p=p+scale_fill_d3()
  }
  p
}

#input: data frame with n columns, all of the columns store continuous data
#rows are observations, cols are features
ncplot<-function(data,heatmap.colors=c("black","white"),
  type="heatmap",cluster.by.row=F,cluster.by.col=F,
  minors.threshold=0.05,observations.sort=NULL,features.sort=NULL){
  if((cluster.by.row&!is.null(observations.sort))|(cluster.by.col&!is.null(features.sort))){
    stop("'cluster.by.row/col' is not compatible with sort.observations/features")
  }
  yorder=rownames(data)
  xorder=colnames(data)
  if(cluster.by.row){
    dy=dist(data,method = "euclidian")
    cy=hclust(dy,method = "ward.D2")
    if(!is.null(cy$labels)){
      yorder=cy$labels[cy$order]
    }else{
      yorder=cy$order
    }
  }
  if(cluster.by.col){
    dx=dist(t(data),method = "euclidian")
    cx=hclust(dx,method = "ward.D2")
    if(!is.null(cx$labels)){
      xorder=cx$labels[cx$order]
    }else{
      xorder=cx$order
    }
  }
  if(!is.null(observations.sort)){
    yorder=rev(observations.sort)
  }
  if(!is.null(features.sort)){
    xorder=features.sort
  }
  if(type=="heatmap"){
    data=TransformContinuousDf(data,minors.threshold=minors.threshold)
    if(!is.null(minors.threshold)){
      if(minors.threshold>0){
        xorder=append(xorder,"others")
      }
    }
    ddcplot(data,colors=heatmap.colors,ob.order=yorder,feat.order=xorder)
  }else if(type=="bar"){
    yorder=rev(yorder)
    data=TransformContinuousDf(data,minors.threshold=minors.threshold)
    if(!is.null(minors.threshold)){
      if(minors.threshold>0){
        xorder=append(xorder,"others")
      }
    }
    ddcplot(data,type="bar",ob.order=yorder,feat.order=xorder)
  }else{
    stop("Unsupported plot type")
  }
}

#if length(ds)==1, use exact name of it
ncmdplot.mod<-function(data,cs,ds,type="heatmap",show.obname=T,show.featname=T,logTrans=F,scale.by="none",cluster.by.y=FALSE,
                   cluster.by.x=FALSE,heatmap.main.colors=c("black","white"),legend.colors=pal_d3("category20")(20),
                   colorbar.width=NULL,font.size=12,font.face="plain",minors.threshold=0.05,
                   observations.sort=NULL,features.sort=NULL){
  if(logTrans|scale.by=="row"|scale.by=="col"){
    minors.threshold=NULL
  }
  #get continuous data into "cdata"
  cdata=as.data.frame(data[,cs])
  rownames(cdata)<-as.character(rownames(cdata))
  #get discrete data into "ddata"
  ddata=as.data.frame(data[ds])
  #initialize plot list with "p"
  p=list()
  #Heatmap branch
  if(type=="heatmap"){
    #if heatmap rows/cols is going to be clustered by corresponding data, then it would not be allowed
    #to sort rows/cols manually or arbitrarily, notice that features.sort only works for features of 
    #continuous data in heatmap
    if((cluster.by.y&!is.null(observations.sort))|(cluster.by.x&!is.null(features.sort))){
      stop("'cluster.by.x/y' is not compatible with sort.observations/features")
    }
    #set color bar width if the parameter was set, provided that the heatmap width was 15, then the default width of
    #a single colorbar would be 1, else a single colorbar width would be 'colorbar.width' that had been set.
    if(is.null(colorbar.width)){
      ratios=c(rep(1,length(ds)),15)
    }else{
      ratios=c(rep(colorbar.width,length(ds)),15)
    }
    if(logTrans){
      cdata=log2(cdata+1)
    }
    #scale data by row/col
    if(scale.by=="row"){
      cdata=t(scale(t(cdata)))
    }else if(scale.by=="col"){
      cdata=scale(cdata)
    }
    #heatmap part of the plot
    p1=ncplot(cdata,heatmap.colors=heatmap.main.colors,cluster.by.row=cluster.by.y,
              minors.threshold=minors.threshold,cluster.by.col=cluster.by.x,observations.sort=observations.sort,
              features.sort=features.sort)+
      theme(axis.text.x=element_text(angle=45,hjust=1,size=font.size,face=font.face),
            axis.title.y=element_blank(),
            legend.text=element_text(size=font.size+6,face=font.face),
            legend.title=element_text(size=font.size+6,face="bold"),
            axis.title.x=element_text(size=font.size+6,face="bold"))
    #determine if object names were going to be demonstrated
    if(!show.obname){
      p1=p1+theme(axis.text.y=element_blank())
    }
    #determine if feature names were going to be demonstrated
    if(!show.featname){
      p1=p1+theme(axis.text.x=element_blank())
    }
    oblevels=levels(p1$data[,1])
    p2legends=list()
    for(i in 1:ncol(ddata)){
      n.group=ddata[i] %>% table() %>% length()
      if(n.group>20){
        message("debug1")
        message(oblevels[1])
        print(head(ddata[i]))
        p[[i]]=ndplot(ddata[i],ob.level=oblevels,flip.coord = T)+
          scale_fill_gradientn(colours=head(legend.colors,2))+
          guides(fill=guide_colorbar(title=colnames(ddata[i])))
        legend.colors=tail(legend.colors,length(legend.colors)-2)
      }else{
        p[[i]]=ndplot(ddata[i],ob.level=oblevels,flip.coord = T)+
          scale_fill_manual(values=head(legend.colors,n.group))+
          guides(fill=guide_legend(title=colnames(ddata[i]),ncol=ceiling(n.group/4)))
        legend.colors=tail(legend.colors,length(legend.colors)-n.group)
      }
      p[[i]]=p[[i]]+theme(axis.title.y=element_blank(),
              axis.title.x=element_blank(),
              axis.text.y=element_blank(),
              axis.text.x=element_text(angle=45,hjust=1,size=font.size+6,face=font.face),
              legend.text=element_text(size=font.size+6,face=font.face),
              legend.title=element_text(size=font.size+6,face="bold"))
      p2legends[[i]]=get_legend(p[[i]])
      p[[i]]=p[[i]]+theme(legend.position="none")
    }
    p[[ncol(ddata)+1]]=p1
    p1legend=get_legend(p[[ncol(ddata)+1]])
    p2legends[[ncol(ddata)+1]]=p1legend
    legends=ggarrange(plotlist=p2legends,nrow=ncol(ddata)+1,align="v")+theme(plot.margin=margin(0,1,2,0.5,"cm"))
    p[[ncol(ddata)+1]]=p[[ncol(ddata)+1]]+theme(legend.position="none")
    main=ggarrange(plotlist=p,ncol=length(p),widths=ratios,align="h")
    ggarrange(plotlist=list(main,legends),ncol=2,widths=c(5,1))+ylab("observations")+
      theme(plot.margin=margin(0.5,0.5,0.5,0.5,"cm"))
  }else if(type=="bar"){
    if(is.null(colorbar.width)){
      ratios=c(rep(2,length(ds)),15)
    }else{
      ratios=c(rep(colorbar.width,length(ds)),15)
    }
    p1=ncplot(cdata,type = "bar",minors.threshold=minors.threshold,
              observations.sort=observations.sort,features.sort=features.sort)
    if(!show.obname){
      p1=p1+theme(axis.text.x=element_blank())
    }
    oblevels=levels(p1$data[,1])
    p2legends=list()
    for(i in 1:ncol(ddata)){
      message(paste0("plotting discrete data ",i))
      n.group=ddata[i] %>% table() %>% length()
      if(n.group>20){
        p[[i]]=ndplot(ddata[i],ob.level=oblevels,flip.coord = F)+
          scale_fill_gradientn(colours=head(legend.colors,2))+
          guides(fill=guide_colorbar(title=colnames(ddata[i]),barwidth=unit(200,"pt")))
        legend.colors=tail(legend.colors,length(legend.colors)-2)
      }else{
        p[[i]]=ndplot(ddata[i],ob.level=oblevels,flip.coord = F)+
          scale_fill_manual(values=head(legend.colors,n.group))+
          guides(fill=guide_legend(title=colnames(ddata[i]),nrow=2))
          legend.colors=tail(legend.colors,length(legend.colors)-n.group)
      }
      p[[i]]=p[[i]]+theme(axis.title.y=element_blank(),
            axis.title.x=element_blank(),
            axis.text.y=element_blank(),
            axis.text.x=element_blank(),
            legend.text=element_text(size=font.size+6,face=font.face),
            legend.title=element_text(size=font.size+6,face="bold"),
            legend.position="top")
      p2legends[[i]]=get_legend(p[[i]])
      p[[i]]=p[[i]]+theme(legend.position="none")
    }
    p[[ncol(ddata)+1]]=p1
#      p1legend=get_legend(p[[ncol(ddata)+1]])
#      p[[ncol(ddata)+1]]=p[[ncol(ddata)+1]]+theme(legend.position="none")
#      p2legends[[ncol(ddata)+1]]=p1legend
#      legends=ggarrange(plotlist=p2legends,ncol=ncol(ddata)+1,align="h")+theme(plot.margin=margin(0.5,0.5,0.5,0.5,"cm"))
      legends=ggarrange(plotlist=p2legends,ncol=ncol(ddata),align="h")+theme(plot.margin=margin(0.5,0.5,0.5,0.5,"cm"))
      main=ggarrange(plotlist=p,nrow=length(p),heights=ratios,align="v")
      ggarrange(plotlist=list(legends,main),nrow=2,heights=c(1,5))+xlab("observations")+
        theme(plot.margin=margin(0.5,0.5,0.5,0.5,"cm"))
  }else{
    stop("Unsupported type")
  }
}


#input: two data frame, one data frame with n continuous cols("cs" parameter) and m("ds" parameter) discrete cols(n!=0&m!=0),
#rows are observations, cols are features, the other data.frame is meta.data for features, with n rows, and the rownames
#should contain all feature names, i.e. the colnames of "cs" in "data".
ncmdldplot<-function(data,cs,ds,data2,type="heatmap",show.obname=T,show.featname=T,logTrans=F,scale.by="none",
                     cluster.by.y=FALSE,cluster.by.x=FALSE,heatmap.main.colors=c("black","white"),
                     legend.colors=pal_d3("category20")(20),colorbar.width=NULL,font.size=12,font.face="plain",
                     minors.threshold=0.05,observations.sort=NULL,features.sort=NULL){
  cdata=as.data.frame(data[,cs])
  ob.meta=as.data.frame(data[,ds])
  feat.meta=as.data.frame(data2)
  if(type=="heatmap"){
    #if heatmap rows/cols is going to be clustered by corresponding data, then it would not be allowed
    #to sort rows/cols manually or arbitrarily, notice that features.sort only works for features of 
    #continuous data in heatmap
    if((cluster.by.y&!is.null(observations.sort))|(cluster.by.x&!is.null(features.sort))){
      stop("'cluster.by.x/y' is not compatible with sort.observations/features")
    }
    #set color bar width if the parameter was set, provided that the heatmap width was 15, then the default width of
    #a single colorbar would be 1, else a single colorbar width would be 'colorbar.width' that had been set.
    if(is.null(colorbar.width)){
      ratios=c(rep(1,length(ds)),15)
    }else{
      ratios=c(rep(colorbar.width,length(ds)),15)
    }
    if(logTrans){
      cdata=log10(cdata+1)
    }
    #scale data by row/col
    if(scale.by=="row"){
      cdata=t(scale(t(cdata)))
    }else if(scale.by=="col"){
      cdata=scale(cdata)
    }
    #heatmap part of the plot
    p.main=ncplot(cdata,heatmap.colors=heatmap.main.colors,cluster.by.row=cluster.by.y,
              minors.threshold=minors.threshold,cluster.by.col=cluster.by.x,observations.sort=observations.sort,
              features.sort=features.sort)+
      theme(axis.text.x=element_text(angle=45,hjust=1,size=font.size,face=font.face),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),axis.title.x=element_blank(),
            legend.text=element_text(size=font.size,face=font.face),
            legend.title=element_text(size=font.size,face=font.face),
            legend.margin=margin(0,0,0,0,"cm"),
            legend.key.width=unit(25,"pt"))
    #determine if object names were going to be demonstrated
    if(!show.obname){
      p.main=p.main+theme(axis.text.y=element_blank())
    }
    #determine if feature names were going to be demonstrated
    if(!show.featname){
      p.main=p.main+theme(axis.text.x=element_blank())
    }
    oblevels=levels(p.main$data[,1])
    featlevels=levels(p.main$data[,2])
    
    #observation meta part of the whole plot
    p.ob.meta=list()
    p.ob.meta.legends=list()
    for(i in 1:ncol(ob.meta)){
      n.group=ob.meta[i] %>% table() %>% length()
      if(n.group>20){
        p.ob.meta[[i]]=ndplot(ob.meta[i],ob.level=oblevels,flip.coord = T)+
          scale_fill_gradientn(colours=head(legend.colors,2))+
          guides(fill=guide_colorbar(title=colnames(ob.meta[i]),barwidth=unit(200,"pt")))
        legend.colors=tail(legend.colors,length(legend.colors)-2)
      }else{
        p.ob.meta[[i]]=ndplot(ob.meta[i],ob.level=oblevels,flip.coord = T)+
          scale_fill_manual(values=head(legend.colors,n.group))+
          guides(fill=guide_legend(title=colnames(ob.meta[i]),nrow=2))
        legend.colors=tail(legend.colors,length(legend.colors)-n.group)
      }
      p.ob.meta[[i]]=p.ob.meta[[i]]+theme(axis.title.y=element_blank(),
              axis.title.x=element_blank(),
              axis.text.y=element_blank(),
              axis.text.x=element_text(angle=45,hjust=1,size=font.size,face=font.face),
              legend.text=element_text(size=font.size,face=font.face),
              legend.title=element_text(size=font.size,face=font.face),
              legend.margin=margin(0,0,0,0,"cm"))+
        guides(fill=guide_legend(title=colnames(ob.meta[i])))
      p.ob.meta.legends[[i]]=get_legend(p.ob.meta[[i]])
      p.ob.meta[[i]]=p.ob.meta[[i]]+theme(legend.position="none")
    }
    
    #feature meta part of the whole plot
    p.feat.meta=list()
    p.feat.meta.legends=list()
    for(i in 1:ncol(feat.meta)){
      n.group=feat.meta[i] %>% table() %>% length()
      if(n.group>10){
        p.feat.meta[[i]]=ndplot(feat.meta[i],ob.level=featlevels,flip.coord = F)+
          scale_fill_gradientn(colours=head(legend.colors,2))+
          guides(fill=guide_colorbar(title=colnames(feat.meta[i]),barwidth=unit(200,"pt")))
        legend.colors=tail(legend.colors,length(legend.colors)-2)
      }else{
        p.feat.meta[[i]]=ndplot(feat.meta[i],ob.level=featlevels,flip.coord = F)+
          scale_fill_manual(values=head(legend.colors,n.group))+
          guides(fill=guide_legend(title=colnames(feat.meta[i]),nrow=2))
        legend.colors=tail(legend.colors,length(legend.colors)-n.group)
      }
      p.feat.meta[[i]]=p.feat.meta[[i]]+theme(axis.title.y=element_blank(),
             axis.title.x=element_blank(),
             axis.text.y=element_blank(),
             axis.text.x=element_blank(),
             legend.text=element_text(size=font.size,face=font.face),
             legend.title=element_text(size=font.size,face=font.face))+
        guides(fill=guide_legend(title=colnames(feat.meta[i]),nrow=1,direction = "horizontal"))
      p.feat.meta.legends[[i]]=get_legend(p.feat.meta[[i]])
      p.feat.meta[[i]]=p.feat.meta[[i]]+theme(legend.position="none")
    }
    p.ob.meta=ggarrange(plotlist=p.ob.meta,ncol=ncol(ob.meta),align="h")
    p.main.legend=get_legend(p.main)
    p.main=p.main+theme(legend.position="none")
    p.ob.meta.legends=ggarrange(plotlist=p.ob.meta.legends,nrow=ncol(ob.meta),align="v")+
      theme(plot.margin=margin(1,0,0,0,"cm"))
    p.main.ob.legends=ggarrange(plotlist=list(p.main.legend,p.ob.meta.legends),
      heights=c(1,length(p.ob.meta.legends)),ncol=2)+
      theme(plot.margin=margin(0,0,3,0,"cm"))
    p.feat.meta=ggarrange(plotlist=p.feat.meta,nrow=ncol(feat.meta),align="v")
    p.feat.meta.legends=ggarrange(plotlist=p.feat.meta.legends,ncol=ncol(feat.meta),align="h")
    
    #arrange the whole plot and add blank plot into a Sudoku
    allplot=list()
    blankplot=ggplot()+geom_blank()+theme_void()
    allplot[[1]]=blankplot
    allplot[[2]]=p.feat.meta.legends
    allplot[[3]]=blankplot
    allplot[[4]]=blankplot
    allplot[[5]]=p.feat.meta
    allplot[[6]]=blankplot
    allplot[[7]]=p.ob.meta
    allplot[[8]]=p.main
    allplot[[9]]=p.main.ob.legends
    
    ggarrange(plotlist=allplot,nrow=3,ncol=3,
      widths=c(ncol(ob.meta)/4,ncol(cdata)*2,1),
      heights=c(1,ncol(feat.meta),nrow(cdata)))+theme(plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))
  }else{
    stop("Unsupported plot type")
  }
}

decorate_ccplot<-function(plot,base.size=13,x.limits=NULL,y.limits=NULL,x.expand=c(T,0,0),y.expand=c(T,0,0),
                          legend.pos=NULL,x.n.breaks = NULL,y.n.breaks = NULL){
  if(x.expand[1]&&y.expand[1]){
    plot=plot+scale_x_continuous(expand=x.expand[2:3],limits=x.limits,n.breaks=x.n.breaks)+
      scale_y_continuous(expand=y.expand[2:3],limits=y.limits,n.breaks=y.n.breaks)
  }else if((!x.expand[1])&&y.expand[1]){
    plot=plot+scale_y_continuous(expand=y.expand[2:3],limits=x.limits,n.breaks=x.n.breaks)
  }else if(x.expand[1]&&!y.expand[1]){
    plot=plot+scale_x_continuous(expand=x.expand[2:3],limits=y.limits,n.breaks=y.n.breaks)
  }else{
    plot=plot
  }
  plot=plot=plot+theme_publication(base.size=base.size,legend.pos=legend.pos)
  plot
}

theme_publication<-function(base.size=12,legend.pos="right",angle=45,hjust=1,x=T,margin.t=0.5,
                            margin.r=0.5,margin.b=0.5,margin.l=0.5){
  if(!x){
    theme=theme(text=element_text(size=base.size),
                axis.text.x=element_blank(),
                axis.line.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.line.y.left = element_line(),
                panel.background=element_rect(fill="white"),
                plot.margin = margin(margin.t,margin.r,margin.b,margin.l,"cm"),legend.position=legend.pos)
    return(theme)
  }
  theme(text=element_text(size=base.size),
    axis.text.x=element_text(angle=angle,hjust=hjust),
    axis.line.x = element_line(size = 0.5),axis.line.y = element_line(size = 0.5),
    panel.background=element_rect(fill="white"),
    plot.margin = margin(margin.t,margin.r,margin.b,margin.l,"cm"),legend.position=legend.pos)
}

dimplot_publication<-function(srt,group.by="seurat_clusters",colors=pal_d3(palette="category20")(20),reduction.name="umap",
                              reduction="umap",arrow_extend_ratio=0.25,label=F,raster=NULL){
  embedding_text=switch(reduction,"umap"="UMAP","tsne"="tSNE","pca"="PC","int"="integrated")
  embedings=srt@reductions[[reduction.name]]@cell.embeddings %>% as.data.frame()
  yintercetpt=min(embedings[,1])-2.5
  xintercept=min(embedings[,2])-2.5
  arrow_extends=ifelse(abs(yintercetpt)>abs(xintercept),abs(yintercetpt)*arrow_extend_ratio,abs(xintercept)*arrow_extend_ratio)
  DimPlot(srt,group.by=group.by,reduction=reduction.name,label=label,raster=raster)+scale_color_manual(values=colors)+
    geom_segment(aes(x=yintercetpt,xend=yintercetpt,y=xintercept,yend=xintercept+arrow_extends),
         arrow=arrow(length=unit(0.4,"cm"),type="closed",angle=15))+
    geom_segment(aes(x=yintercetpt,xend=yintercetpt+arrow_extends,y=xintercept,yend=xintercept),
         arrow=arrow(length=unit(0.4,"cm"),type="closed",angle=15))+
    xlab(paste0(embedding_text,1))+ylab(paste0(embedding_text,2))+ggtitle("")+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title=element_text(hjust=0.1,vjust=0.1,size=10))
}

theme_dimplot<-function(){
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())
}

generateColor<-function(n,rgb.max=1,rgb.min=0.1,seed=2015,col.dist.min=0.5,alpha=1,maxColorValue=1){
  set.seed(seed)
  details=data.frame(r=c(0,1),g=c(0,0),b=c(1,0))
  i=2
  colors=c(rgb(0,0,1,alpha=alpha,maxColorValue = maxColorValue),rgb(1,0,0,alpha=alpha,maxColorValue=maxColorValue))
  j=1
  while(i<n){
    rgb=runif(3,min=rgb.min,max=rgb.max)
    if(nrow(details)==0){
      detail=data.frame(r=rgb[1],g=rgb[2],b=rgb[3])
      colors=c(colors,rgb(rgb[1],rgb[2],rgb[3],alpha=alpha,maxColorValue=maxColorValue))
      details=rbind(details,detail)
      i=i+1
    }else{
      dist=details
      dist$r=dist$r-rgb[1]
      dist$g=dist$g-rgb[2]
      dist$b=dist$b-rgb[3]
      dist$dist=abs(dist$r)+abs(dist$g)+abs(dist$b)
      if(min(dist$dist)>col.dist.min){
        detail=data.frame(r=rgb[1],g=rgb[2],b=rgb[3])
        colors=c(colors,rgb(rgb[1],rgb[2],rgb[3],alpha=alpha,maxColorValue=maxColorValue))
        details=rbind(details,detail)
        i=i+1
      }else{
        j=j+1
        if(j>100000){
          stop("col.dist.min too large to obtain enough colors, set it smaller")
        }else{
          next
        }
      }
    }
  }
  #details$order=1:nrow(details)
  #details=arrange(details,r,g,b)
  return(colors)
}

