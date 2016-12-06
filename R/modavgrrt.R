######################################################################
#'  @ Function calculating model averaged parameters for rrt models
#'  @ Created by Marco Girardello 19/05/2016
#'  @ The function takes the following arguments
#'  @ intable=table from model selection results 
#'  @ index=model selection index to be used i.e. BIC or AIC 
#'####################################################################


modavgrrt<-function(intable,index){
  if((index=="AIC")==TRUE){
    index=c("weightAIC")
  }
  if((index=="BIC")==TRUE){
    index=("weightBIC")
  }
  # 95% cumsum 
  intable$cumsum<-cumsum(intable[,index])
  rown<-subset(intable,cumsum >= 0.95)
  rown1<-min(as.numeric(row.names(rown)))
  intable_95<-intable[1:rown1,]
  # exclude first and last six columns
  totcols<-1:dim(intable)[2]
  totcols1<-totcols[2:(length(totcols)-7)]
  intable_95a<-intable_95[,totcols1]
  # vector to store results
  res.avg<-vector("list",length(totcols1))  
  # loop through each variable
  for (i in 1:length(res.avg)){
    # get columns with specific variable, weight and modelID
    n<-c(names(intable_95a)[i],index,"modID") 
    tmp<-data.frame(intable_95a[i],intable_95[,n[2:3]])
    colnames(tmp)[1]<-names(intable_95a)[i]
    tmp<-na.exclude(tmp)
    # get model object
    mod.l<-get(load("models"))
    # results for model objects
    ml<-as.character(tmp$modID)
    res.ml<-vector("list",length(ml))
    # loop through each model
    for (y in 1:length(ml)){
      mymod<-mod.l[ml[y]] # extract model
      options(warn = -1)
      tmpdf<-data.frame(par=summary(mymod[[1]])$coefficients[,1],se=summary(mymod[[1]])$coefficients[,2])
      if(dim(tmpdf)[1]==0){
       next
      } else {
      tmpdf$name<-row.names(tmpdf)
      namedf<-data.frame(name=names(intable_95a)[i])
      pardf<-merge(namedf,tmpdf)
      pardf<-data.frame(pardf,modID=ml[y])
      res.ml[[y]]<-pardf
      }
    }
    options(warn = -0)
    res.ml1<-do.call("rbind",res.ml) # bind results for parameter estimate
    findf<-merge(res.ml1[2:4],tmp[2:3]) # merge parameters data frame with model weight data frame
    findf1<-par.avg(x=findf$par,se=findf$se,weight=findf[,index]) # model averaging (uses MuMin par.avg)
    findf1<-data.frame(variable=names(intable_95a)[i],t(as.data.frame(findf1))) # insert variable name
    res.avg[[i]]<-findf1 # store results
  }
  
  table<-do.call("rbind",res.avg)
  row.names(table)<-1:dim(table)[1]
  table1<-table[,-c(4)]
  return(table1)
}
  

