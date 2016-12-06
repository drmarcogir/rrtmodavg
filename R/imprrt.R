
######################################################################
#'  @ Function calculating variable importance for rrt models
#'  @ Created by Marco Girardello 19/05/2016
#'  @ The function takes the following arguments
#'  @ intable=table from model selection results 
#'  @ index=model selection index to be used i.e. BIC or AIC 
#'####################################################################


imprrt<-function(intable,index){
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
    resml<-vector("list",length(ml))
    # loop through each model
    for (y in 1:length(ml)){
      mymod<-mod.l[ml[y]] # extract model
      options(warn = -1)
      # get t-scores
      coefst<-summary(mymod[[1]])$coefficients[,3]
      # only for the parameter of interest
      cint<-coefst[c(names(intable_95a)[i])]
      # relative importance
      rat<-abs(cint)/max(abs(summary(mymod[[1]])$coefficients[,3][-1]))  
      # get model weight too
      df<-data.frame(rat,w=subset(tmp,modID==ml[y])$weightBIC)
      # store results
      resml[[y]]<-df
    }
    options(warn = -0)
    resml1<-do.call("rbind",resml)
    resml1<-na.exclude(resml1)
    imp<-sum(resml1$rat*resml1$w)/sum(resml1$w)
    impdf<-data.frame(var=names(intable_95a)[i],imp)
    res.avg[[i]]<-impdf
  }
  
  importance.final<-do.call("rbind",res.avg)
  return(importance.final)
}

    



