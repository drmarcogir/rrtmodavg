
######################################################################
#'  @ Function calculating model averaged predictions for rrt models
#'  @ Created by Marco Girardello 19/05/2016
#'  @ The function takes the following arguments
#'  @ intable=table from model selection results 
#'  @ index=model selection index to be used i.e. BIC or AIC 
#'####################################################################


predrrt<-function(intable=NULL,index=NULL){
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
  # get models
  mod.l<-get(load("models"))
  # get predictions
  pred.vec<-vector("list",dim(intable_95)[1])
  for (i in 1:dim(intable_95)[1]){
  mymod<-mod.l[intable_95[i,]$modID] # extract model
  options(warn = -1)
  tmpdf<-try(data.frame(ID=1:length(predict(mymod[[1]])[,1]),pred=predict(mymod[[1]])[,1],modID=intable_95[i,]$modID),silent=TRUE) 
  if(class(tmpdf)=="try-error"){
    next
  } else {
  pred.vec[[i]]<-tmpdf
  }
}
options(warn = -0)
  # bind results together
  preddf<-do.call("rbind",pred.vec)
  # arrange model order so that it is same when computing averages
  preddf1<-spread(preddf,modID, pred)
  preddf2<-preddf1[ , order(names(preddf1))]
  intable_95new<-arrange(intable_95,modID)
  intable_95new<-intable_95new[intable_95new$modID %in% names(preddf2)[-1],]
  # calculated model-averaged predictions 
  p.res<-NULL
  for (x in 1:dim(preddf2)[1]){
    p<-weighted.mean(preddf2[x,][-1],intable_95new$weightBIC)      
    p.res[x]<-p
  }
  
return(p.res)
  
}


