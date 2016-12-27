######################################################################
#'  @ Function calculating model averaged parameters for rrt models
#'  @ Created by Marco Girardello 19/05/2016
#'  @ The function takes the following arguments
#'  @ intable=table from model selection results 
#'  @ index=model selection index to be used i.e. BIC or AIC 
#'####################################################################
modavgrrt<-function(intable=NULL,index=NULL){
    if((index=="AIC")==TRUE){
        index=c("weightAIC")
    }
    if((index=="BIC")==TRUE){
        index=("weightBIC")
    }
    intable$cumsum<-cumsum(intable[,index])
    rown<-subset(intable,cumsum >= 0.95)
    rown1<-min(as.numeric(row.names(rown)))
    intable_95<-intable[1:rown1,]
    totcols<-1:dim(intable)[2]
    totcols1<-totcols[2:(length(totcols)-7)]
    intable_95a<-intable_95[,totcols1]
    # get only data for specific variables
    nums <- sapply(intable_95a, is.numeric)
    intable_95a<-intable_95a[ , nums]
    res.avg<-vector("list",dim(intable_95a)[2])  
    for (i in 1:length(res.avg)){
        n<-c(names(intable_95a)[i],index,"modID") 
        tmp<-data.frame(intable_95a[i],intable_95[,n[2:3]])
        colnames(tmp)[1]<-names(intable_95a)[i]
        tmp<-na.exclude(tmp)
        mod.l<-get(load("models"))
        ml<-as.character(tmp$modID)
        res.ml<-vector("list",length(ml))
        for (y in 1:length(ml)){
            mymod<-mod.l[ml[y]] 
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
        res.ml1<-do.call("rbind",res.ml) 
        findf<-merge(res.ml1[2:4],tmp[2:3]) 
        findf1<-par.avg(x=findf$par,se=findf$se,weight=findf[,index]) 
        findf1<-data.frame(variable=names(intable_95a)[i],t(as.data.frame(findf1))) 
        res.avg[[i]]<-findf1 
    }
    
    table<-do.call("rbind",res.avg)
    row.names(table)<-1:dim(table)[1]
    table1<-table[,-c(4)]
    return(table1)
}


