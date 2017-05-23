######################################################################
#'  @ Function calculating model averaged parameters for rrt models
#'  @ Created by Marco Girardello 19/05/2016
#'  @ The function takes the following arguments
#'  @ intable=table from model selection results 
#'  @ index=model selection index to be used i.e. BIC or AIC 
#'####################################################################
modavgrrt<-function(intable=NULL,index=NULL,bin.factors=TRUE,delta4=NULL){
    if((index=="AIC")==TRUE){
        index=c("weightAIC")
    }
    if((index=="BIC")==TRUE){
        index=("weightBIC")
    }
    if(delta4==FALSE){
        intable$cumsum<-cumsum(intable[,index])
        rown<-subset(intable,cumsum >= 0.95)
        rown1<-min(as.numeric(row.names(rown)))
        intable_95<-intable[1:rown1,]
        totcols<-1:dim(intable)[2]
        totcols1<-totcols[2:(length(totcols)-7)]
        intable_95a<-intable_95[,totcols1]    
    }
    if(delta4==TRUE){
        index1<-paste("delta",substring(index, 7, 10),sep="")
        intable_95<-intable[intable[,index1] <= 4,]
        totcols<-1:dim(intable)[2]
        totcols1<-totcols[2:(length(totcols)-7)]
        intable_95a<-intable_95[,totcols1]  
    }
    # get only data for specific variables
#   nums <- sapply(intable_95a, is.numeric)
#   intable_95a<-intable_95a[ , nums]
    res.avg<-vector("list",dim(intable_95a)[2])  
    for (i in 1:length(res.avg)){
        if(i==1){
        colnames(intable_95a)[1]<-c("(Intercept)")    
        }
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
                totlength<-length(names(intable_95a)[i])
                fact <- sapply(intable_95a, is.factor) # extract factors
                if((fact[names(intable_95a)[i]])==TRUE){
                dd<-stri_detect_fixed(tmpdf$name,names(intable_95a)[i]) # partial string matching
                pardf<-tmpdf[dd,]
                pardf<-data.frame(pardf,modID=ml[y])
                res.ml[[y]]<-pardf
                } else {
                pardf<-tmpdf[names(intable_95a)[i],]
                pardf<-data.frame(pardf,modID=ml[y])
                res.ml[[y]]<-pardf
                }
            }
        }
        options(warn = -0)
        res.ml1<-do.call("rbind",res.ml) 
        findf<-merge(res.ml1,tmp[2:3]) 
        findf1<-par.avg(x=findf$par,se=findf$se,weight=findf[,index])
        findf2<-data.frame(variable=unique(findf$name),as.data.frame(findf1),type=names(findf1))
        findf3<-spread(findf2,type,findf1)
        res.avg[[i]]<-findf3 
    }
    
    table<-do.call("rbind",res.avg)
    table1<-table[,-c(2)]
    table2<-table1[,c(1,2,4,3,5)]
    return(table2)
}


