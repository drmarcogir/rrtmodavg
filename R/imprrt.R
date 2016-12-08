#' Calculate variable importance for rrt models
#' 
#' @ The function takes the following arguments
#' @ intable=table from model selection results 
#' @ index=model selection index to be used i.e. BIC or AIC 
#'

imprrt<-function(intable=NULL,index=NULL){
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
    res.avg<-vector("list",length(totcols1))  
    for (i in 1:length(res.avg)){
        n<-c(names(intable_95a)[i],index,"modID") 
        tmp<-data.frame(intable_95a[i],intable_95[,n[2:3]])
        colnames(tmp)[1]<-names(intable_95a)[i]
        tmp<-na.exclude(tmp)
        mod.l<-get(load("models"))
        ml<-as.character(tmp$modID)
        resml<-vector("list",length(ml))
        for (y in 1:length(ml)){
            mymod<-mod.l[ml[y]] 
            options(warn = -1)
            coefst<-summary(mymod[[1]])$coefficients[,3]
            cint<-coefst[c(names(intable_95a)[i])]
            rat<-abs(cint)/max(abs(summary(mymod[[1]])$coefficients[,3][-1]))  
            df<-data.frame(rat,w=subset(tmp,modID==ml[y])$weightBIC)
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
