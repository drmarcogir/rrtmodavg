#' Calculate variable importance for rrt models
#' 
#' @ intable=table from model selection results 
#' @ index=model selection index to be used i.e. BIC or AIC 
#'

imprrt<-function(intable=NULL,index=NULL,method=NULL){
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
    res.avg<-NULL  # change into append!
    for (i in 1:length(res.avg)){
        n<-c(names(intable_95a)[i],index,"modID") 
        tmp<-data.frame(intable_95a[i],intable_95[,n[2:3]])
        colnames(tmp)[1]<-names(intable_95a)[i]
        tmp<-na.exclude(tmp)
        mod.l<-get(load("models"))
        ml<-as.character(tmp$modID)
        resml<-NULL
        for (y in 1:length(ml)){
            mymod<-mod.l[ml[y]] # get model
            options(warn = -1)
            coefst<-summary(mymod[[1]])$coefficients[,3] # summary with coefficients ALL coefficients here
            if(length(coefst)==0){
                next
            }  else {              
                coefnm<-names(coefst) # extract names
                dd<-stri_detect_fixed(coefnm, names(intable_95a)[i]) # partial string matching
                dd1<-coefst[dd] # extract specific coefficients (including factors)
                rat<-abs(dd1)/max(abs(summary(mymod[[1]])$coefficients[,3][-1])) # calculation of relative importance
                df<-data.frame(varname=names(rat),rat,w=subset(tmp,modID==ml[y])[,index]) # put together with model weight
                resml<-rbind(df,resml) 
            }
        }
        } # end of i loop
        options(warn = -0)
        resml1<-na.exclude(resml)
        var.l<-unique(resml1[,c("varname")])
        for (p in 1:length(var.l)){
             resml2<-resml1[resml1$varname %in% var.l[p],]
        imp<-sum(resml2$rat*resml2$w)/sum(resml2$w)
        impdf<-data.frame(var=var.l[p],imp)
        res.avg<-rbind(impdf,res.avg)
    }
        ############################## burnand
    return(res.avg)
    }

