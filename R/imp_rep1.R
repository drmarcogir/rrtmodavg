#' Calculate variable importance for rrt models
#' 
#' @ intable=table from model selection results 
#' @ index=model selection index to be used i.e. BIC or AIC 
#'

imp_rep1<-function(intable=NULL,index=NULL,modfilt=NULL,runid=NULL){
    mod.l<-get(load("models"))
    results<-NULL                  # store results of predictions
    iter.l<-unique(intable[,runid])  # list ID for iterations
    for (x in 1:length(iter.l)){      # loop through randomized model iterations 
        tmp<-intable[intable[,runid] %in% iter.l[x],]
        tmp$delta<-tmp[,index] - min(tmp[,index])
        tmp$weight<-exp(-tmp$delta/ 2) / sum(exp(-tmp$delta / 2)) 
        if(modfilt==1){
            tmp$cumsum<-cumsum(tmp$weight)
            rown<-subset(tmp,cumsum >= 0.95)
            modsubset<-intable[intable$modID %in% rown$modID,]
        } 
        if(modfilt==2){
            tmp1<-subset(tmp,delta < 4)
            modsubset<-intable[intable$modID %in% tmp1$modID,]   
        } 
        
        tmp2<-melt(tmp1[,c(1,3:15,26)],id.vars=c("modID","weight"))
        tmp3<-na.exclude(tmp2)
        impfinal<-aggregate(weight~variable,data=tmp3,FUN=sum)
        impfinal1<-data.frame(impfinal,iteration=x)
        results<-rbind(impfinal1,results)
    }
    results1<-aggregate(weight~variable,data=results,FUN=mean)
    colnames(results1)[2]<-c("importance")
    return(results1)
}
