######################################################################
#'  Function calculating importance of variables for rrt models 
#'  when multiple replications are present
#'  The function takes the following arguments
#'  @ intable=table from model selection results 
#'  @ index=model selection index to be used i.e. BIC or AIC 
#'####################################################################

imp_rep<-function(intable=NULL,index=NULL,modfilt=NULL,runid=NULL){
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
        if(modfilt==3){
            modsubset<-intable 
        } 
        tmp2<-melt(modsubset[,c(1,3:15)],id.vars=c("modID"))
        tmp3<-subset(tmp2,!is.na(value))
        tmp4<-ddply(tmp3,.(variable),nrow)
        colnames(tmp4)[2]<-c("count")
        tmp4$count<-tmp4$count/dim(tmp1)[1]
        results<-rbind(tmp4,results)
        
    }  
    results1<-aggregate(count~variable,FUN=mean,data=results)
return(results1)
}
    