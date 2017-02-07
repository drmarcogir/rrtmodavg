
######################################################################
#'  @ Function calculating model averaged predictions for rrt models
#'  @ Created by Marco Girardello 19/05/2016
#'  @ The function takes the following arguments
#'  @ intable=table from model selection results 
#'  @ index=model selection index to be used i.e. BIC or AIC 
#'####################################################################


predrrt_rep<-function(intable=NULL,index=NULL,modfilt=NULL,runid=NULL){
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
        # get models
        
        # get predictions
        pred.vec<-vector("list",dim(modsubset)[1])
        modsubset$modID1<-paste(modsubset$modID,".fitted",sep="")
        
        for (i in 1:dim(modsubset)[1]){
            pred<-mod.l[modsubset[i,]$modID1] # extract model
            options(warn = -1)
            tmpdf<-try(data.frame(iteration=modsubset[i,]$iteration,pred=pred[[1]],modID=as.character(modsubset[i,]$modID))) 
            if(class(tmpdf)=="try-error"){
                next
            } else {
                pred.vec[[i]]<-tmpdf
            }
        }
        options(warn = -0)
        # bind results together
        preddf<-do.call("rbind",pred.vec)
        modid.l<-unique(preddf$modID)
        preddf1<-NULL
        for(b in 1:length(modid.l)){
            tmp1<-subset(preddf,modID==modid.l[b])
            tmp1$siteid<-1:dim(tmp1)[1]
            preddf1<-rbind(tmp1,preddf1)
        }
        
        preddf2<-merge(preddf1,tmp[,c("modID","weight")])
        site.l<-unique(preddf2$siteid)
        finalpred<-NULL
        for(m in 1:length(site.l)){
            tmp2<-subset(preddf2,siteid==site.l[m])
            p<-weighted.mean(tmp2$pred,tmp2$weight)         
            finalpred<-rbind(data.frame(p,siteid=site.l[m],iteration=x),finalpred)
        }
        
        results<-rbind(finalpred,results)
        
    }
    results1<-aggregate(p~siteid,data=results,FUN=mean)
    return(results1)
}