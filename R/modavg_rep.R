######################################################################
#'  @ Function calculating model averaged coefficients for rrt models
#'  when multiple runs with different subsets of data were used
#'  @ Created by Marco Girardello 19/05/2016
#'  @ The function takes the following arguments
#'  @ intable=table from model selection results 
#'  @ index=model selection index to be used i.e. BIC or AIC 
#'####################################################################

modavg_rep<-function(intable=NULL,index=NULL,modfilt=NULL,runid=NULL,cols=NULL){
    coefs.l<-get(load("coefsl"))
    vcov.l<-get(load("vcovl"))
    results<-NULL                  # store results of predictions
    iter.l<-unique(intable[,runid])  # list ID for iterations    
    for (x in 1:length(iter.l)){      # loop through randomized model iterations 
        print(x)
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
            modsubset<-tmp
        } 
        
        res.avg<-vector("list",14) 
        var.l<-names(modsubset)[cols]
        var.l[1]<-c("(Intercept)")   
        for (i in 1:length(res.avg)){ # for each variable
            print(var.l[i])
            if(i==1){
                colnames(modsubset)[2]<-c("(Intercept)")    
            }
                n<-c(var.l[i],index,"modID") # columns to select
                tmp<-data.frame(modsubset[,var.l[i]],modsubset[,n[2:3]]) # select columns for variable of interest
                colnames(tmp)[1]<-names(modsubset[var.l[i]]) # rename intercept column
                tmp<-na.exclude(tmp) # exclude NAs i.e. models not fitted
                ml<-as.character(tmp$modID) # model list to loop through
                res.ml<-vector("list",length(ml)) # store parameter values for model averaging
                for (y in 1:length(ml)){
                    mlcoef<-paste(ml[y],".coefficients",sep="") # create names to subset coefficient values
                    coefs<-coefs.l[mlcoef]  # subset coefficients
                    mlse<-paste(ml[y],".vcov",sep="") # create names to subset standard error values
                    sevcov<-vcov.l[mlse] # subset standard errors
                    se <- sqrt(abs(diag(sevcov[[1]]))) # calculate standard errors
                    if(class(se)=="matrix"){
                        next
                    } else {
                        options(warn = -1)
                        tmpdf<-data.frame(par=coefs[[1]],se=se) # data frame containing parameter estimates & se for model of interest!
                        if(dim(tmpdf)[1]==0){
                            next
                        } else {
                            tmpdf$name<-row.names(tmpdf) # add column with name
                            namedf<-data.frame(name=modsubset[var.l[i]])  
                            totlength<-length(names(modsubset[var.l[i]]))
                            fact <- sapply(modsubset[,cols], is.factor) # extract factors
                            fact[names(modsubset[var.l[i]])]
                            if(fact[names(modsubset[var.l[i]])]==FALSE){
                                #dd<-stri_detect_fixed(tmpdf$name,names(intable_95a)[i])
                                dd<-(tmpdf$name==names(modsubset[var.l[i]])) # partial string matching
                                pardf<-tmpdf[dd,]
                                pardf<-data.frame(pardf,modID=ml[y])
                                res.ml[[y]]<-pardf
                            } else {
                                pardf<-tmpdf[names(modsubset[var.l[i]]),]
                                pardf<-data.frame(pardf,modID=ml[y])
                                res.ml[[y]]<-pardf
                            }
                        }
                    }
                    
                }
                options(warn = -0)
                res.ml1<-do.call("rbind",res.ml) # contains parameter estimates for variable of interest
                if(is.null(res.ml1)){
                next    
                } else {
                findf<-merge(res.ml1,tmp[2:3]) # merge with BIC or AIC values
                findf1<-par.avg(x=findf$par,se=findf$se,weight=findf[,index]) # performs model averaging
                findf2<-data.frame(variable=unique(findf$name),as.data.frame(findf1),type=names(findf1))
                findf3<-spread(findf2,type,findf1)
                res.avg[[i]]<-findf3
                }
            }
            
            table<-do.call("rbind",res.avg)
            table1<-table[,-c(2)]
            table2<-table1[,c(1,2,4,3,5)]
            table2$iteration<-iter.l[x]
            results<-rbind(table2,results)
    }
    
    return(table2)           
}

