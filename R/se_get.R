######################################################################
#'  @ Function for retrieving a table of standard errors from rrt
#'  @ models
#'  @ The function takes the following arguments
#'  @ intable=table from model selection results 
#'####################################################################
se_get<-function(intable){
    ml<-as.character(intable$modID)
    mod.l<-get(load("models"))
    resml<-NULL
    for (m in 1:length(ml)){
        mymod<-mod.l[ml[m]] # get model
        options(warn = -1)
        sest<-summary(mymod[[1]])$coefficients[,2] # summary with coefficients ALL coefficients here
        if(length(sest)==0){
            next
        } 
        resdf<-data.frame(modID=ml[m],varname=names(sest),se=sest)
        resml<-rbind(resdf,resml)
    }
    results<-spread(data=resml,key=varname,value=se)
    return(results)
}

