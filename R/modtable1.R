#'  Function for computing all possible model combinations using the RRlog function
#' 
#'  
#'  @ y=character vector containing the name of response 
#'  @ x= character vector containing the name of the predictors  
#'  @ df=name of dataframe 
#'  @ combos= number of variabiable combinations 2-variable models, 3-variable etc. if all variables 
#'  @ combinations are desired simply insert "all"
#'  @ inp = numerical vector containing randomization probability/probabilities
#'

modtable<- function(y=NULL, x=NULL, df=NULL, combos=NULL,inp=NULL,rep=NULL) {
  # if combos all do all models
  if(combos=="all"){
    combos=length(x)
  }
   if(is.null(inp)){
        inp=c(0.1, 0.1)
    } 
  results.final<-NULL # object where to store model selection table
  mods.final<-list() # object where to store all the models
  totlength <- length(x) # total length of predictors
  f <- sapply(df, is.factor)  # find out which variables are factors

  for (h in 1:rep){
      
print(paste("random sample set ",h,sep=""))  
ones<-subset(df,RRT2_all==1)
zeroes<-subset(df,RRT2_all==0)
row.names(zeroes)<-1:dim(zeroes)[1]
rr<-sample(1:dim(zeroes)[1],dim(ones)[1])
zeroes1<-zeroes[rr,]
comb<-rbind(ones,zeroes1)      

  for (i in 1:combos) {
    print(paste(i,"-variable models",sep=""))
    tmp <- combinations(totlength, i, x) # create model combinations
     # split combinations into list
    tmpl<-split(tmp, seq(nrow(tmp)))
    names(tmpl)<-paste("mod_",i,"_",1:dim(tmp)[1],"_",h,sep="") # labels for models
    ######################################################################
    fm<-function(predin){   # function for model fitting
      restmp <- matrix(nrow = 1, ncol = totlength + 3)  # matrix to store resuls
     colnames(restmp) <- c("(Intercept)", x, "AIC","BIC")
      predin1<-paste(predin,collapse="+")       # create formula
      formtmp <- as.formula(paste(y, "~", predin1, collapse = ""))
      mod <- RRlog(formtmp, data = df, model = "FR", p =inp, LR.test  =TRUE,fit.n = 1) # fit model
      if (class(mod) == "try-error") { # what to do if model did not converge
        failpar <- rep("fail", length(predin))
        names(failpar) <- predin
        cols <- match(names(failpar),colnames(restmp)) # match with matrix
        restmp[1, cols] <- failpar
      } else {
        n = length(predin) + 1
        aic = -2 * mod$logLik + 2 * n  # AIC
        bic=-2 * mod$logLik + log(dim(df)[1])*n
        resv <- c(mod$coefficients)  # get coefficients
        allpredf <- names(f[f == TRUE])  # coefficients for factors
        coef.fac <- allpredf[allpredf %in% predin]
        coefres <- rep("+", length(coef.fac))
        names(coefres) <- coef.fac
        coef.cont <- predin[!predin %in% coef.fac]  # other coeffcients
        resv1 <- c(resv[1], resv[coef.cont], coefres, aic,bic)
        names(resv1)[length(names(resv1))] <- "BIC"  # rename last element
        names(resv1)[length(names(resv1))-1] <- "AIC"   # rename element before last
        cols <- match(names(resv1),colnames(restmp))  # match with matrix
        restmp[1, cols] <- resv1
      }  # close else statement
      return(list(mod,restmp))
    }
    ############################
#   sfInit(parallel=TRUE, cpus=6)
#   sfLibrary('RRreg',character.only=TRUE)
#   sfExport('totlength',"y","x","tmpl","df","inp","f")
#     sfExport('x')
#     sfExport('y')
#     sfExport('tmpl')
#     sfExport('df')
#     sfExport('inp')
#    res<-sfLapply(tmpl, fm)
#    sfStop( nostop=FALSE )
    
   res<-lapply(tmpl,fm) # fit models usint list of formulae
    mods<-sapply(res, `[`,1)
    pars<-sapply(res, `[`,2)
    pars1<-do.call("rbind",pars) # create data frame from list
    pars2<-cbind(modID=names(pars),pars1)
    results.final<-rbind(pars2,results.final) # bind results to data frame
    mods.final<-c(mods,mods.final)
  } # end of i loop
} # end of h loop
  # save model objects
  mods.final1<-unlist(mods.final, recursive = FALSE)
  save(mods.final1,file="models")
  # calculate AIC and BIC weights
  write.csv(results.final,file="tmp.csv",row.names=F)
  results.final<-read.csv("tmp.csv")
  removed<-file.remove("tmp.csv")
  results.final1<-subset(results.final,!is.na(BIC)) 
  results.final1$deltaBIC<-results.final1$BIC - min(results.final1$BIC)   # delta BIC
  results.final1$weightBIC<-exp(-results.final1$deltaBIC / 2) / sum(exp(-results.final1$deltaBIC / 2))   # BIC weights
  results.final1<-subset(results.final1,!is.na(AIC))   # AIC weights
  results.final1$deltaAIC<- results.final1$AIC - min(results.final1$AIC)   # delta AIC
  results.final1$weightAIC<-exp(- results.final1$deltaAIC / 2) / sum(exp(- results.final1$deltaAIC / 2))   # AIC weights
  # change second column (always Intercept!)
  colnames(results.final1)[2]<-"(Intercept)"
  # change column names for factorial variables (same as the ones spat out by the model)
  # change names of columns for factors (match those of the model)
#  fv<-paste(names(f[f==TRUE]))
#  for(z in 1:length(fv)){
#   formtmp <- as.formula(paste(y, "~", fv[z], collapse = ""))
#   mod <- try(RRlog(formtmp, data = df, model = "FR", p = c(0.1, 0.1), LR.test  =TRUE,fit.n = 1))
#   modname<-names(summary(mod)$coefficients[,1][-1])
#   cols <- match(fv[x],colnames(results.final1))  # match with matrix
#   colnames(results.final1)[cols] <- fv[x]
#  }
  return(results.final1)
}
    


