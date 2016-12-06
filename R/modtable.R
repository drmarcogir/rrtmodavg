##################################################################################
#'  @ Function for computing all possible model combinations using the RRlog function
#'  @ Created by Marco Girardello 19/05/2016
#'  @ The function takes the following arguments
#'  @ y=character vector containing the name of response 
#'  @ x= character vector containing the name of the predictors  
#'  @ df=name of dataframe 
#'  @ combos= number of variabiable combinations 2-variable models, 3-variable etc.
#'################################################################################


modtable<- function(y, x, df, combos) {
  # if combos all do all models
  if(combos=="all"){
    combos=length(x)
  }
  results.final<-NULL # object where to store model selection table
  mods.final<-vector("list",combos) # object where to store all the models
  totlength <- length(x) # total length of predictors
  f <- sapply(df, is.factor)  # find out which variables are factors 
  for (i in 1:combos) {
    print(paste(i,"-variable models",sep=""))
    tmp <- combinations(totlength, i, x) # create model combinations
     # split combinations into list
    tmpl<-split(tmp, seq(nrow(tmp)))
    names(tmpl)<-paste("mod_",i,"_",1:dim(tmp)[1],sep="") # labels for models
    fm<-function(predin){   # function for model fitting
      restmp <- matrix(nrow = 1, ncol = totlength + 3)  # matrix to store resuls
      colnames(restmp) <- c("(Intercept)", x, "AIC","BIC")
      predin1<-paste(predin,collapse="+")       # create formula
      formtmp <- as.formula(paste(y, "~", predin1, collapse = ""))
      mod <- try(RRlog(formtmp, data = df, model = "FR", p = c(0.1, 0.1), LR.test  =TRUE,fit.n = 1)) # fit model
      if (class(mod) == "try-error") { # what to do if model did not converge
        failpar <- rep("fail", length(namesest))
        names(failpar) <- namesest
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
    res<-lapply(tmpl,fm) # fit models usint list of formulae
    mods<-sapply(res, `[`,1)
    pars<-sapply(res, `[`,2)
    pars1<-do.call("rbind",pars) # create data frame from list
    pars2<-cbind(modID=names(pars),pars1)
    results.final<-rbind(pars2,results.final) # bind results to data frame
    mods.final[[i]]<-mods
  }
  # save model objects
  mods.final1<-unlist(mods.final, recursive = FALSE)
  save(mods.final1,file="models")
  # calculate AIC and BIC weights
  write.csv(results.final,file="tmp.csv",row.names=F)
  results.final<-read.csv("tmp.csv")
  system("rm tmp.csv")
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
  fv<-paste(names(f[f==TRUE]))
  for(x in 1:length(fv)){
    formtmp <- as.formula(paste(y, "~", fv[x], collapse = ""))
    mod <- try(RRlog(formtmp, data = df, model = "FR", p = c(0.1, 0.1), LR.test  =TRUE,fit.n = 1))
    modname<-names(summary(mod)$coefficients[,1][-1])
    cols <- match(fv[x],colnames(results.final1))  # match with matrix
    colnames(results.final1)[cols] <- modname
  }
  return(results.final1)
}
    


