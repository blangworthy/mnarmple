#'Output the MNAR MPLE likelihood
#'
#'#'This takes dat for competing risks survival times, where failure type is missing for some observations
#'@param dat A data table that has observed time (failure or censoring), values of covariates in the missingness model, values of covariates in the Cox model, basis variables for the baseline hazard ratio spline model, a censoring indicator names 'cens' that takes value 0 if censored, and a cause variable named 'cause' that has a numeric value indicating the cause and takes 0 if cause is unknown
#'@param params The starting values for the missingness model and Cox model coefficients associated with the covariates in the dat data frame. Missingness model currently assumed to be logistic.
#'@param mnarparams The fixed values for the effect of missingness type on probability of missingness in a logistic model. This controls whether data are assumed to be MNAR
#'@param timevar The name of the time variable
#'@param xvarsmiss The names of the variables to be included in the missingness model
#'@param xvarscox The names of the variables to be included in the Cox model
#'@param splinevars The names of the variables that are the basis function for the spline model
#'@param ncause The number of different falure types (not including unknown failure type)

likesplinefunctionmnarfixed <- function(dat,params,mnarparams,timevar,xvarsmiss,xvarscox,splinevars,ncause){
  ###Order based on timevar if so we can use cumsum for denominator
  dat <- dat[order(-dat[[timevar]])]
  ties <- sum(duplicated(dat[[timevar]]))>0
  dat <- dummy_cols(dat,select_columns = "cause")
  ###parameters for misingness model
  ax <- params[grep("ax",names(params))]
  ay <- mnarparams
  ###parameters for cox model
  for(i in 1:ncause){
    assign(paste0("b",i),params[grep(paste0("b",i),names(params))])
  }
  ###parameters for baseline hazard ratio
  for(i in 2:ncause){
    assign(paste0("bas",i),params[grep(paste0("bas",i),names(params))])
  }
  b <- do.call(cbind,mget(paste0("b",1:ncause)))
  bas <- do.call(cbind,mget(paste0("bas",2:ncause)))
  ###Cox portion of model
  coxbeta <- exp(tcrossprod(as.matrix(dat[,..xvarscox]),t(b)))
  ###Cumulative sum of cox model for portion for denominator
  coxbetasum <- apply(coxbeta,2,cumsum)
  colnames(coxbeta) <- paste0("beta",1:ncause)
  colnames(coxbetasum) <- paste0("betasum",1:ncause)
  ###Baseline hazard ratio
  basrat <- cbind(tcrossprod(as.matrix(dat[,..splinevars]),t(bas)))
  colnames(basrat) <- paste0("eta",2:ncause)
  ###Missingness model
  obsx <- tcrossprod(as.matrix(dat[,..xvarsmiss]),t(ax))
  pobs <- inv.logit(sapply(ay,function(k) obsx + k))
  colnames(pobs) <- paste0("pobs",1:ncause)
  dat <- data.table(dat,coxbeta,coxbetasum,basrat,pobs)
  ###Denominator of likelihood
  dat[,denom :=  eval(parse(text=paste0("betasum1+", paste(paste0("betasum",2:ncause),paste0("eta",2:ncause),sep="*",collapse = "+"))))]

  if(ties == TRUE){
    dat[,denom:= max(denom),by=timevar]
  }
  ###Numerator of likelihood
  dat[,num := eval(parse(text=paste0("pobs1*beta1*cause_1 +",paste(paste0("pobs",2:ncause),paste0("eta",2:ncause),paste0("beta",2:ncause),paste0("cause_",2:ncause),sep="*",collapse = "+"),
                                     "+(1-pobs1)*beta1*cause_0+",paste(paste0("(1-pobs",2:ncause,")"),paste0("eta",2:ncause),paste0("beta",2:ncause),paste0("cause_0"),sep="*",collapse = "+"))))]
  ###Log likelihood
  dat[,llik := log((num/denom)*cens + 1*(1-cens))]

  return(-1*sum(dat$llik))

}
