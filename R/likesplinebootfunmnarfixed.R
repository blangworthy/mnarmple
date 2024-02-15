#'Output the MNAR MPLE likelihood to be used with boot function
#'
#'#'This takes dat for competing risks survival times, where failure type is missing for some observations and is compatible with boot() function
#'@param dat A data table that has observed time (failure or censoring), values of covariates in the missingness model, values of covariates in the Cox model, basis variables for the baseline hazard ratio spline model, a censoring indicator names 'cens' that takes value 0 if censored, and a cause variable named 'cause' that has a numeric value indicating the cause and takes 0 if cause is unknown
#'@param params The starting values for the missingness model and Cox model coefficients associated with the covariates in the dat data frame. Missingness model currently assumed to be logistic.
#'@param mnarparams The fixed values for the effect of missingness type on probability of missingness in a logistic model. This controls whether data are assumed to be MNAR
#'@param optmethod Method to be used in optim() function
#'@param timevar The name of the time variable
#'@param xvarsmiss The names of the variables to be included in the missingness model
#'@param xvarscox The names of the variables to be included in the Cox model
#'@param splinevars The names of the variables that are the basis function for the spline model
#'@param ncause The number of different falure types (not including unknown failure type)
#'@param indices Indices to be used for boot() function
#'@param nmaxit Maximum iterations for optim()
likesplinebootfunmnarfixed <- function(data,params,mnarparams,optmethod,timevar,xvarsmiss,xvarscox,splinevars,ncause,indices,nmaxit=500){

  d <- data[indices,]
  ests <- optim(params,method = optmethod,fn=likesplinefunctionmnarfixed,dat=d,mnarparams=mnarparams,timevar=timevar,xvarsmiss=xvarsmiss,xvarscox=xvarscox,splinevars=splinevars,ncause=ncause,control=list(maxit=nmaxit))
  return(ests$par)

}
