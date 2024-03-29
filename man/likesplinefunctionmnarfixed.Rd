% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likesplinefunctionmnarfixed.R
\name{likesplinefunctionmnarfixed}
\alias{likesplinefunctionmnarfixed}
\title{Output the MNAR MPLE likelihood}
\usage{
likesplinefunctionmnarfixed(
  dat,
  params,
  mnarparams,
  timevar,
  xvarsmiss,
  xvarscox,
  splinevars,
  ncause,
  ties = FALSE
)
}
\arguments{
\item{dat}{A data table that has observed time (failure or censoring), values of covariates in the missingness model, values of covariates in the Cox model, basis variables for the baseline hazard ratio spline model, a censoring indicator names 'cens' that takes value 0 if censored, and a cause variable named 'cause' that has a numeric value indicating the cause and takes 0 if cause is unknown}

\item{params}{The starting values for the missingness model and Cox model coefficients associated with the covariates in the dat data frame. Missingness model currently assumed to be logistic.}

\item{mnarparams}{The fixed values for the effect of missingness type on probability of missingness in a logistic model. This controls whether data are assumed to be MNAR}

\item{timevar}{The name of the time variable}

\item{xvarsmiss}{The names of the variables to be included in the missingness model}

\item{xvarscox}{The names of the variables to be included in the Cox model}

\item{splinevars}{The names of the variables that are the basis function for the spline model}

\item{ncause}{The number of different falure types (not including unknown failure type)}

\item{ties}{TRUE indicates there are tied failure times}
}
\description{
#'This takes dat for competing risks survival times, where failure type is missing for some observations
}
