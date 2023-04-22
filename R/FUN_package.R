#' Fit negative-binomial regression to estimate unrecognized deaths attributable to
#' certain infectious disease
#'
#' This function fits negative binomial regression models using measures of virus transmission
#' as covariates to estimate unrecognized deaths attributable to infections caused by those virus
#' in epidemiological studies with a time-series design.
#'
#' @param formula  an object of class "formula" specifying the model to be fitted
#' using a symbolic description.
#' @param data a data frame containing variables used for fitting models.
#' @param offset an optional argument specifying an \emph{a prior} known component to be
#' included in the linear predictor. For example, log of at-risk populations.
#' @param niter a number specifying the number of MCMC iterations, the default is 5000.
#' @param burn_in a number specifying the number of iterations that is discarded, the default is 2500.
#' @param rand.int a logical value indicating whether time-specific random effects are included.
#' @param str.int a list containing two character strings when \code{rand.int} = \code{TRUE}.
#' The first character string specifies the variable indicating the clusters of
#' time-specific random effects, this variable must be included in the \code{data};
#' the second character string specifies the structure of time-specific random effects,
#' \code{AR1} and \code{AR2} fit random walk model of order 1 or 2,
#' and \code{exch} specifies exchangeable time-specific random effects.
#' @param covariate.AD a vector containing character strings specifying which variables need
#' to be set to 0 to estimate attributable deaths.
#' @param countsonly a logical value indicating whether posterior samples of random effects are
#' returned when \code{rand.int} = \code{TRUE}, the default is \code{FALSE}.
#' @param ID.spacetime an optional vector containing variable names specifying location and time indicators of 
#' observations, those indicators have to be included in the \code{data}, the default is \code{NULL}. 
#' Specifying location and time indicators is important if the user want to obtain estimates of 
#' attributable deaths across a certain time period and locations. 
#'
#' @return This function returns a list containing:\tabular{ll}{
#'    \code{beta} \tab a matrix containing posterior samples of fixed effects and mean of time-varying 
#'    coefficients when \code{slope.tvarying} = \code{TRUE}. \cr
#'    \tab \cr
#'    \code{parm} \tab a matrix containing posterior samples of the over-dispersion parameter \code{xi}, 
#'    and parameters controlling temporal dependency of time-specific random effects when \code{rand.int} = \code{TRUE}. 
#'    When \code{exch} is specified for time-specific random effects, \code{tau2.t} is included; 
#'    when \code{AR1} or \code{AR2} are specified, \code{tau2.t} and \code{rho.t} are included. \cr
#'    \tab \cr
#'    \code{pred.counts} \tab a matrix containing posterior samples with burn-in samples are discarded 
#'    of observed counts. \cr
#'    \tab \cr
#'    \code{pred.null.counts} \tab a matrix containing posterior samples with burn-in samples are discarded 
#'    of estimated counts assuming variables specified using \code{covariate.AD} equal to 0. \cr
#'    \tab \cr
#'    \code{AD} \tab a matrix containing posterior samples with burn-in samples are discarded 
#'    of attributable deaths 
#'    (i.e., point-wise difference between \code{pred.counts} and \code{pred.null.counts}). \cr
#'    \tab \cr
#'    \code{accept.rate} \tab acceptance rate for updating the over-dispersion parameter 
#'    in Metropolis-Hasting algorithm. \cr
#'    \tab \cr
#'    \code{WAIC} \tab the widely applicable information criterion (WAIC) 
#'    in deviance scale that measures how the model fit to the observed data. \cr
#'    \tab \cr
#'    \code{lppd} \tab the log point-wise predictive probability of the 
#'    observed data under the fitted model. \cr
#'    \tab \cr
#'    \code{pWAIC} \tab the approximated effective number of parameters; \cr
#' }
#' when \code{countsonly} = \code{FALSE}, this function also returns the following: \tabular{ll}{
#'    \code{rand.int.t} \tab a matrix containing posterior samples of time-specific random effects
#'     when \code{rand.int} = \code{TRUE}. \cr
#'    \tab \cr
#' }    
#' 
#'
#' @references
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive 
#' information criteria for Bayesian models. 
#' \emph{Statistics and computing, 24}(6), 997-1016.
#' 
#'
#' @examples
#' data(toydata)
#' 
#' # NO time-specific random effects 
#' fit.1 <- fit.NB.timeseries(formula = y ~ lag1 + lag2, 
#' data = toydata, offset = logPop, 
#' niter = 5000, burn_in = 2500,
#' rand.int = FALSE,
#' covariate.AD = c("lag1", "lag2"),
#' countsonly = TRUE) 
#' 
#' # Assume time-specific random effects modeled 
#' # via random walk model of order 1
#' fit.2 <- fit.NB.timeseries(formula = y ~ lag1 + lag2, 
#' data = toydata, offset = logPop, 
#' niter = 5000, burn_in = 2500,
#' rand.int = TRUE, str.int = list("week", "AR1"),
#' covariate.AD = c("lag1", "lag2"),
#' countsonly = TRUE) 
#' 
#' @export
#' @import Rcpp
#' @import RcppArmadillo
#' @import lme4
#' @import Matrix
#' @import MASS
#' @importFrom igraph get.adjacency graph.edgelist
#' @import dplyr
#' @import truncnorm
#' @importFrom BayesLogit rpg
#' @import splines
### A function to fit negative-binomial regression to
### estimate attributable deaths (AD) with time-series designs
fit.NB.timeseries <- function(formula, data, offset = NULL,
                              niter = 5000, burn_in = 2500,
                              rand.int = TRUE, str.int = NULL,
                              covariate.AD,
                              countsonly = FALSE,
                              ID.spacetime = NULL, ... ){

  ### INPUTS:
  ## formula: y ~ x1 + x2
  ## data: a dataframe containing data
  ## offset: specify a variable that will be used as an offset term in the model
  ## niter: the number of MCMC iterations
  ## burn_in: the number of burn-in samples
  ## rand.int: a logic value indicating whether time-specific intercepts are specified
  ## str.int: a list list containing: (1) a variable name indicating
  ###         the clusters of time-specific intercepts; (2) structure: exchangeable/AR-1/AR-2
  ## covariate.AD: a vector of strings which specifies variable names corresponding
  ### covariates that need to set to zero when computing attributable deaths
  ## countsonly: a logic value indicating whether posterior samples of
  ###   random effects (i.e., time-specific intercepts or time-varying coefs) are returned



  pars.all <- as.list(match.call()[-1]) # remove the function name

  outcome.name <- all.vars(formula)[1]
  pred.names <- all.vars(formula)[-1]

  offset.name <- as.character(pars.all$offset)
  if(length(offset.name) != 0){
    offset.vec <- as.matrix(data[,offset.name])
    offset.vec <- as.numeric(offset.vec)
  }else{
    offset.vec <- NULL
  }
  outcome.vec <- as.matrix(data[,outcome.name])
  outcome.vec <- as.numeric(outcome.vec)
  num.obs <- length(outcome.vec)

  if(!all(covariate.AD %in% pred.names)){
    stop("variables must be included in the formula")
  }

  if(!is.null(ID.spacetime)){
    if(!all(ID.spacetime %in% colnames(data))){
      stop("space and time indicators of observations must be inlucded in the dataset")
    }
  }

  if(is.null(ID.spacetime)){
    ID.spacetime.df <- NULL
  }else{
    ID.spacetime.df <- data[,ID.spacetime]
  }


  ###### SCENARIO 1: no random intercept and no random slopes ######
  if(rand.int == FALSE){
    X.design.fix <- model.matrix(formula, data = data)
    re0 <- fit.NB.time.s4(outcome.vec = outcome.vec,
                          X.design.fix = X.design.fix,
                          offset.vec = offset.vec,
                          niter = niter, burn_in = burn_in,
                          covariate.AD = covariate.AD,
                          countsonly = countsonly,
                          ID.spacetime.df = ID.spacetime.df)
  }



  ###### SCENARIO 2: random intercept and no random slopes ######
  if(rand.int == TRUE){
    warning.text = paste0("please provide a list containing",
                          " a variable name indicating",
                          " the clusters of random intercept",
                          " and the structure of the random intercept")
    warning.text1 = paste0("the structure can only be specified as ",
                           "exch, AR1, or AR2")
    warning.text2 = paste0("the variable name indicating",
                           " the clusters of random intercept",
                           " must be inlcuded in the dataset")
    if(is.null(str.int)){
      stop(warning.text)
    }
    if(!is.null(str.int)){
      if(!is.list(str.int)){
        stop(warning.text)
      }
      ll = length(str.int)
      if(ll != 2){
        stop(warning.text)
      }
      if(ll == 2){
        if(!(str.int[[2]] %in% c("exch", "AR1", "AR2"))){
          stop(warning.text1)
        }
        if(!(str.int[[1]] %in% colnames(data))){
          stop(warning.text2)
        }
      }
    }
    id.int.name = str.int[[1]]
    ntimes.int = length(unique(as.matrix(data[,id.int.name])))
    if(ntimes.int == num.obs){
      X.design.fix <- model.matrix(formula, data = data)
      X.design.rand.int <- Diagonal( x = rep(1, num.obs))
    }else if(ntimes.int < num.obs){
      comp.vec <- paste0("(1|", id.int.name, ")")
      formula.final <- as.formula(paste(c(deparse(formula), comp.vec),
                                        collapse = "+"))
      X.all <- lme4::lFormula(formula = formula.final, data = data)
      X.design.rand.int <- t(X.all$reTrms$Zt)
      X.design.fix <- X.all$X
    }

    re0 <- fit.NB.time.s2(outcome.vec = outcome.vec,
                          X.design.fix = X.design.fix,
                          X.design.rand.int = X.design.rand.int,
                          offset.vec = offset.vec,
                          rand.int.str = str.int[[2]],
                          ntimes = ntimes.int,
                          niter = niter, burn_in = burn_in,
                          covariate.AD = covariate.AD,
                          countsonly = countsonly,
                          ID.spacetime.df = ID.spacetime.df)
  } # END SCENARIO 2
  
  return(re0)
  
}





#' Fit space-time additive negative-binomial regression to 
#' estimate unrecognized deaths attributable to certain infectious disease
#' based on spatial-temporal data
#'
#' This function fits space-time additive negative binomial regression models 
#' using measures of virus transmission
#' as covariates to estimate unrecognized deaths attributable to infections caused 
#' by those virus
#' in epidemiological studies collecting data with spatial-temporal structures.
#'
#' @param formula  an object of class "formula" specifying the model to be fitted
#' using a symbolic description.
#' @param data a data frame containing variables used for fitting models.
#' @param offset an optional argument specifying an \emph{a prior} known component to be
#' included in the linear predictor. For example, log of at-risk populations.
#' @param niter a number specifying the number of MCMC iterations, the default is 5000.
#' @param burn_in a number specifying the number of iterations that is discarded, the default is 2500.
#' @param rand.int a logical value indicating whether time-specific random effects are included.
#' @param str.int a list containing two character strings when \code{rand.int} = \code{TRUE}.
#' The first character string specifies the variable indicating the clusters of
#' time-specific random effects, this variable must be included in the \code{data};
#' the second character string specifies the structure of time-specific random effects,
#' \code{AR1} and \code{AR2} fit random walk model of order 1 or 2,
#' and \code{exch} specifies exchangeable time-specific random effects.
#' @param slope.tvarying a logical value indicating whether time-varying coefficients are included.
#' @param str.tvarying a list containing vectors including three character strings
#' when \code{slope.tvarying} = \code{TRUE}.
#' The first character string specifies which variable having time-varying regression coefficient;
#' the second character string specifies the variable indicating the clusters of
#' time-varying coefficient, this variable must be included in the \code{data};
#' the third character string specifies the structure of time-varying coefficient,
#' \code{AR1} and \code{AR2} fit random walk model of order 1 or 2,
#' and \code{exch} specifies exchangeable random effects for time-varying coefficients.
#' @param covariate.AD a vector containing character strings specifying which variables need
#' to be set to 0 to estimate attributable deaths.
#' @param countsonly a logical value indicating whether posterior samples of random effects are
#' returned when \code{rand.int} = \code{TRUE} and/or \code{slope.tvarying} = \code{TRUE}, 
#' the default is \code{FALSE}.
#' @param spatial.str a list containing two character strings. The 
#' first character string specifies the variable indicating the clusters of
#' location-specific random effects, this variable must be included in the \code{data};
#' the second character string specifies the structure of location-specific random effects,
#' \code{exch} specifies exchangeable location-specific random effects, 
#' and \code{CAR} fits conditional autoregressive model (CAR).
#' @param adj.mat.s an optional matrix, when the CAR model is fitted, 
#' an between-location adjacency matrix must to be provided. 
#'
#' @return This function returns a list containing:\tabular{ll}{
#'    \code{beta} \tab a matrix containing posterior samples of fixed effects and mean of time-varying 
#'    coefficients when \code{slope.tvarying} = \code{TRUE}. \cr
#'    \tab \cr
#'    \code{parm} \tab a matrix containing posterior samples of the over-dispersion parameter \code{xi}, 
#'    parameters controlling temporal dependency of time-specific random effects 
#'    when \code{rand.int} = \code{TRUE}, and parameters controlling spatial dependency of 
#'    location-specific random effects. When \code{exch} is specified for time-specific random effects, 
#'    \code{tau2.t} is included; when \code{AR1} or \code{AR2} are specified for time-specific random effects, 
#'    \code{tau2.t} and \code{rho.t} are included. When \code{exch} is specified for 
#'    location-specific random effects, \code{tau2.s} is included; 
#'    when \code{CAR} is specified for location-specific random effects, 
#'    \code{tau2.s} and \code{rho.s} are incldued.\cr
#'    \tab \cr
#'    \code{pred.counts} \tab a matrix containing posterior samples with burn-in samples are discarded 
#'    of observed counts. \cr
#'    \tab \cr
#'    \code{pred.null.counts} \tab a matrix containing posterior samples with burn-in samples are discarded 
#'    of estimated counts assuming variables specified using \code{covariate.AD} equal to 0. \cr
#'    \tab \cr
#'    \code{AD} \tab a matrix containing posterior samples with burn-in samples are discarded 
#'    of attributable deaths 
#'    (i.e., point-wise difference between \code{pred.counts} and \code{pred.null.counts}). \cr
#'    \tab \cr
#'    \code{accept.rate} \tab acceptance rate for updating the over-dispersion parameter 
#'    in Metropolis-Hasting algorithm. \cr
#'    \tab \cr
#'    \code{WAIC} \tab the widely applicable information criterion (WAIC) 
#'    in deviance scale that measures how the model fit to the observed data. \cr
#'    \tab \cr
#'    \code{lppd} \tab the log point-wise predictive probability of the 
#'    observed data under the fitted model. \cr
#'    \tab \cr
#'    \code{pWAIC} \tab the approximated effective number of parameters; \cr
#' }
#' when \code{countsonly} = \code{FALSE}, this function also returns the following: \tabular{ll}{
#'    \code{rand.int.t} \tab a matrix containing posterior samples of time-specific random effects
#'     when \code{rand.int} = \code{TRUE}. \cr
#'    \tab \cr
#'    \code{rand.int.s} \tab a matrix containing posterior samples of 
#'    location-specific random effects. \cr
#'    \tab \cr
#'    \code{beta.tvarying} \tab a list containing matrices of posterior samples 
#'    of time-varying coefficients when \code{slope.tvarying} = \code{TRUE}. \cr
#'    \tab \cr
#'    \code{parm.beta.tvarying} \tab a list containing matrices of posterior samples of 
#'    parameters controlling temporal dependency of time-varying coefficients when 
#'    \code{slope.tvarying} = \code{TRUE}; when \code{exch} is specified, \code{sigma2.X} is included; 
#'    when \code{AR1} or \code{AR2} are specified, \code{sigma2.X} and \code{rho.X} are included, 
#'    where \code{X} is replaced by variable names specified via \code{covariate.AD}. \cr
#' }    
#' 
#'
#' @references
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive 
#' information criteria for Bayesian models. 
#' \emph{Statistics and computing, 24}(6), 997-1016.
#' 
#'
#' @examples
#' data(toydata)
#' 
#' # Assume time-specific random effects and 
#' # assume time-varying coefficients
#' fit.st <- fit.NB.st(formula = y ~ lag1 + lag2,
#' data = dat.fit.NB,
#' offset = logPop,
#' niter = 5000, burn_in = 2500,
#' rand.int = TRUE,
#' str.int = list("week", "exch"),
#' slope.tvarying = TRUE, 
#' str.tvarying = list(c("lag1", "week", "exch"),
#'                     c("lag2", "week", "exch")),
#' covariate.AD = c("lag1", "lag2"),
#' spatial.str = list("states", "exch"), adj.mat.s = NULL)
#' 
#' @export
### A function to fit space-time additive negative-binomial regression model
### to estimate attributable deaths (AD) based on spatial-temporal data
fit.NB.st <- function(formula, data, offset = NULL, niter = 5000, 
                      burn_in = 2500,
                      rand.int = TRUE, str.int = NULL,
                      slope.tvarying = FALSE, str.tvarying = NULL,
                      covariate.AD,
                      countsonly = FALSE,
                      spatial.str, adj.mat.s = NULL,
                      ID.spacetime = NULL, ...){

  ### INPUTS:
  ## formula: y ~ x1 + x2
  ## data: a dataframe containing data
  ## offset: specify a variable that will be used as an offset term in the model
  ## niter: the number of MCMC iterations
  ## burn_in: the number of burn-in samples
  ## str.int: a list containing: (1) a variable name indicating
  ###         the clusters of time-specific intercepts; (2) structure: exchangeable/AR-1/AR-2
  ## slope.tvarying: a logic value indicating whether time-varying coefs are specifieds
  ## str.tvarying: a list list containing: (1) a variable name indicating
  ###         the clusters of time-varying coefs;
  ### (2) structure: exchangeable/AR-1/AR-2 can only be "exch", "AR1", "AR2"
  ## covariate.AD: a vector of strings which specifies variable names corresponding
  ### covariates that need to set to zero when computing attributable deaths
  ## countsonly: a logic value indicating whether posterior samples of
  ###   random effects (i.e., time-specific intercepts or time-varying coefs) are returned
  ## spatial.str: a list containing: (1) a variable name indicating
  ### the clusters of location-specific intercepts;
  ### (2) structure: exchangeable/CAR can only be "exch" or "CAR"
  ## adj.mat.s: an adjacency matrix for location-specific intercepts, must be provided when
  ## spatial.str = c("", "CAR")



  pars.all <- as.list(match.call()[-1]) # remove the function name

  outcome.name <- all.vars(formula)[1]
  pred.names <- all.vars(formula)[-1]

  offset.name <- as.character(pars.all$offset)
  if(length(offset.name) != 0){
    offset.vec <- as.matrix(data[,offset.name])
    offset.vec <- as.numeric(offset.vec)
  }else{
    offset.vec <- NULL
  }
  outcome.vec <- as.matrix(data[,outcome.name])
  outcome.vec <- as.numeric(outcome.vec)
  num.obs <- length(outcome.vec)

  if(!all(covariate.AD %in% pred.names)){
    stop("variables must be included in the formula")
  }

  if(spatial.str[[2]] == "CAR" & is.null(adj.mat.s)){
    stop("an adjacency matrix for location-specific random effects must be specified")
  }

  if(!is.null(ID.spacetime)){
    if(!all(ID.spacetime %in% colnames(data))){
      stop("space and time indicators of observations must be inlucded in the dataset")
    }
  }

  if(is.null(ID.spacetime)){
    ID.spacetime.df <- NULL
  }else{
    ID.spacetime.df <- data[,ID.spacetime]
  }

  warning.text.a = paste0("please provide a list containing",
                          " a variable name indicating",
                          " the clusters of location-specific random effects",
                          " and the structure of location-specific random effects")
  warning.text1.a = paste0("the structure can only be specified as ",
                           "exch, or CAR")
  warning.text2.a = paste0("the variable name indicating",
                           " the clusters of location-specific random effects",
                           " must be inlcuded in the dataset")

  if(is.null(spatial.str)){
    stop(warning.text.a)
  }
  if(!is.null(spatial.str)){
    if(!is.list(spatial.str)){
      stop(warning.text.a)
    }
    ll = length(spatial.str)
    if(ll != 2){
      stop(warning.text.a)
    }
    if(ll == 2){
      if(!(spatial.str[[2]] %in% c("exch", "CAR"))){
        stop(warning.text1.a)
      }
      if(!(spatial.str[[1]] %in% colnames(data))){
        stop(warning.text2.a)
      }
    }
  }
  id.loc.name = spatial.str[[1]]
  nlocs = length(unique(as.matrix(data[,id.loc.name])))
  if(nlocs == num.obs){
    X.design.spatial <- Diagonal( x = rep(1, num.obs))
  }else if(nlocs < num.obs){
    comp.vec <- paste0("(1|", id.loc.name, ")")
    formula.final <- as.formula(paste(c(deparse(formula), comp.vec),
                                      collapse = "+"))
    X.all <- lme4::lFormula(formula = formula.final, data = data)
    X.design.spatial <- t(X.all$reTrms$Zt)
  }

  ###### SCENARIO 1: no temporal random effects and no time-varying coefs ######
  if(rand.int == FALSE & slope.tvarying == FALSE){
    X.design.fix <- model.matrix(formula, data = data)
    re0 <- fit.NB.st.add.s4(outcome.vec = outcome.vec,
                            X.design.fix = X.design.fix,
                            offset.vec = offset.vec,
                            niter = niter, burn_in = burn_in,
                            covariate.AD = covariate.AD,
                            countsonly = countsonly,
                            X.design.spatial = X.design.spatial,
                            spatial.str = spatial.str[[2]],
                            W.s = adj.mat.s,
                            ID.spacetime.df = ID.spacetime.df)
  } # END SCENARIO 1


  ###### SCENARIO 2: no temporal random effects and time-varying coefs ######
  if(rand.int == FALSE & slope.tvarying == TRUE){
    warning.text.a = paste0("please provide a list containing",
                            " vectors that specify: ",
                            " (1) variable names of time-varying coefs,",
                            " (2) variable names indicating",
                            " the clusters of time-varying coefs",
                            " , and (3) the structure of time-varying coefs")
    warning.text.b = paste0("the variable name indicating",
                            " the clusters of time-varying coefs",
                            " must be inlcuded in the dataset")
    warning.text.c = paste0("variables having",
                            " time-varying coefs",
                            " must be inlcuded in the model")

    if(is.null(str.tvarying)){
      stop(warning.text.a)
    }
    if(!is.null(str.tvarying)){
      if(!is.list(str.tvarying)){
        stop(warning.text.a)
      }
      ll.vec = sapply(str.tvarying, length)
      if(!all(ll.vec == 3)){
        stop(warning.text.a)
      }
      check <- sapply(str.tvarying, function(x) x[1])
      check.1 <- sapply(str.tvarying, function(x) x[2])
      check.2 <- sapply(str.tvarying, function(x) x[3])
      if(!all(check %in% pred.names)){
        stop(warning.text.c)
      }
      if(!all(check.1 %in% colnames(data))){
        stop(warning.text.b)
      }
      if(!all(check.2 %in% c("exch", "AR1", "AR2"))){
        stop(warning.text1)
      }
    }

    X.design.fix <- model.matrix(formula, data = data)
    var.slop.name = sapply(str.tvarying, function(x) x[1])
    id.slop.name = sapply(str.tvarying, function(x) x[2])
    str.slop.name = sapply(str.tvarying, function(x) x[3])

    X.design.rand.slop <- list()
    for(lll in 1:length(var.slop.name)){
      ntimes.slop.lll = length(unique(as.matrix(data[,id.slop.name[lll]])))
      if(ntimes.slop.lll == num.obs){
        X.design.rand.slop[[lll]] <- Diagonal( x = rep(1, num.obs))
      }else if(ntimes.slop.lll < num.obs){
        comp.vec <- paste0("(0+", var.slop.name[lll] ,"|", id.slop.name[lll], ")")
        formula.lll <- as.formula(paste(c(deparse(formula), comp.vec),
                                        collapse = "+"))
        X.inter <- lme4::lFormula(formula = formula.lll, data = data)
        X.design.rand.slop[[lll]] <- t(X.inter$reTrms$Zt)
      }
    }
    names(X.design.rand.slop) <- var.slop.name
    var.fix.name = colnames(X.design.fix)

    re0 <- fit.NB.st.add.s3(outcome.vec = outcome.vec,
                            X.design.fix = X.design.fix,
                            X.design.rand.slop = X.design.rand.slop,
                            offset.vec = offset.vec,
                            rand.slop.str = str.slop.name,
                            niter = niter, burn_in = burn_in,
                            covariate.AD = covariate.AD,
                            countsonly = countsonly,
                            X.design.spatial = X.design.spatial,
                            spatial.str = spatial.str[[2]],
                            W.s = adj.mat.s,
                            ID.spacetime.df = ID.spacetime.df)
  }


  ## SCENARIO 3: random intercept and no time-varying coefs
  if(rand.int == TRUE & slope.tvarying == FALSE){

    warning.text = paste0("please provide a list containing",
                          " a variable name indicating",
                          " the clusters of time-specific random effects",
                          " and the structure of time-specific random effects")
    warning.text1 = paste0("the structure can only be specified as ",
                           "exch, AR1, or AR2")
    warning.text2 = paste0("the variable name indicating",
                           " the clusters of time-specific random effects",
                           " must be inlcuded in the dataset")

    if(is.null(str.int)){
      stop(warning.text)
    }
    if(!is.null(str.int)){
      if(!is.list(str.int)){
        stop(warning.text)
      }
      ll = length(str.int)
      if(ll != 2){
        stop(warning.text)
      }
      if(ll == 2){
        if(!(str.int[[2]] %in% c("exch", "AR1", "AR2"))){
          stop(warning.text1)
        }
        if(!(str.int[[1]] %in% colnames(data))){
          stop(warning.text2)
        }
      }
    }
    id.int.name = str.int[[1]]
    ntimes.int = length(unique(as.matrix(data[,id.int.name])))
    if(ntimes.int == num.obs){
      X.design.fix <- model.matrix(formula, data = data)
      X.design.rand.int <- Diagonal( x = rep(1, num.obs))
    }else if(ntimes.int < num.obs){
      comp.vec <- paste0("(1|", id.int.name, ")")
      formula.final <- as.formula(paste(c(deparse(formula), comp.vec),
                                        collapse = "+"))
      X.all <- lme4::lFormula(formula = formula.final, data = data)
      X.design.rand.int <- t(X.all$reTrms$Zt)
      X.design.fix <- X.all$X
    }

    re0 <- fit.NB.st.add.s2(outcome.vec = outcome.vec,
                            X.design.fix = X.design.fix,
                            X.design.rand.int = X.design.rand.int,
                            offset.vec = offset.vec,
                            rand.int.str = str.int[[2]],
                            niter = niter, burn_in = burn_in,
                            covariate.AD = covariate.AD,
                            countsonly = countsonly,
                            X.design.spatial = X.design.spatial,
                            spatial.str = spatial.str[[2]],
                            W.s = adj.mat.s,
                            ID.spacetime.df = ID.spacetime.df)

  } # END SCENARIO 2


  ###### SCENARIO 4: random intercept and random slopes ######
  if(rand.int == TRUE & slope.tvarying == TRUE){

    warning.text = paste0("please provide a list containing",
                          " a variable name indicating",
                          " the clusters of time-specific random effects",
                          " and the structure of time-specific random effects")
    warning.text1 = paste0("the structure can only be specified as ",
                           "exch, AR1, or AR2")
    warning.text2 = paste0("the variable name indicating",
                           " the clusters of time-specific random effects",
                           " must be inlcuded in the dataset")
    if(is.null(str.int)){
      stop(warning.text)
    }
    if(!is.null(str.int)){
      if(!is.list(str.int)){
        stop(warning.text)
      }
      ll = length(str.int)
      if(ll != 2){
        stop(warning.text)
      }
      if(ll == 2){
        if(!(str.int[[2]] %in% c("exch", "AR1", "AR2"))){
          stop(warning.text1)
        }
        if(!(str.int[[1]] %in% colnames(data))){
          stop(warning.text2)
        }
      }
    }

    warning.text.a = paste0("please provide a list containing",
                            " vectors that specify: ",
                            " (1) variable names of time-varying coefs,",
                            " (2) variable names indicating",
                            " the clusters of time-varying coefs",
                            " , and (3) the structure of time-varying coefs")
    warning.text.b = paste0("the variable name indicating",
                            " the clusters of time-varying coefs",
                            " must be inlcuded in the dataset")
    warning.text.c = paste0("variables having",
                            " time-varying coefs",
                            " must be inlcuded in the model")

    if(is.null(str.tvarying)){
      stop(warning.text.a)
    }
    if(!is.null(str.tvarying)){
      if(!is.list(str.tvarying)){
        stop(warning.text.a)
      }
      ll.vec = sapply(str.tvarying, length)
      if(!all(ll.vec == 3)){
        stop(warning.text.a)
      }
      check <- sapply(str.tvarying, function(x) x[1])
      check.1 <- sapply(str.tvarying, function(x) x[2])
      check.2 <- sapply(str.tvarying, function(x) x[3])
      if(!all(check %in% pred.names)){
        stop(warning.text.c)
      }
      if(!all(check.1 %in% colnames(data))){
        stop(warning.text.b)
      }
      if(!all(check.2 %in% c("exch", "AR1", "AR2"))){
        stop(warning.text1)
      }
    }

    warning.text.a = paste0("please provide a list containing",
                            " a variable name indicating",
                            " the clusters of location-specific random effects",
                            " and the structure of location-specific random effects")
    warning.text1.a = paste0("the structure can only be specified as ",
                             "exch, or CAR")
    warning.text2.a = paste0("the variable name indicating",
                             " the clusters of location-specific random effects",
                             " must be inlcuded in the dataset")

    if(is.null(spatial.str)){
      stop(warning.text.a)
    }
    if(!is.null(spatial.str)){
      if(!is.list(spatial.str)){
        stop(warning.text.a)
      }
      ll = length(spatial.str)
      if(ll != 2){
        stop(warning.text.a)
      }
      if(ll == 2){
        if(!(spatial.str[[2]] %in% c("exch", "CAR"))){
          stop(warning.text1.a)
        }
        if(!(spatial.str[[1]] %in% colnames(data))){
          stop(warning.text2.a)
        }
      }
    }

    id.int.name = str.int[[1]]
    ntimes.int = length(unique(as.matrix(data[,id.int.name])))
    if(ntimes.int == num.obs){
      X.design.fix <- model.matrix(formula, data = data)
      X.design.rand.int <- Diagonal( x = rep(1, num.obs))
    }else if(ntimes.int < num.obs){
      comp.vec <- paste0("(1|", id.int.name, ")")
      formula.final <- as.formula(paste(c(deparse(formula), comp.vec),
                                        collapse = "+"))
      X.all <- lme4::lFormula(formula = formula.final, data = data)
      X.design.rand.int <- t(X.all$reTrms$Zt)
      X.design.fix <- X.all$X
    }

    var.slop.name = sapply(str.tvarying, function(x) x[1])
    id.slop.name = sapply(str.tvarying, function(x) x[2])
    str.slop.name = sapply(str.tvarying, function(x) x[3])

    X.design.rand.slop <- list()
    for(lll in 1:length(var.slop.name)){
      ntimes.slop.lll = length(unique(as.matrix(data[,id.slop.name[lll]])))
      if(ntimes.slop.lll == num.obs){
        X.design.rand.slop[[lll]] <- Diagonal( x = rep(1, num.obs))
      }else if(ntimes.slop.lll < num.obs){
        comp.vec <- paste0("(0+", var.slop.name[lll] ,"|", id.slop.name[lll], ")")
        formula.lll <- as.formula(paste(c(deparse(formula), comp.vec),
                                        collapse = "+"))
        X.inter <- lme4::lFormula(formula = formula.lll, data = data)
        X.design.rand.slop[[lll]] <- t(X.inter$reTrms$Zt)
      }
    }
    names(X.design.rand.slop) <- var.slop.name
    var.fix.name = colnames(X.design.fix)

    id.loc.name = spatial.str[[1]]
    nlocs = length(unique(as.matrix(data[,id.loc.name])))
    if(nlocs == num.obs){
      X.design.spatial <- Diagonal( x = rep(1, num.obs))
    }else if(nlocs < num.obs){
      comp.vec <- paste0("(1|", id.loc.name, ")")
      formula.final <- as.formula(paste(c(deparse(formula), comp.vec),
                                        collapse = "+"))
      X.all <- lme4::lFormula(formula = formula.final, data = data)
      X.design.spatial <- t(X.all$reTrms$Zt)
    }

    re0 <- fit.NB.st.add.s1(outcome.vec = outcome.vec,
                            X.design.fix = X.design.fix,
                            X.design.rand.int = X.design.rand.int,
                            X.design.rand.slop = X.design.rand.slop,
                            offset.vec = offset.vec,
                            rand.int.str = str.int[[2]],
                            rand.slop.str = str.slop.name,
                            niter = niter, burn_in = burn_in,
                            covariate.AD = covariate.AD,
                            countsonly = countsonly,
                            X.design.spatial = X.design.spatial,
                            spatial.str = spatial.str[[2]],
                            W.s = adj.mat.s,
                            ID.spacetime.df = ID.spacetime.df)
  } # END SCENARIO 4 (space-time additive)

  return(re0)
  ### OUTPUTS: a list contains:
  ## (1) "beta": posterior samples of fixed effects or
  ###    mean of time-varying coefficients
  ## (2) "parm": posterior samples of parameters controlling temporal dependency
  ###  of time-specific random effects,
  ###   parameters controlling spatial dependency of location-specific random effects,
  ###   and the over-dispersion parameter
  ## (3) "rand.int.t": posterior samples of time-specific random effects
  ###     (when countsonly = FALSE)
  ## (4) "rand.int.s":  posterior samples of location-specific effects
  ###     (when countsonly = FALSE)
  ## (5) "beta.tvarying": a list contains posterior samples of
  ### time-varying coefs (when slope.tvarying = TRUE and countsonly = FALSE)
  ## (6) "parm.beta.tvarying": a list contains posterior samples of
  ### parameters controlling temporal dependency
  ###  of time-varying coefs (when slope.tvarying = TRUE and countsonly = FALSE)
  ## (7) "pred.counts": posterior samples (burn-in samples are discarded) of
  ### predicted outcome counts
  ## (8) "pred.null.counts": posterior samples (burn-in samples are discarded) of
  ### predicted deaths assuming no certain infections
  ## (9) "AD": posterior samples (burn-in samples are discarded) of
  ### attributable deaths
  ## (10): "accept.rate": acceptance rate of M-H algorithm applied for updating
  ### the over-dispersion parameter
  ## (11): "WAIC": the WAIC measuring how the model fits to the data
}




#' Summarize results from estimating attributable deaths
#'
#' This function summarizes results of negative-binomial models fitted to 
#' the observed data to obtain posterior mean, posterior standard deviation, and 
#' 95\% credible interval of estimated counts  across time-period of interest 
#' and at locations of interest.
#'
#' @param counts.post  a matrix of posterior samples of counts, this matrix also
#' has to include location and time index of each counts in a wide format.
#' @param ID.aggre a character string specifying the grouping variable.
#' @param timeperiod a vector specifying the time period of interest.
#' @param locs a vector specifying locations of interest.
#' @param ID.time a character string specifying the time index of posterior samples.
#' @param ID.loc a character string specifying the location index of posterior samples.
#'
#' @return This function returns a data frame containing:\tabular{ll}{
#'    \code{est} \tab posterior mean of estimated counts. \cr
#'    \tab \cr
#'    \code{sd} \tab posterior standard deviation of estimated counts.\cr
#'    \tab \cr
#'    \code{lci} \tab lower bound of 95\% credible interval for estimated counts. \cr
#'    \tab \cr
#'    \code{uci} \tab upper bound of 95\% credible interval for estimated counts. \cr
#' }
#'
#' @examples
#' data(toydata)
#' 
#' ### Note that this example is just for the illustration, 
#' ### the model doesn't fit the data well
#' fit <- fit.NB.timeseries(formula = y ~ lag1 + lag2, 
#' data = toydata, 
#' offset = logPop, niter = 5000, burn_in = 2500,
#' rand.int = FALSE, 
#' covariate.AD = c("lag1", "lag2"),
#' countsonly = TRUE,
#' ID.spacetime = c("week", "state"))
#' 
#' ### obtain posterior mean/sd and credible intervals 
#' ### for observed counts across week 20-40 
#' ### at states CA and NY 
#' sum.counts(counts.post = fit$pred.counts,
#' ID.aggre = "state", timeperiod = 20:40, locs = c("CA", "NY"),
#' ID.time = "week", ID.loc = "state")
#' 
#' @export
### A function to provide summarized results of predicted counts ###
sumcounts <- function(counts.post, ID.time, ID.loc,
                      ID.aggre, timeperiod, locs){
  ## INPUTs:
  # counts.post: posterior samples (burn-in samples are discarded) of
  ## predicted counts obtained from applying functions fit.NB.st/fit.NB.timeseries
  # ID.aggre: a string indicates the variable that posterior samples are aggregated by
  # timeperiod: time-period of interest
  # locs: locations of interest
  # ID.time: a string indicates the time indicator of posterior samples
  # ID.loc: a string indicates the space indicator of posterior samples

  ID.all <- c(ID.time, ID.loc)
  col.names.post <- colnames(counts.post)
  if(!(ID.aggre %in% col.names.post)){
    stop("aggregation index must be included in the dataframe containing posterior samples")
  }
  if(!all(ID.all %in% col.names.post)){
    stop("variables indicating space and time of posterior samples must be included in the dataframe")
  }
  if(!(ID.aggre %in% ID.all)){
    stop("aggregation index must be one of space/time indicators")
  }
  index1 = counts.post[,ID.time] %in% timeperiod
  index2 = counts.post[,ID.loc] %in% locs
  post.sub <- subset(counts.post, index1 & index2)

  ID.notaggre <- ID.all[!(ID.all %in% ID.aggre)]
  post.sub1 <- post.sub %>% dplyr::select(-ID.notaggre) %>%
    group_by(across(ID.aggre)) %>% summarise_all(sum)

  get.sumres <- function(df){
    re = data.frame(est = rowMeans(df),
                    sd = apply(df, 1, sd),
                    lci = apply(df, 1, quantile, 0.025),
                    uci = apply(df, 1, quantile, 0.975))
    return(re)
  }

  post.sub2 <- post.sub1[,!colnames(post.sub1) %in% ID.aggre]
  re0 <- get.sumres(post.sub2)
  return(data.frame(post.sub1[,ID.aggre], re0))
}
