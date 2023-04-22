## A function to fit additive CAR for health data (death counts) ##

get.pg.parm2.s3 <- function(offset.vec,
                            X.design.fix.0, beta.mat,
                            beta.fix.names.0, nbeta.fix.0,
                            X.design.rand.slop, beta.rand.list, nvar.rand,
                            curr.index){

  comp.fix =  as.numeric(X.design.fix.0 %*%
                           matrix(beta.mat[curr.index, beta.fix.names.0],
                                  ncol = nbeta.fix.0))
  comp.rand.slop.list = lapply(1:nvar.rand,
                               function(x)
                                 as.numeric(X.design.rand.slop[[x]] %*%
                    matrix(beta.rand.list[[x]][curr.index, ], ncol = 1)))

  if(!is.null(offset.vec)){
    offset.vec <- as.numeric(offset.vec)
    re <- offset.vec + comp.fix + Reduce("+", comp.rand.slop.list)
  }else{
    re <- comp.fix + Reduce("+", comp.rand.slop.list)
  }
  return(re)
}

get.pg.parm2 <- function(offset.vec,
                         X.design.fix.0, beta.mat,
                         beta.fix.names.0, nbeta.fix.0,
                         X.design.rand.int, e.mat,
                         X.design.rand.slop, beta.rand.list, nvar.rand,
                         curr.index){

  comp.fix =  as.numeric(X.design.fix.0 %*%
                           matrix(beta.mat[curr.index, beta.fix.names.0],
                                  ncol = nbeta.fix.0))
  comp.rand.int = as.numeric(X.design.rand.int %*%
                               matrix(e.mat[curr.index, ], ncol = 1))
  comp.rand.slop.list = lapply(1:nvar.rand,
                               function(x)
                                 as.numeric(X.design.rand.slop[[x]] %*%
                                              matrix(beta.rand.list[[x]][curr.index, ], ncol = 1)))

  if(!is.null(offset.vec)){
    offset.vec <- as.numeric(offset.vec)
    re <- offset.vec + comp.fix + comp.rand.int + Reduce("+", comp.rand.slop.list)
  }else{
    re <- comp.fix + comp.rand.int + Reduce("+", comp.rand.slop.list)
  }
  return(re)
}


get.randcovariate.list <- function(X.design.rand.slop,
                                   beta.rand.list,
                                   curr.index){
  num.ele = length(beta.rand.list)
  comp.rand.slop.list = lapply(1:num.ele,
                               function(x)
                                 as.numeric(X.design.rand.slop[[x]] %*%
                                              matrix(beta.rand.list[[x]][curr.index, ], ncol = 1)))
  # re <- Reduce("+", comp.rand.slop.list)
  return(comp.rand.slop.list)
}

compute.adjmat <- function(ntimes, str){
  if(str == "AR1"){
    adj.list = as.matrix(cbind(c(1:(ntimes-1)), c(2:ntimes)))
    W.t = igraph::get.adjacency(
      igraph::graph.edgelist(adj.list, directed=FALSE))
  }else if(str == "AR2"){
    adj.list <- as.matrix(rbind(cbind(c(1:(ntimes-1)), c(2:ntimes)),
                                cbind(c(1:(ntimes-2)), c(3:ntimes))))
    adj.list <- adj.list[order(adj.list[,1]), ]
    W.t = igraph::get.adjacency(
      igraph::graph.edgelist(adj.list, directed=FALSE))
  }
  return(W.t)
}


### random intercept (baseline temporal trend) and constant slopes ###
fit.NB.time.s2 <- function(outcome.vec,
                           X.design.fix,
                           X.design.rand.int,
                           offset.vec = NULL,
                           rand.int.str,
                           ntimes,
                           niter, burn_in,
                           covariate.AD,
                           countsonly,
                           ID.spacetime.df){
  ## INPUTs:
  ## niter (# MCMC iterations), nstates (# locations), nweeks (# time points)

  num.obs <- length(outcome.vec)
  nbeta.fix <- ncol(X.design.fix)
  outcome.vec <- as.numeric(outcome.vec)

  # ------------ priors ----------- #
  # priors #
  ## beta
  if(nbeta.fix == 1){
    c.beta0 = matrix(rep(0, nbeta.fix), ncol = 1)
    C.beta0.pre = matrix(1/100, ncol = 1)
  }else{
    c.beta0 = matrix(rep(0, nbeta.fix), ncol = 1)
    C.beta0.pre = diag(rep(1/100, nbeta.fix))
  }
  ## tau2.t, inverse-Gamma
  a1 = 0.1; b1 = 0.1

  # tuning parameters #
  xi.tun = 0.2

  ## initialize data frame to store posterior samples ##
  omega.mat <- matrix(NA, ncol = num.obs, nrow = niter + 1)
  beta.mat <- matrix(NA, ncol = nbeta.fix, nrow = niter + 1)
  e.mat <- matrix(NA, ncol = ntimes, nrow = niter + 1)
  beta.mat[1, ] <- rep(0, nbeta.fix)
  e.mat[1, ] <- rnorm(ntimes, 0, 1)

  if(rand.int.str == "exch"){
    parm.mat <- as.data.frame(matrix(NA, ncol = 2, nrow = niter + 1))
    colnames(parm.mat) <- c("tau2.t", "xi")
    parm.mat[1, ] <- c(1, round(mean(outcome.vec)/2))
  }else if(rand.int.str %in% c("AR1", "AR2")){
    parm.mat <- as.data.frame(matrix(NA, ncol = 3, nrow = niter + 1))
    colnames(parm.mat) <- c("tau2.t", "rho.t", "xi")
    parm.mat[1, ] <- c(1, 0.5, round(mean(outcome.vec)/2))
  }

  ## initial values ##
  accept.mat <- matrix(NA, ncol = 1, nrow = niter + 1)

  if(!is.null(offset.vec)){
    pg.parm1 <- outcome.vec + parm.mat[1, "xi"]
    offset.vec <- as.numeric(offset.vec)
    pg.parm2 <- offset.vec + X.design.fix%*%matrix(beta.mat[1, ], ncol = 1) +
      X.design.rand.int %*% matrix(e.mat[1, ], ncol = 1)
    pg.parm2 <- as.numeric(pg.parm2)
    ntotal = num.obs
    omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  }else{
    pg.parm1 <- outcome.vec + parm.mat[1, "xi"]
    pg.parm2 <- X.design.fix%*%matrix(beta.mat[1, ], ncol = 1) +
      X.design.rand.int %*% matrix(e.mat[1, ], ncol = 1)
    pg.parm2 <- as.numeric(pg.parm2)
    ntotal = num.obs
    omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  }

  if(rand.int.str == "AR1"){
    adj.list = as.matrix(cbind(c(1:(ntimes-1)), c(2:ntimes)))
    W.t = get.adjacency(graph.edgelist(adj.list, directed=FALSE))
    Dw.t <- Diagonal(x = apply(W.t, 1, sum))

    lambda.t = eigen(solve(Dw.t)%*%W.t, only.values = TRUE)$values
    rho.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho.1 = sapply(rho.prior.val, function(x) 0.5*sum(log(1-x*lambda.t)), simplify = TRUE)
  }else if(rand.int.str == "AR2"){
    adj.list <- as.matrix(rbind(cbind(c(1:(ntimes-1)), c(2:ntimes)),
                                cbind(c(1:(ntimes-2)), c(3:ntimes))))
    adj.list <- adj.list[order(adj.list[,1]), ]
    W.t = get.adjacency(graph.edgelist(adj.list, directed=FALSE))
    Dw.t <- Diagonal(x = apply(W.t, 1, sum))

    lambda.t = eigen(solve(Dw.t)%*%W.t, only.values = TRUE)$values
    rho.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho.1 = sapply(rho.prior.val, function(x) 0.5*sum(log(1-x*lambda.t)), simplify = TRUE)
  }

  for(i in 1:niter){
    # 1. update omega
    omega.vec <- omega.mat[i, ]
    Omega <- Diagonal(x = omega.vec)

    # 2. update beta (conjugate)
    ## compute latent normal variables ##
    z.vec <- (outcome.vec - parm.mat[i, "xi"])/(2*omega.vec)
    ## posterior mean and variance of beta ##
    M <- solve(crossprod(X.design.fix*sqrt(omega.vec)) + C.beta0.pre)
    e.vec.i <- as.numeric(X.design.rand.int %*% matrix(e.mat[i, ], ncol = 1))
    if(!is.null(offset.vec)){
      res.vec.i <- z.vec - e.vec.i - offset.vec
    }else{
      res.vec.i <- z.vec - e.vec.i
    }
    m <- M%*%(C.beta0.pre%*%c.beta0 +
                t(sqrt(omega.vec)*X.design.fix)%*%
                (sqrt(omega.vec)*matrix(res.vec.i, ncol = 1)))
    beta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )

    # 3. update e (temporal random effects)
    ## posterior mean and variance of e ##
    fix.eff <- X.design.fix %*% matrix(beta.mat[i+1, ], ncol = 1)
    if(rand.int.str == "exch"){
      pre.t.0 <- Diagonal(x = rep((1/parm.mat[i, "tau2.t"]), ntimes))
      pre.t.post <- t(X.design.rand.int)%*%Omega%*%X.design.rand.int + pre.t.0
      Sigma.t.post <- solve(pre.t.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff  - offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff
      }
      mu.t.post <- Sigma.t.post%*%(t(X.design.rand.int)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      e.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.t.post, Sigma = Sigma.t.post) )

      # 4. update tau2.t
      parm.mat[i+1, "tau2.t"] <- 1/rgamma(1, ntimes/2 + a1,
                                          b1 + as.numeric(sum(e.mat[i+1, ]^2))/2)

    }else if(rand.int.str %in% c("AR1", "AR2")){

      pre.t.0 <- (1/parm.mat[i, "tau2.t"])*(Dw.t - parm.mat[i, "rho.t"]*W.t)
      pre.t.post <- t(X.design.rand.int)%*%Omega%*%X.design.rand.int + pre.t.0
      Sigma.t.post <- solve(pre.t.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff  - offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff
      }
      mu.t.post <- Sigma.t.post%*%(t(X.design.rand.int)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      e.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.t.post, Sigma = Sigma.t.post) )

      # 4. update tau22
      parm.mat[i+1, "tau2.t"] <-
        1/rgamma(1, ntimes/2 + a1,
                 b1 + as.numeric(t(matrix(e.mat[i+1, ], ncol = 1))%*%
                                   (Dw.t - parm.mat[i, "rho.t"]*W.t)%*%
                                   matrix(e.mat[i+1, ], ncol = 1))/2)


      # 5. update rho.t (M-H)
      inter = as.numeric( t(matrix(e.mat[i+1, ], ncol = 1))%*%
                            W.t%*%matrix(e.mat[i+1, ], ncol = 1) )
      ll.rho = ll.rho.1 + rho.prior.val/(2*parm.mat[i+1, "tau2.t"]^2)*inter
      parm.mat[i+1, "rho.t"] <- sample(x = rho.prior.val, size = 1,
                                       prob = exp(ll.rho - max(ll.rho)))


    }

    # 9. update xi (over dispersion)
    if(!is.null(offset.vec)){
      eta.vec = offset.vec + X.design.fix%*%matrix(beta.mat[i+1, ], ncol = 1) +
        as.numeric(X.design.rand.int%*%e.mat[i+1, ], ncol = 1)
      eta.vec <- as.numeric(eta.vec)
    }else{
      eta.vec = X.design.fix%*%matrix(beta.mat[i+1, ], ncol = 1) +
        as.numeric(X.design.rand.int%*%e.mat[i+1, ], ncol = 1)
      eta.vec <- as.numeric(eta.vec)
    }
    q.nb = 1/(1+exp(eta.vec)) # 1 - Pr(success)
    xi.star = rtruncnorm(1, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun)

    ll.star = sum(dnbinom(outcome.vec, size = xi.star, prob = q.nb, log = TRUE))
    ll.curr = sum(dnbinom(outcome.vec, size = parm.mat[i, "xi"], prob = q.nb, log = TRUE))
    q.star = log( dtruncnorm(parm.mat[i, "xi"], a = 0, mean = xi.star, sd = xi.tun) )
    q.curr = log( dtruncnorm(xi.star, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun) )

    ratio = min(1, exp(ll.star + q.star - ll.curr - q.curr))
    if(ratio > runif(1)){
      parm.mat[i+1, "xi"] <- xi.star
      accept.mat[i+1, 1] = 1
    }else{
      parm.mat[i+1, "xi"] <- parm.mat[i, "xi"]
      accept.mat[i+1, 1] = 0
    }

    omega.mat[i+1, ] <- rpg(ntotal,
                            outcome.vec + parm.mat[i+1, "xi"], eta.vec)

    ## tuning parameters ##
    if(i <= 4000 & i%%100 == 0){
      accept.rate <- mean(accept.mat[2:i, 1])
      if(accept.rate > 0.6){xi.tun = 1.1*xi.tun}
      if(accept.rate < 0.2){xi.tun = 0.8*xi.tun}
    }

  } # end MCMC

  ## compute WAIC ##
  ## WAIC = -2*(lppd - pWAIC);
  ## lppd (log pointwise predictive density);
  ## pWAIC (effective number of parameters)
  parm.post.comb <- cbind(beta.mat[-c(1:burn_in), ],
                          e.mat[-c(1:burn_in), ])
  X.design.comb <- cbind(X.design.fix,
                         X.design.rand.int)
  X.design.comb <- as.matrix(X.design.comb)
  pred.val.mat <- get_predvals(parm_post = parm.post.comb,
                               X_all = X.design.comb,
                               ntotal = ntotal)

  compq <- function(pred.val.vec, offset.vec){
    if(!is.null(offset.vec)){
      q = 1/(1+exp(offset.vec + pred.val.vec))
    }else{
      q = 1/(1+exp(pred.val.vec))
    }
    return(q)
  }
  q.list <- lapply(1:(ncol(pred.val.mat)),
                   function(y) compq(pred.val.vec = pred.val.mat[,y],
                                     offset.vec = offset.vec))
  q.vec <- do.call(cbind, q.list)
  compll <- function(x, xi.post){
    ll.i = dnbinom(x[1], size = xi.post, prob = x[-1], log = TRUE)
    return(ll.i)
  }
  log.pd = t( apply(cbind(outcome.vec, q.vec), 1,
                    compll, xi.post = parm.mat[-c(1:burn_in),"xi"]) )

  lppd = sum(log(apply(exp(log.pd), 1, mean)))
  pWAIC = sum(apply(log.pd, 1, var))

  X.design.fix.null <- X.design.fix
  X.design.fix.null[,covariate.AD] <- 0
  X.design.null <- cbind(X.design.fix.null,
                         X.design.rand.int)
  X.design.null <- as.matrix(X.design.null)
  pred.val.null.mat <- get_predvals(parm_post = parm.post.comb,
                                    X_all = X.design.null,
                                    ntotal = ntotal)
  
  xi.post = parm.mat[-c(1:burn_in),"xi"]
  
  if(!is.null(offset.vec)){
    
    pred.counts.mat <- get_count(offset = offset.vec,
                                 pred_mat = pred.val.mat,
                                 overdisp_post =
                                   as.numeric(parm.mat[-c(1:burn_in), "xi"]))
    pred.counts.null.mat <- get_count(offset = offset.vec,
                                      pred_mat = pred.val.null.mat,
                                      overdisp_post =
                                        as.numeric(parm.mat[-c(1:burn_in), "xi"]))
    
  }else{
    
    pred.counts.mat <- get_count(offset = rep(0, ntotal),
                                 pred_mat = pred.val.mat,
                                 overdisp_post =
                                   as.numeric(parm.mat[-c(1:burn_in), "xi"]))
    pred.counts.null.mat <- get_count(offset = rep(0, ntotal),
                                      pred_mat = pred.val.null.mat,
                                      overdisp_post =
                                        as.numeric(parm.mat[-c(1:burn_in), "xi"]))
    
  }
  
  ## get predictive distribution of outcomes and baseline counts ###
  pred.counts.nb.mat <- 
    do.call(cbind, lapply(1:length(xi.post), 
                          function(x) 
                            stats::rnbinom(n = ntotal, 
                                    size = xi.post[x], 
                                    mu = pred.counts.mat[,x])))
  pred.counts.null.nb.mat <- 
    do.call(cbind, lapply(1:length(xi.post), 
                          function(x) 
                            stats::rnbinom(n = ntotal, 
                                    size = xi.post[x], 
                                    mu = pred.counts.null.mat[,x])))
  
  AD.mat <- pred.counts.nb.mat - pred.counts.null.nb.mat
  if(!is.null(ID.spacetime.df)){
    AD.mat <- data.frame(ID.spacetime.df, AD.mat)
    pred.counts.nb.mat <- data.frame(ID.spacetime.df, 
                                     pred.counts.nb.mat)
    pred.counts.null.nb.mat <- data.frame(ID.spacetime.df, 
                                          pred.counts.null.nb.mat)
  }

  if(countsonly){
    re <- list(beta = beta.mat,
               parm = parm.mat,
               pred.counts = pred.counts.nb.mat,
               pred.null.counts = pred.counts.null.nb.mat,
               AD = AD.mat,
               WAIC = -2*(lppd - pWAIC),
               lppd = lppd, pWAIC = pWAIC)
  }else{
    re <- list(beta = beta.mat,
               parm = parm.mat,
               rand.int.t = e.mat,
               pred.counts = pred.counts.nb.mat,
               pred.null.counts = pred.counts.null.nb.mat,
               AD = AD.mat,
               accept.rate = mean(accept.mat[-c(1:burn_in), 1]),
               WAIC = -2*(lppd - pWAIC),
               lppd = lppd, pWAIC = pWAIC)
  }

  return(re)

}

### NO random intercept (baseline temporal trend) and NO time-varying slopes ###
### no random intercept and random slope ###
fit.NB.time.s4 <- function(outcome.vec,
                           X.design.fix,
                           offset.vec = NULL,
                           niter, burn_in,
                           covariate.AD,
                           countsonly,
                           ID.spacetime.df){
  ## INPUTs:
  ## niter (# MCMC iterations), nstates (# locations), nweeks (# time points)

  num.obs <- length(outcome.vec)
  nbeta.fix <- ncol(X.design.fix)

  # ------------ priors ----------- #
  # priors #
  ## beta
  if(nbeta.fix == 1){
    c.beta0 = matrix(rep(0, nbeta.fix), ncol = 1)
    C.beta0.pre = matrix(1/100, ncol = 1)
  }else{
    c.beta0 = matrix(rep(0, nbeta.fix), ncol = 1)
    C.beta0.pre = diag(rep(1/100, nbeta.fix))
  }

  # tuning parameters #
  xi.tun = 0.2

  ## initialize data frame to store posterior samples ##
  omega.mat <- matrix(NA, ncol = num.obs, nrow = niter + 1)
  beta.mat <- matrix(NA, ncol = nbeta.fix, nrow = niter + 1)
  parm.mat <- as.data.frame(matrix(NA, ncol = 1, nrow = niter + 1))
  colnames(parm.mat) <- c("xi")

  ## initial values ##
  accept.mat <- matrix(NA, ncol = 1, nrow = niter + 1)

  beta.mat[1, ] <- rep(0, nbeta.fix)
  parm.mat[1, ] <- c(round(mean(outcome.vec)/2))

  pg.parm1 <- as.numeric(outcome.vec) + parm.mat[1, "xi"]
  if(is.null(offset.vec)){
    pg.parm2 <- X.design.fix%*%matrix(beta.mat[1,], ncol = 1)
  }else{
    pg.parm2 <- offset.vec + X.design.fix%*%matrix(beta.mat[1,], ncol = 1)
  }

  ntotal = num.obs
  omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)

  for(i in 1:niter){
    # 1. update omega
    omega.vec <- omega.mat[i, ]
    Omega <- Diagonal(x = omega.vec)

    # 2. update beta (conjugate)
    z.vec <- (outcome.vec - parm.mat[i, "xi"])/(2*omega.vec)
    # 2.1 fixed beta
    ## posterior mean and variance of beta0 ##
    M <- solve(crossprod(X.design.fix*sqrt(omega.vec)) + C.beta0.pre)
    if(!is.null(offset.vec)){
      res.beta0 <- matrix(z.vec - offset.vec, ncol = 1)
    }else{
      res.beta0 <- matrix(z.vec, ncol = 1)
    }
    m <- M%*%(C.beta0.pre%*%c.beta0 +
                t(sqrt(omega.vec)*X.design.fix)%*%(sqrt(omega.vec)*res.beta0))
    beta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )

    # 3. update xi (over dispersion)
    if(!is.null(offset.vec)){
      eta.vec.i1 = offset.vec + X.design.fix%*%matrix(c(beta.mat[i+1, ]), ncol = 1)
    }else{
      eta.vec.i1 = X.design.fix%*%matrix(c(beta.mat[i+1, ]), ncol = 1)
    }
    q.nb = 1/(1+exp(eta.vec.i1)) # 1 - Pr(success)
    xi.star = rtruncnorm(1, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun)

    ll.star = sum(dnbinom(outcome.vec, size = xi.star, prob = q.nb, log = TRUE))
    ll.curr = sum(dnbinom(outcome.vec, size = parm.mat[i, "xi"], prob = q.nb, log = TRUE))
    q.star = log( dtruncnorm(parm.mat[i, "xi"], a = 0, mean = xi.star, sd = xi.tun) )
    q.curr = log( dtruncnorm(xi.star, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun) )

    ratio = min(1, exp(ll.star + q.star - ll.curr - q.curr))
    if(ratio > runif(1)){
      parm.mat[i+1, "xi"] <- xi.star
      accept.mat[i+1, 1] = 1
    }else{
      parm.mat[i+1, "xi"] <- parm.mat[i, "xi"]
      accept.mat[i+1, 1] = 0
    }

    omega.mat[i+1, ] <- rpg(num.obs, outcome.vec + parm.mat[i+1, "xi"], eta.vec.i1)

    ## tuning parameters ##
    if(i <= 4000 & i%%100 == 0){
      accept.rate <- mean(accept.mat[2:i, 1])
      if(accept.rate > 0.6){xi.tun = 1.1*xi.tun}
      if(accept.rate < 0.2){xi.tun = 0.8*xi.tun}
    }

  } # end MCMC

  ## compute WAIC ##
  ## WAIC = -2*(lppd - pWAIC);
  ## lppd (log pointwise predictive density);
  ## pWAIC (effective number of parameters)
  pred.val.mat <- get_predvals(parm_post =
                                 as.matrix(beta.mat[-c(1:burn_in), ]),
                               X_all =
                                 as.matrix(X.design.fix),
                               ntotal = ntotal)
  compq <- function(pred.val.vec, offset.vec){
    if(!is.null(offset.vec)){
      q = 1/(1+exp(offset.vec + pred.val.vec))
    }else{
      q = 1/(1+exp(pred.val.vec))
    }
    return(q)
  }
  q.list <- lapply(1:(ncol(pred.val.mat)),
                   function(y) compq(pred.val.vec = pred.val.mat[,y],
                                     offset.vec = offset.vec))
  q.vec <- do.call(cbind, q.list)
  compll <- function(x, xi.post){
    ll.i = dnbinom(x[1], size = xi.post, prob = x[-1], log = TRUE)
    return(ll.i)
  }
  log.pd = t( apply(cbind(outcome.vec, q.vec), 1,
                    compll, xi.post = parm.mat[-c(1:burn_in),"xi"]) )

  lppd = sum(log(apply(exp(log.pd), 1, mean)))
  pWAIC = sum(apply(log.pd, 1, var))

  X.design.fix.null <- X.design.fix
  X.design.fix.null[,covariate.AD] <- 0
  X.design.null <- as.matrix(X.design.fix.null)
  pred.val.null.mat <- get_predvals(parm_post =
                                      as.matrix(beta.mat[-c(1:burn_in), ]),
                                    X_all = X.design.null,
                                    ntotal = ntotal)

  if(!is.null(offset.vec)){
    pred.counts.mat <- get_count(offset = offset.vec,
                                 pred_mat = pred.val.mat,
                                 overdisp_post =
                                   as.numeric(parm.mat[-c(1:burn_in), "xi"]))
    pred.counts.null.mat <- get_count(offset = offset.vec,
                                      pred_mat = pred.val.null.mat,
                                      overdisp_post =
                                        as.numeric(parm.mat[-c(1:burn_in), "xi"]))
  }else{
    pred.counts.mat <- get_count(offset = rep(0, ntotal),
                                 pred_mat = pred.val.mat,
                                 overdisp_post =
                                   as.numeric(parm.mat[-c(1:burn_in), "xi"]))
    pred.counts.null.mat <- get_count(offset = rep(0, ntotal),
                                      pred_mat = pred.val.null.mat,
                                      overdisp_post =
                                        as.numeric(parm.mat[-c(1:burn_in), "xi"]))
  }
  
  xi.post = parm.mat[-c(1:burn_in),"xi"]
  ## get predictive distribution of outcomes and baseline counts ###
  pred.counts.nb.mat <- 
    do.call(cbind, lapply(1:length(xi.post), 
                          function(x) 
                            stats::rnbinom(n = ntotal, 
                                           size = xi.post[x], 
                                           mu = pred.counts.mat[,x])))
  pred.counts.null.nb.mat <- 
    do.call(cbind, lapply(1:length(xi.post), 
                          function(x) 
                            stats::rnbinom(n = ntotal, 
                                           size = xi.post[x], 
                                           mu = pred.counts.null.mat[,x])))
  
  AD.mat <- pred.counts.nb.mat - pred.counts.null.nb.mat
  if(!is.null(ID.spacetime.df)){
    AD.mat <- data.frame(ID.spacetime.df, AD.mat)
    pred.counts.nb.mat <- data.frame(ID.spacetime.df, 
                                     pred.counts.nb.mat)
    pred.counts.null.nb.mat <- data.frame(ID.spacetime.df, 
                                          pred.counts.null.nb.mat)
  }
  
  if(countsonly){
    re <- list(beta = beta.mat,
               parm = parm.mat,
               pred.counts = pred.counts.nb.mat,
               pred.null.counts = pred.counts.null.nb.mat,
               AD = AD.mat,
               WAIC = -2*(lppd - pWAIC),
               lppd = lppd, pWAIC = pWAIC)
  }else{
    re <- list(beta = beta.mat,
               parm = parm.mat,
               pred.counts = pred.counts.nb.mat,
               pred.null.counts = pred.counts.null.nb.mat,
               AD = AD.mat,
               accept.rate = mean(accept.mat[-c(1:burn_in), 1]),
               WAIC = -2*(lppd - pWAIC),
               lppd = lppd, pWAIC = pWAIC)
  }

  return(re)

}


