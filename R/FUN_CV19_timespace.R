## A function to fit additive CAR for spatial-temporal health data (death counts) ##
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
    W.t = get.adjacency(graph.edgelist(adj.list, directed=FALSE))
  }else if(str == "AR2"){
    adj.list <- as.matrix(rbind(cbind(c(1:(ntimes-1)), c(2:ntimes)),
                                cbind(c(1:(ntimes-2)), c(3:ntimes))))
    adj.list <- adj.list[order(adj.list[,1]), ]
    W.t = get.adjacency(graph.edgelist(adj.list, directed=FALSE))
  }
  return(W.t)
}


### space-time additive model
### time-varying coefficients and time-specific random effects ###
fit.NB.st.add.s1 <- function(outcome.vec,
                             X.design.fix, 
                             X.design.rand.int,
                             X.design.rand.slop,
                             offset.vec = NULL,
                             rand.int.str,
                             rand.slop.str,
                             niter, burn_in,
                             covariate.AD,
                             countsonly,
                             X.design.spatial,
                             spatial.str, W.s = NULL, 
                             ID.spacetime.df){
  ## INPUTs:
  ## niter (# MCMC iterations), nstates (# locations), nweeks (# time points)
  
  num.obs <- length(outcome.vec)
  nbeta.fix <- ncol(X.design.fix)
  beta.fix.names <- colnames(X.design.fix)
  
  outcome.vec <- as.numeric(outcome.vec)
  
  nbeta.rand.vec <- sapply(X.design.rand.slop, ncol)
  beta.rand.names <- names(X.design.rand.slop)
  
  
  # ------------ priors ----------- #
  # priors #
  ## beta 
  beta.fix.names.0 <- beta.fix.names[! beta.fix.names %in% beta.rand.names]
  nbeta.fix.0 <- length(beta.fix.names.0)
  if(nbeta.fix.0 == 1){
    c.beta0 = matrix(0, ncol = 1)
    C.beta0.pre = matrix(1/100, ncol = 1)
  }else{
    c.beta0 = matrix(rep(0, nbeta.fix.0), ncol = 1)
    C.beta0.pre = diag(rep(1/100, nbeta.fix.0))
  }
  
  C.beta.bar <- matrix(1/100, ncol = 1)
  c.beta.bar <- matrix(0, ncol = 1)
  
  ## tau2.t, inverse-Gamma 
  a1 = 0.1; b1 = 0.1
  
  # tuning parameters #
  xi.tun = 0.2
  
  ## initialize data frame to store posterior samples ##
  omega.mat <- matrix(NA, ncol = num.obs, nrow = niter + 1)
  
  beta.mat <- matrix(NA, ncol = nbeta.fix, nrow = niter + 1)
  colnames(beta.mat) <- beta.fix.names
  beta.mat[1, ] <- rep(0, nbeta.fix)
  
  ntimes.int <- ncol(X.design.rand.int)
  e.mat <- matrix(NA, ncol = ntimes.int, nrow = niter + 1)
  e.mat[1, ] <- rnorm(ntimes.int, 0, 1)
  
  ntimes.slop.vec <- as.numeric(nbeta.rand.vec)
  beta.rand.list <- list()
  parm.beta.list <- list()
  W.t.slop <- Dw.t.slop <- lambda.t.slop <- ll.rho.slop <- list()
  nvar.rand <- length(nbeta.rand.vec)
  for(lll in 1:nvar.rand){
    inter.lll <- matrix(NA, nrow = niter + 1,
                        ncol = nbeta.rand.vec[lll])
    inter.lll[1, ] <- rnorm(nbeta.rand.vec[lll], 0, 1)
    beta.rand.list[[lll]] <- inter.lll
    
    if(rand.slop.str[lll] == "exch"){
      parm.beta.list[[lll]] <- as.data.frame(matrix(NA, ncol = 1, 
                                                    nrow = niter + 1))
      colnames(parm.beta.list[[lll]]) <- c(paste0("sigma2.", 
                                                  beta.rand.names[lll]))
      parm.beta.list[[lll]][1, ] <- 1
    }else if(rand.slop.str[lll] %in% c("AR1", "AR2")){
      
      parm.beta.list[[lll]] <- as.data.frame(matrix(NA, ncol = 2, 
                                                    nrow = niter + 1))
      colnames(parm.beta.list[[lll]]) <- c(paste0(c("sigma2.", "rho."), 
                                                  beta.rand.names[lll]))
      parm.beta.list[[lll]][1, ] <- c(1, 0.5)
      
      W.t.slop[[lll]] = compute.adjmat(ntimes = ntimes.slop.vec[lll], 
                                       str = rand.slop.str[lll])
      Dw.t.slop[[lll]] <- Diagonal(x = apply(W.t.slop[[lll]], 1, sum))
      lambda.t.slop[[lll]] = eigen(solve(Dw.t.slop[[lll]])%*%W.t.slop[[lll]], 
                                   only.values = TRUE)$values
      rho.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
      ll.rho.slop[[lll]] = sapply(rho.prior.val, 
                                  function(x) 0.5*sum(log(1-x*lambda.t.slop[[lll]])), 
                                  simplify = TRUE)
      
    }
  }
  names(beta.rand.list) <- beta.rand.names
  names(parm.beta.list) <- beta.rand.names
  # names(W.t.slop) <- names(Dw.t.slop) <- beta.rand.names
  # names(lambda.t.slop) <- names(ll.rho.slop) <- beta.rand.names
  
  if(rand.int.str == "exch"){
    parm.mat <- as.data.frame(matrix(NA, ncol = 2, nrow = niter + 1))
    colnames(parm.mat) <- c("tau2.t", "xi")
    parm.mat[1, ] <- c(1, round(mean(outcome.vec)/2))
  }else if(rand.int.str %in% c("AR1", "AR2")){
    
    parm.mat <- as.data.frame(matrix(NA, ncol = 3, nrow = niter + 1))
    colnames(parm.mat) <- c("tau2.t", "rho.t", "xi")
    parm.mat[1, ] <- c(1, 0.5, round(mean(outcome.vec)/2))
    
    W.t.int = compute.adjmat(ntimes = ntimes.int, str = rand.int.str)
    Dw.t.int <- Diagonal(x = apply(W.t.int, 1, sum))
    lambda.t.int = eigen(solve(Dw.t.int)%*%W.t.int, only.values = TRUE)$values
    rho.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho.int.1 = sapply(rho.prior.val, 
                          function(x) 0.5*sum(log(1-x*lambda.t.int)), simplify = TRUE)
    
  }
  
  nlocs = ncol(X.design.spatial)
  theta.mat <- matrix(NA, ncol = nlocs, nrow = niter + 1)
  theta.mat[1, ] <- rnorm(nlocs, 0, 1)
  if(spatial.str == "exch"){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
  }else if(spatial.str != "exch" & !is.null(W.s)){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
    parm.mat[["rho.s"]] <- c(0.5, rep(NA, niter))
    Dw.s <- Diagonal(x = apply(W.s, 1, sum))
    
    lambda.s = eigen(solve(Dw.s)%*%W.s, only.values = TRUE)$values
    rho1.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho.s.1 = sapply(rho1.prior.val, function(x) 0.5*sum(log(1-x*lambda.s)), simplify = TRUE)
    
  }
  
  ## initial values ##
  accept.mat <- matrix(NA, ncol = 1, nrow = niter + 1)
  
  X.design.fix.0 <- matrix(X.design.fix[,beta.fix.names.0], ncol = nbeta.fix.0)
  colnames(X.design.fix.0) <- beta.fix.names.0
  
  pg.parm1 <- outcome.vec + parm.mat[1, "xi"]
  
  comp.fix =  as.numeric(X.design.fix.0 %*%
                           matrix(beta.mat[1, beta.fix.names.0],
                                  ncol = nbeta.fix.0))
  comp.rand.int = as.numeric(X.design.rand.int %*%
                               matrix(e.mat[1, ], ncol = 1))
  comp.rand.slop.list = lapply(1:nvar.rand,
                               function(x)
                                 as.numeric(X.design.rand.slop[[x]] %*%
                                              matrix(beta.rand.list[[x]][1, ], ncol = 1)))
  comp.spatial = as.numeric(X.design.spatial %*%
                              matrix(theta.mat[1, ], ncol = 1))
  if(!is.null(offset.vec)){
    offset.vec <- as.numeric(offset.vec)
    pg.parm2 <- offset.vec + comp.fix + comp.rand.int + 
      Reduce("+", comp.rand.slop.list) + comp.spatial
  }else{
    pg.parm2 <- comp.fix + comp.rand.int + 
      Reduce("+", comp.rand.slop.list) + comp.spatial
  }
  pg.parm2 <- as.numeric(pg.parm2)
  ntotal = num.obs
  omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  
  for(i in 1:niter){
    # 1. update omega 
    omega.vec <- omega.mat[i, ]
    Omega <- Diagonal(x = omega.vec)
    
    # 2. update beta (conjugate)
    ## compute latent normal variables ##
    z.vec <- (outcome.vec - parm.mat[i, "xi"])/(2*omega.vec)
    
    # 2.1 update constant beta coefficients 
    ## posterior mean and variance of beta ##
    M <- solve(crossprod(X.design.fix.0*sqrt(omega.vec)) + C.beta0.pre)
    e.vec.i <- as.numeric(X.design.rand.int %*% matrix(e.mat[i, ], ncol = 1))
    theta.vec.i <- as.numeric(X.design.spatial %*% matrix(theta.mat[i, ], ncol = 1))
    cov.rand.list.i <- get.randcovariate.list(X.design.rand.slop = X.design.rand.slop,
                                              beta.rand.list = beta.rand.list,
                                              curr.index = i)
    cov.rand.vec.i <- Reduce("+", cov.rand.list.i)
    if(!is.null(offset.vec)){
      res.vec.i <- z.vec - e.vec.i - theta.vec.i - cov.rand.vec.i - offset.vec
    }else{
      res.vec.i <- z.vec - e.vec.i - theta.vec.i - cov.rand.vec.i
    }
    m <- M%*%(C.beta0.pre%*%c.beta0 + 
                t(sqrt(omega.vec)*X.design.fix.0)%*%
                (sqrt(omega.vec)*matrix(res.vec.i, ncol = 1)))
    beta.mat[i+1, beta.fix.names.0] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )
    
    fix.eff.0 <- as.numeric(X.design.fix.0 %*% 
                              matrix(beta.mat[i+1, beta.fix.names.0], ncol = 1))
    
    # 2.2 update time-varying beta coefficients  
    # 2.3 update mean of time-varying coefficients 
    # 2.4 update parameters controlling smoothness of random slopes
    update.beta <- function(X.design.beta, C.pre, c.pri, 
                            omega.vec, res.vec.i){
      M.beta <- solve(crossprod(X.design.beta*sqrt(omega.vec)) + C.pre)
      m.beta <- M.beta%*%(C.pre%*%c.pri + 
                            t(sqrt(omega.vec)*X.design.beta)%*%
                            (sqrt(omega.vec)*res.vec.i))
      return(list(M = M.beta, m = m.beta))
    }
    
    if(!is.null(offset.vec)){
      res.beta.i0 <- z.vec - fix.eff.0 - e.vec.i - theta.vec.i - offset.vec
    }else{
      res.beta.i0 <- z.vec - fix.eff.0 - e.vec.i - theta.vec.i
    }
    
    for(lll in 1:nvar.rand){
      
      ## 2.2 update time-specific coefficients
      if(rand.slop.str[lll] == "exch"){
        parm.beta.lll = as.numeric(parm.beta.list[[lll]][i, ])
        C.beta.pre.lll <- Diagonal(x = rep(1/parm.beta.lll[1], 
                                           ntimes.slop.vec[lll]))
        c.beta = matrix(rep(beta.mat[i, beta.rand.names[lll]], 
                            nbeta.rand.vec[lll]), 
                        ncol = 1)
      }else if(rand.slop.str[lll] %in% c("AR1", "AR2")){
        parm.beta.lll = as.numeric(parm.beta.list[[lll]][i, ])
        C.beta.pre.lll <- as.matrix((1/parm.beta.lll[1])*
                                      (Dw.t.slop[[lll]] - parm.beta.lll[2]*W.t.slop[[lll]]))
        c.beta = matrix(rep(beta.mat[i, beta.rand.names[lll]], 
                            nbeta.rand.vec[lll]), 
                        ncol = 1)
      }
      
      if(lll == 1 | lll == nvar.rand){
        
        if(nvar.rand == 1){
          res.beta.i <- res.beta.i0
        }else{
          index.iter <- ifelse(lll==1, i, i+1)
          cov.rand.vec.lll = 
            get.randcovariate.list(X.design.rand.slop = 
                                     X.design.rand.slop[-lll],
                                   beta.rand.list = beta.rand.list[-lll],
                                   curr.index = index.iter)
          res.beta.i <- res.beta.i0 - Reduce("+", cov.rand.vec.lll)
        }
        res.beta.i <- matrix(res.beta.i, ncol = 1)
        mM.list <- update.beta(X.design.beta = 
                                 X.design.rand.slop[[lll]], 
                               C.pre = C.beta.pre.lll, 
                               c.pri = c.beta, 
                               omega.vec = omega.vec, 
                               res.vec.i = res.beta.i)
        beta.rand.list[[lll]][i+1, ] <- as.numeric( 
          mvrnorm(1, mu = mM.list$m, Sigma = mM.list$M) )
        
      }else{
        
        cov.rand.vec.lll.0 = 
          get.randcovariate.list(X.design.rand.slop = 
                                   X.design.rand.slop[(1:(lll-1))],
                                 beta.rand.list = beta.rand.list[(1:(lll-1))],
                                 curr.index = i+1)
        cov.rand.vec.lll.1 = 
          get.randcovariate.list(X.design.rand.slop = 
                                   X.design.rand.slop[((lll+1):nvar.rand)],
                                 beta.rand.list = 
                                   beta.rand.list[((lll+1):nvar.rand)],
                                 curr.index = i)
        res.beta.i <- res.beta.i0 - 
          Reduce("+", cov.rand.vec.lll.0) - 
          Reduce("+", cov.rand.vec.lll.1)
        res.beta.i <- matrix(res.beta.i, ncol = 1)
        
        mM.list <- update.beta(X.design.beta = 
                                 X.design.rand.slop[[lll]], 
                               C.pre = C.beta.pre.lll, 
                               c.pri = c.beta, 
                               omega.vec = omega.vec, 
                               res.vec.i = res.beta.i)
        
        beta.rand.list[[lll]][i+1, ] <- as.numeric( 
          mvrnorm(1, mu = mM.list$m, Sigma = mM.list$M) )
        
      }
      
      ## 2.3 update mean of time-varying coefficients
      X.design.betabar = matrix(rep(1, as.numeric(nbeta.rand.vec[lll])), 
                                ncol = 1)
      pre.multi.beta = matrix(colSums(C.beta.pre.lll), nrow = 1)
      M <- solve(pre.multi.beta%*%X.design.betabar + C.beta.bar)
      m <- M%*%(C.beta.bar%*%c.beta.bar + 
                  pre.multi.beta%*%matrix(beta.rand.list[[lll]][i+1, ], 
                                          ncol = 1))
      beta.mat[i+1, beta.rand.names[lll]] <- 
        as.numeric( mvrnorm(1, mu = m, Sigma = M) )
      
      ## 2.4 update parameter controlling smoothness of random slopes
      ncluster = as.numeric(nbeta.rand.vec[lll])
      if(rand.slop.str[lll] == "exch"){
        res.vec.i <- as.numeric(beta.rand.list[[lll]][i+1, ] - 
                                  beta.mat[i+1, beta.rand.names[lll]])
        
        # paste0("sigma2.", beta.rand.names[lll])
        parm.beta.list[[lll]][i+1, 1] <- 
          1/rgamma(1, ncluster/2 + a1, b1 + as.numeric(sum(res.vec.i^2))/2)
      }else if(rand.slop.str[lll] %in% c("AR1", "AR2")){
        
        parm.beta.lll = as.numeric(parm.beta.list[[lll]][i, ])
        res.vec.i <- as.numeric(beta.rand.list[[lll]][i+1, ] - 
                                  beta.mat[i+1, beta.rand.names[lll]])
        
        # paste0("sigma2.", beta.rand.names[lll])
        parm.beta.list[[lll]][i+1, 1] <- 
          1/rgamma(1, ncluster/2 + a1, 
                   b1 + as.numeric(t(matrix(res.vec.i, ncol = 1))%*%
                                     (Dw.t.slop[[lll]] - parm.beta.lll[2]*W.t.slop[[lll]])%*%
                                     matrix(res.vec.i, ncol = 1))/2)
        
        # paste0("rho.", beta.rand.names[lll]) (M-H)
        inter = as.numeric( t(matrix(res.vec.i, ncol = 1))%*%
                              W.t.slop[[lll]]%*%matrix(res.vec.i, ncol = 1) )
        ll.rho = ll.rho.slop[[lll]] + rho.prior.val/(2*parm.beta.list[[lll]][i+1, 1]^2)*inter
        parm.beta.list[[lll]][i+1, 2] <- sample(x = rho.prior.val, size = 1, 
                                                prob = exp(ll.rho - max(ll.rho)))
        
      } 
      ## end 2.4
    } # end loop over time-varying coefficients 
    
    # 3. update e (temporal random effects; random intercept) 
    ## posterior mean and variance of e ##
    rand.eff.i <- get.randcovariate.list(X.design.rand.slop = 
                                           X.design.rand.slop,
                                         beta.rand.list = beta.rand.list,
                                         curr.index = i+1)
    rand.eff.i <- Reduce("+", rand.eff.i)
    
    if(rand.int.str == "exch"){
      pre.t.0 <- Diagonal(x = rep((1/parm.mat[i, "tau2.t"]), ntimes.int))
      pre.t.post <- t(X.design.rand.int)%*%Omega%*%X.design.rand.int + pre.t.0
      Sigma.t.post <- solve(pre.t.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff.0  - rand.eff.i - theta.vec.i - offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff.0 - rand.eff.i - theta.vec.i
      }
      mu.t.post <- Sigma.t.post%*%
        (t(X.design.rand.int)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      e.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.t.post, 
                                          Sigma = Sigma.t.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.t"] <- 1/rgamma(1, ntimes.int/2 + a1, 
                                          b1 + as.numeric(sum(e.mat[i+1, ]^2))/2)
      
    }else if(rand.int.str %in% c("AR1", "AR2")){
      
      pre.t.0 <- (1/parm.mat[i, "tau2.t"])*(Dw.t.int - parm.mat[i, "rho.t"]*W.t.int)
      pre.t.post <- t(X.design.rand.int)%*%Omega%*%X.design.rand.int + pre.t.0
      Sigma.t.post <- solve(pre.t.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff.0  - rand.eff.i - theta.vec.i - offset.vec
      }else{
        res.vec.i <-  z.vec - fix.eff.0  - rand.eff.i - theta.vec.i
      }
      mu.t.post <- Sigma.t.post%*%(t(X.design.rand.int)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      e.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.t.post, Sigma = Sigma.t.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.t"] <- 
        1/rgamma(1, ntimes.int/2 + a1, 
                 b1 + as.numeric(t(matrix(e.mat[i+1, ], ncol = 1))%*%
                                   (Dw.t.int - parm.mat[i, "rho.t"]*W.t.int)%*%
                                   matrix(e.mat[i+1, ], ncol = 1))/2)
      
      
      # 5. update rho.t (M-H)
      inter = as.numeric( t(matrix(e.mat[i+1, ], ncol = 1))%*%
                            W.t.int%*%matrix(e.mat[i+1, ], ncol = 1) )
      ll.rho = ll.rho.int.1 + rho.prior.val/(2*parm.mat[i+1, "tau2.t"]^2)*inter
      parm.mat[i+1, "rho.t"] <- sample(x = rho.prior.val, size = 1, 
                                       prob = exp(ll.rho - max(ll.rho)))
      
    }
    
    
    # 4. update theta (spatial random effects) 
    ## posterior mean and variance of e ##
    e.vec.i <- as.numeric(X.design.rand.int%*%e.mat[i+1, ], ncol = 1)
    if(spatial.str == "exch"){
      pre.s.0 <- Diagonal(x = rep((1/parm.mat[i, "tau2.s"]), nlocs))
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff.0  - rand.eff.i - e.vec.i - offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff.0 - rand.eff.i - e.vec.i
      }
      mu.s.post <- Sigma.s.post%*%
        (t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                          Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 1/rgamma(1, nlocs/2 + a1, 
                                          b1 + as.numeric(sum(theta.mat[i+1, ]^2))/2)
      
    }else if(spatial.str == "CAR"){
      
      pre.s.0 <- (1/parm.mat[i, "tau2.s"])*(Dw.s - parm.mat[i, "rho.s"]*W.s)
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff.0  - rand.eff.i - e.vec.i - offset.vec
      }else{
        res.vec.i <-  z.vec - fix.eff.0  - rand.eff.i - e.vec.i
      }
      mu.s.post <- Sigma.s.post%*%(t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 
        1/rgamma(1, nlocs/2 + a1, 
                 b1 + as.numeric(t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                                   (Dw.s - parm.mat[i, "rho.s"]*W.s)%*%
                                   matrix(theta.mat[i+1, ], ncol = 1))/2)
      
      
      # 5. update rho.t (M-H)
      inter = as.numeric( t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                            W.s%*%matrix(theta.mat[i+1, ], ncol = 1) )
      ll.rho = ll.rho.s.1 + rho1.prior.val/(2*parm.mat[i+1, "tau2.s"]^2)*inter
      parm.mat[i+1, "rho.s"] <- sample(x = rho1.prior.val, size = 1, 
                                       prob = exp(ll.rho - max(ll.rho)))
      
    }
    
    
    # update xi (over dispersion) 
    if(!is.null(offset.vec)){
      eta.vec = offset.vec + fix.eff.0 + rand.eff.i + 
        e.vec.i + as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
      eta.vec <- as.numeric(eta.vec)
    }else{
      eta.vec = fix.eff.0 + rand.eff.i + 
        e.vec.i + as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
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
  parm.post.comb <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                          do.call(cbind, beta.rand.list)[-c(1:burn_in), ],
                          e.mat[-c(1:burn_in), ],
                          theta.mat[-c(1:burn_in), ])
  X.design.comb <- cbind(X.design.fix.0, 
                         do.call(cbind, X.design.rand.slop),
                         X.design.rand.int,
                         X.design.spatial)
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
  
  if(all(covariate.AD %in% names(X.design.rand.slop))){
    index.rand <- !(names(X.design.rand.slop) %in% covariate.AD)
    parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                            do.call(cbind, beta.rand.list[index.rand])[-c(1:burn_in),],
                            e.mat[-c(1:burn_in), ],
                            theta.mat[-c(1:burn_in), ])
    X.design.null <- cbind(X.design.fix.0, 
                           do.call(cbind, X.design.rand.slop[index.rand]),
                           X.design.rand.int,
                           X.design.spatial)
    X.design.null <- as.matrix(X.design.null)
    pred.val.null.mat <- get_predvals(parm_post = parm.post.null,
                                      X_all = X.design.null,
                                      ntotal = ntotal)
  }else if(any(covariate.AD %in% names(X.design.rand.slop))){
    index.curr <- covariate.AD %in% names(X.design.rand.slop)
    index.fix <- covariate.AD[!index.curr]
    index.rand <- !names(X.design.rand.slop)%in%covariate.AD
    X.design.fix.null <- X.design.fix.0
    X.design.fix.null[,index.fix] <- 0
    beta.rand.comb <- beta.rand.list[index.rand]
    if(length(beta.rand.comb) == 0){
      parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                              e.mat[-c(1:burn_in), ],
                              theta.mat[-c(1:burn_in), ])
      X.design.null <- cbind(X.design.fix.0,
                             X.design.rand.int,
                             X.design.spatial)
    }else{
      parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                              do.call(cbind, beta.rand.comb)[-c(1:burn_in), ],
                              e.mat[-c(1:burn_in), ],
                              theta.mat[-c(1:burn_in), ])
      X.design.null <- cbind(X.design.fix.0,
                             do.call(cbind, 
                                     X.design.rand.slop[index.rand]),
                             X.design.rand.int,
                             X.design.spatial)
    }
    X.design.null <- as.matrix(X.design.null)
    pred.val.null.mat <- get_predvals(parm_post = parm.post.null,
                                      X_all = X.design.null,
                                      ntotal = ntotal)
  }else if(!any(covariate.AD %in% names(X.design.rand.slop))){
    X.design.fix.null <- X.design.fix.0
    X.design.fix.null[,covariate.AD] <- 0
    X.design.null <- cbind(X.design.fix.null, 
                           do.call(cbind, X.design.rand.slop),
                           X.design.rand.int,
                           X.design.spatial)
    parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                            do.call(cbind, beta.rand.list)[-c(1:burn_in), ],
                            e.mat[-c(1:burn_in), ],
                            theta.mat[-c(1:burn_in), ])
    X.design.null <- as.matrix(X.design.null)
    pred.val.null.mat <- get_predvals(parm_post = parm.post.null,
                                      X_all = X.design.null,
                                      ntotal = ntotal)
  }
  
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
               parm.beta.tvarying = parm.beta.list, 
               pred.counts = pred.counts.nb.mat,
               pred.null.counts = pred.counts.null.nb.mat,
               AD = AD.mat, 
               WAIC = -2*(lppd - pWAIC), 
               lppd = lppd, pWAIC = pWAIC)
  }else{
    re <- list(beta = beta.mat, 
               parm = parm.mat, 
               rand.int.t = e.mat, 
               rand.int.s = theta.mat,
               beta.tvarying = beta.rand.list,
               parm.beta.tvarying = parm.beta.list,
               pred.counts = pred.counts.nb.mat,
               pred.null.counts = pred.counts.null.nb.mat,
               AD = AD.mat,
               accept.rate = mean(accept.mat[-c(1:burn_in), 1]), 
               WAIC = -2*(lppd - pWAIC), 
               lppd = lppd, pWAIC = pWAIC)
  }
  return(re)
  
  
}


### space-time additive model
### random intercept (baseline temporal trend) and constant slopes ###
fit.NB.st.add.s2 <- function(outcome.vec,
                             X.design.fix, 
                             X.design.rand.int,
                             offset.vec = NULL,
                             rand.int.str,
                             niter, burn_in,
                             covariate.AD,
                             countsonly,
                             X.design.spatial,
                             spatial.str, W.s = NULL,
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
  
  ntimes = ncol(X.design.rand.int)
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
  
  nlocs = ncol(X.design.spatial)
  theta.mat <- matrix(NA, ncol = nlocs, nrow = niter + 1)
  theta.mat[1, ] <- rnorm(nlocs, 0, 1)
  if(spatial.str == "exch"){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
  }else if(spatial.str != "exch" & !is.null(W.s)){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
    parm.mat[["rho.s"]] <- c(0.5, rep(NA, niter))
    Dw.s <- Diagonal(x = apply(W.s, 1, sum))
    
    lambda.s = eigen(solve(Dw.s)%*%W.s, only.values = TRUE)$values
    rho1.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho.s.1 = sapply(rho1.prior.val, function(x) 0.5*sum(log(1-x*lambda.s)), simplify = TRUE)
    
  }
  
  ## initial values ##
  accept.mat <- matrix(NA, ncol = 1, nrow = niter + 1)
  
  if(!is.null(offset.vec)){
    pg.parm1 <- outcome.vec + parm.mat[1, "xi"]
    offset.vec <- as.numeric(offset.vec)
    pg.parm2 <- offset.vec + X.design.fix%*%matrix(beta.mat[1, ], ncol = 1) + 
      X.design.rand.int %*% matrix(e.mat[1, ], ncol = 1) + 
      X.design.spatial %*% matrix(theta.mat[1, ], ncol = 1)
    pg.parm2 <- as.numeric(pg.parm2)
    ntotal = num.obs
    omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  }else{
    pg.parm1 <- outcome.vec + parm.mat[1, "xi"]
    pg.parm2 <- X.design.fix%*%matrix(beta.mat[1, ], ncol = 1) + 
      X.design.rand.int %*% matrix(e.mat[1, ], ncol = 1) + 
      X.design.spatial %*% matrix(theta.mat[1, ], ncol = 1)
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
    theta.vec.i <- as.numeric(X.design.spatial %*% matrix(theta.mat[i, ], ncol = 1))
    if(!is.null(offset.vec)){
      res.vec.i <- z.vec - e.vec.i - theta.vec.i - offset.vec
    }else{
      res.vec.i <- z.vec - e.vec.i - theta.vec.i
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
        res.vec.i <- z.vec - fix.eff - theta.vec.i -  offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff - theta.vec.i
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
        res.vec.i <- z.vec - fix.eff - theta.vec.i - offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff - theta.vec.i
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
    
    
    # 4. update theta (spatial random effects) 
    ## posterior mean and variance of e ##
    e.vec.i <- as.numeric(X.design.rand.int%*%e.mat[i+1, ], ncol = 1)
    if(spatial.str == "exch"){
      pre.s.0 <- Diagonal(x = rep((1/parm.mat[i, "tau2.s"]), nlocs))
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff - e.vec.i - offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff - e.vec.i
      }
      mu.s.post <- Sigma.s.post%*%
        (t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 1/rgamma(1, nlocs/2 + a1, 
                                          b1 + as.numeric(sum(theta.mat[i+1, ]^2))/2)
      
    }else if(spatial.str == "CAR"){
      
      pre.s.0 <- (1/parm.mat[i, "tau2.s"])*(Dw.s - parm.mat[i, "rho.s"]*W.s)
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff  - e.vec.i - offset.vec
      }else{
        res.vec.i <-  z.vec - fix.eff  - e.vec.i
      }
      mu.s.post <- Sigma.s.post%*%(t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 
        1/rgamma(1, nlocs/2 + a1, 
                 b1 + as.numeric(t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                                   (Dw.s - parm.mat[i, "rho.s"]*W.s)%*%
                                   matrix(theta.mat[i+1, ], ncol = 1))/2)
      
      
      # 5. update rho.t (M-H)
      inter = as.numeric( t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                            W.s%*%matrix(theta.mat[i+1, ], ncol = 1) )
      ll.rho = ll.rho.s.1 + rho1.prior.val/(2*parm.mat[i+1, "tau2.s"]^2)*inter
      parm.mat[i+1, "rho.s"] <- sample(x = rho1.prior.val, size = 1, 
                                       prob = exp(ll.rho - max(ll.rho)))
      
    }
    
    
    # 9. update xi (over dispersion) 
    if(!is.null(offset.vec)){
      eta.vec = offset.vec + X.design.fix%*%matrix(beta.mat[i+1, ], ncol = 1) + 
        e.vec.i + as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
      eta.vec <- as.numeric(eta.vec)
    }else{
      eta.vec = X.design.fix%*%matrix(beta.mat[i+1, ], ncol = 1) + 
        e.vec.i + as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
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
                          e.mat[-c(1:burn_in), ],
                          theta.mat[-c(1:burn_in), ])
  X.design.comb <- cbind(X.design.fix, 
                         X.design.rand.int,
                         X.design.spatial)
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
                         X.design.rand.int,
                         X.design.spatial)
  X.design.null <- as.matrix(X.design.null)
  pred.val.null.mat <- get_predvals(parm_post = parm.post.comb,
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
               rand.int.t = e.mat,
               rand.int.s = theta.mat,
               pred.counts = pred.counts.nb.mat,
               pred.null.counts = pred.counts.null.nb.mat,
               AD = AD.mat, 
               accept.rate = mean(accept.mat[-c(1:burn_in), 1]), 
               WAIC = -2*(lppd - pWAIC), 
               lppd = lppd, pWAIC = pWAIC)
  }
  
  return(re)
  
}




### time-varying coefficients and NO time-specific random effects ###
fit.NB.st.add.s3 <- function(outcome.vec,
                             X.design.fix, 
                             X.design.rand.slop,
                             offset.vec = NULL,
                             rand.slop.str,
                             niter, burn_in,
                             covariate.AD,
                             countsonly,
                             X.design.spatial,
                             spatial.str, W.s = NULL,
                             ID.spacetime.df){
  ## INPUTs:
  ## niter (# MCMC iterations), nstates (# locations), nweeks (# time points)
  
  num.obs <- length(outcome.vec)
  nbeta.fix <- ncol(X.design.fix)
  beta.fix.names <- colnames(X.design.fix)
  
  outcome.vec <- as.numeric(outcome.vec)
  
  nbeta.rand.vec <- sapply(X.design.rand.slop, ncol)
  beta.rand.names <- names(X.design.rand.slop)
  
  
  # ------------ priors ----------- #
  # priors #
  ## beta 
  beta.fix.names.0 <- beta.fix.names[! beta.fix.names %in% beta.rand.names]
  nbeta.fix.0 <- length(beta.fix.names.0)
  if(nbeta.fix.0 == 1){
    c.beta0 = matrix(0, ncol = 1)
    C.beta0.pre = matrix(1/100, ncol = 1)
  }else{
    c.beta0 = matrix(rep(0, nbeta.fix.0), ncol = 1)
    C.beta0.pre = diag(rep(1/100, nbeta.fix.0))
  }
  
  C.beta.bar <- matrix(1/100, ncol = 1)
  c.beta.bar <- matrix(0, ncol = 1)
  
  ## tau2.t, inverse-Gamma 
  a1 = 0.1; b1 = 0.1
  
  # tuning parameters #
  xi.tun = 0.2
  
  ## initialize data frame to store posterior samples ##
  omega.mat <- matrix(NA, ncol = num.obs, nrow = niter + 1)
  
  beta.mat <- matrix(NA, ncol = nbeta.fix, nrow = niter + 1)
  colnames(beta.mat) <- beta.fix.names
  beta.mat[1, ] <- rep(0, nbeta.fix)

  ntimes.slop.vec <- as.numeric(nbeta.rand.vec)
  beta.rand.list <- list()
  parm.beta.list <- list()
  W.t.slop <- Dw.t.slop <- lambda.t.slop <- ll.rho.slop <- list()
  nvar.rand <- length(nbeta.rand.vec)
  for(lll in 1:nvar.rand){
    inter.lll <- matrix(NA, nrow = niter + 1,
                        ncol = nbeta.rand.vec[lll])
    inter.lll[1, ] <- rnorm(nbeta.rand.vec[lll], 0, 1)
    beta.rand.list[[lll]] <- inter.lll
    
    if(rand.slop.str[lll] == "exch"){
      parm.beta.list[[lll]] <- as.data.frame(matrix(NA, ncol = 1, 
                                                    nrow = niter + 1))
      colnames(parm.beta.list[[lll]]) <- c(paste0("sigma2.", 
                                                  beta.rand.names[lll]))
      parm.beta.list[[lll]][1, ] <- 1
    }else if(rand.slop.str[lll] %in% c("AR1", "AR2")){
      
      parm.beta.list[[lll]] <- as.data.frame(matrix(NA, ncol = 2, 
                                                    nrow = niter + 1))
      colnames(parm.beta.list[[lll]]) <- c(paste0(c("sigma2.", "rho."), 
                                                  beta.rand.names[lll]))
      parm.beta.list[[lll]][1, ] <- c(1, 0.5)
      
      W.t.slop[[lll]] = compute.adjmat(ntimes = ntimes.slop.vec[lll], 
                                       str = rand.slop.str[lll])
      Dw.t.slop[[lll]] <- Diagonal(x = apply(W.t.slop[[lll]], 1, sum))
      lambda.t.slop[[lll]] = eigen(solve(Dw.t.slop[[lll]])%*%W.t.slop[[lll]], 
                                   only.values = TRUE)$values
      rho.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
      ll.rho.slop[[lll]] = sapply(rho.prior.val, 
                                  function(x) 0.5*sum(log(1-x*lambda.t.slop[[lll]])), 
                                  simplify = TRUE)
      
    }
  }
  names(beta.rand.list) <- beta.rand.names
  names(parm.beta.list) <- beta.rand.names
  # names(W.t.slop) <- names(Dw.t.slop) <- beta.rand.names
  # names(lambda.t.slop) <- names(ll.rho.slop) <- beta.rand.names
  
  
  nlocs = ncol(X.design.spatial)
  theta.mat <- matrix(NA, ncol = nlocs, nrow = niter + 1)
  theta.mat[1, ] <- rnorm(nlocs, 0, 1)
  if(spatial.str == "exch"){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
  }else if(spatial.str != "exch" & !is.null(W.s)){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
    parm.mat[["rho.s"]] <- c(0.5, rep(NA, niter))
    Dw.s <- Diagonal(x = apply(W.s, 1, sum))
    
    lambda.s = eigen(solve(Dw.s)%*%W.s, only.values = TRUE)$values
    rho1.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho.s.1 = sapply(rho1.prior.val, function(x) 0.5*sum(log(1-x*lambda.s)), simplify = TRUE)
  }
  
  ## initial values ##
  accept.mat <- matrix(NA, ncol = 1, nrow = niter + 1)
  
  X.design.fix.0 <- matrix(X.design.fix[,beta.fix.names.0], ncol = nbeta.fix.0)
  colnames(X.design.fix.0) <- beta.fix.names.0
  
  pg.parm1 <- outcome.vec + parm.mat[1, "xi"]
  
  comp.fix =  as.numeric(X.design.fix.0 %*%
                           matrix(beta.mat[1, beta.fix.names.0],
                                  ncol = nbeta.fix.0))
  comp.rand.slop.list = lapply(1:nvar.rand,
                               function(x)
                                 as.numeric(X.design.rand.slop[[x]] %*%
                                              matrix(beta.rand.list[[x]][1, ], ncol = 1)))
  comp.spatial = as.numeric(X.design.spatial %*%
                              matrix(theta.mat[1, ], ncol = 1))
  if(!is.null(offset.vec)){
    offset.vec <- as.numeric(offset.vec)
    pg.parm2 <- offset.vec + comp.fix + 
      Reduce("+", comp.rand.slop.list) + comp.spatial
  }else{
    pg.parm2 <- comp.fix + 
      Reduce("+", comp.rand.slop.list) + comp.spatial
  }
  pg.parm2 <- as.numeric(pg.parm2)
  ntotal = num.obs
  omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  
  for(i in 1:niter){
    # 1. update omega 
    omega.vec <- omega.mat[i, ]
    Omega <- Diagonal(x = omega.vec)
    
    # 2. update beta (conjugate)
    ## compute latent normal variables ##
    z.vec <- (outcome.vec - parm.mat[i, "xi"])/(2*omega.vec)
    
    # 2.1 update constant beta coefficients 
    ## posterior mean and variance of beta ##
    M <- solve(crossprod(X.design.fix.0*sqrt(omega.vec)) + C.beta0.pre)
    theta.vec.i <- as.numeric(X.design.spatial %*% matrix(theta.mat[i, ], ncol = 1))
    cov.rand.list.i <- get.randcovariate.list(X.design.rand.slop = X.design.rand.slop,
                                              beta.rand.list = beta.rand.list,
                                              curr.index = i)
    cov.rand.vec.i <- Reduce("+", cov.rand.list.i)
    if(!is.null(offset.vec)){
      res.vec.i <- z.vec - theta.vec.i - cov.rand.vec.i - offset.vec
    }else{
      res.vec.i <- z.vec - theta.vec.i - cov.rand.vec.i
    }
    m <- M%*%(C.beta0.pre%*%c.beta0 + 
                t(sqrt(omega.vec)*X.design.fix.0)%*%
                (sqrt(omega.vec)*matrix(res.vec.i, ncol = 1)))
    beta.mat[i+1, beta.fix.names.0] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )
    
    fix.eff.0 <- as.numeric(X.design.fix.0 %*% 
                              matrix(beta.mat[i+1, beta.fix.names.0], ncol = 1))
    
    # 2.2 update time-varying beta coefficients  
    # 2.3 update mean of time-varying coefficients 
    # 2.4 update parameters controlling smoothness of random slopes
    update.beta <- function(X.design.beta, C.pre, c.pri, 
                            omega.vec, res.vec.i){
      M.beta <- solve(crossprod(X.design.beta*sqrt(omega.vec)) + C.pre)
      m.beta <- M.beta%*%(C.pre%*%c.pri + 
                            t(sqrt(omega.vec)*X.design.beta)%*%
                            (sqrt(omega.vec)*res.vec.i))
      return(list(M = M.beta, m = m.beta))
    }
    
    if(!is.null(offset.vec)){
      res.beta.i0 <- z.vec - fix.eff.0 - theta.vec.i - offset.vec
    }else{
      res.beta.i0 <- z.vec - fix.eff.0 - theta.vec.i
    }
    
    for(lll in 1:nvar.rand){
      
      ## 2.2 update time-specific coefficients
      if(rand.slop.str[lll] == "exch"){
        parm.beta.lll = as.numeric(parm.beta.list[[lll]][i, ])
        C.beta.pre.lll <- Diagonal(x = rep(1/parm.beta.lll[1], 
                                           ntimes.slop.vec[lll]))
        c.beta = matrix(rep(beta.mat[i, beta.rand.names[lll]], 
                            nbeta.rand.vec[lll]), 
                        ncol = 1)
      }else if(rand.slop.str[lll] %in% c("AR1", "AR2")){
        parm.beta.lll = as.numeric(parm.beta.list[[lll]][i, ])
        C.beta.pre.lll <- as.matrix((1/parm.beta.lll[1])*
                                      (Dw.t.slop[[lll]] - parm.beta.lll[2]*W.t.slop[[lll]]))
        c.beta = matrix(rep(beta.mat[i, beta.rand.names[lll]], 
                            nbeta.rand.vec[lll]), 
                        ncol = 1)
      }
      
      if(lll == 1 | lll == nvar.rand){
        
        if(nvar.rand == 1){
          res.beta.i <- res.beta.i0
        }else{
          index.iter <- ifelse(lll==1, i, i+1)
          cov.rand.vec.lll = 
            get.randcovariate.list(X.design.rand.slop = 
                                     X.design.rand.slop[-lll],
                                   beta.rand.list = beta.rand.list[-lll],
                                   curr.index = index.iter)
          res.beta.i <- res.beta.i0 - Reduce("+", cov.rand.vec.lll)
        }
        res.beta.i <- matrix(res.beta.i, ncol = 1)
        mM.list <- update.beta(X.design.beta = 
                                 X.design.rand.slop[[lll]], 
                               C.pre = C.beta.pre.lll, 
                               c.pri = c.beta, 
                               omega.vec = omega.vec, 
                               res.vec.i = res.beta.i)
        beta.rand.list[[lll]][i+1, ] <- as.numeric( 
          mvrnorm(1, mu = mM.list$m, Sigma = mM.list$M) )
        
      }else{
        
        cov.rand.vec.lll.0 = 
          get.randcovariate.list(X.design.rand.slop = 
                                   X.design.rand.slop[(1:(lll-1))],
                                 beta.rand.list = beta.rand.list[(1:(lll-1))],
                                 curr.index = i+1)
        cov.rand.vec.lll.1 = 
          get.randcovariate.list(X.design.rand.slop = 
                                   X.design.rand.slop[((lll+1):nvar.rand)],
                                 beta.rand.list = 
                                   beta.rand.list[((lll+1):nvar.rand)],
                                 curr.index = i)
        res.beta.i <- res.beta.i0 - 
          Reduce("+", cov.rand.vec.lll.0) - 
          Reduce("+", cov.rand.vec.lll.1)
        res.beta.i <- matrix(res.beta.i, ncol = 1)
        
        mM.list <- update.beta(X.design.beta = 
                                 X.design.rand.slop[[lll]], 
                               C.pre = C.beta.pre.lll, 
                               c.pri = c.beta, 
                               omega.vec = omega.vec, 
                               res.vec.i = res.beta.i)
        
        beta.rand.list[[lll]][i+1, ] <- as.numeric( 
          mvrnorm(1, mu = mM.list$m, Sigma = mM.list$M) )
        
      }
      
      ## 2.3 update mean of time-varying coefficients
      X.design.betabar = matrix(rep(1, as.numeric(nbeta.rand.vec[lll])), 
                                ncol = 1)
      pre.multi.beta = matrix(colSums(C.beta.pre.lll), nrow = 1)
      M <- solve(pre.multi.beta%*%X.design.betabar + C.beta.bar)
      m <- M%*%(C.beta.bar%*%c.beta.bar + 
                  pre.multi.beta%*%matrix(beta.rand.list[[lll]][i+1, ], 
                                          ncol = 1))
      beta.mat[i+1, beta.rand.names[lll]] <- 
        as.numeric( mvrnorm(1, mu = m, Sigma = M) )
      
      ## 2.4 update parameter controlling smoothness of random slopes
      ncluster = as.numeric(nbeta.rand.vec[lll])
      if(rand.slop.str[lll] == "exch"){
        res.vec.i <- as.numeric(beta.rand.list[[lll]][i+1, ] - 
                                  beta.mat[i+1, beta.rand.names[lll]])
        
        # paste0("sigma2.", beta.rand.names[lll])
        parm.beta.list[[lll]][i+1, 1] <- 
          1/rgamma(1, ncluster/2 + a1, b1 + as.numeric(sum(res.vec.i^2))/2)
      }else if(rand.slop.str[lll] %in% c("AR1", "AR2")){
        
        parm.beta.lll = as.numeric(parm.beta.list[[lll]][i, ])
        res.vec.i <- as.numeric(beta.rand.list[[lll]][i+1, ] - 
                                  beta.mat[i+1, beta.rand.names[lll]])
        
        # paste0("sigma2.", beta.rand.names[lll])
        parm.beta.list[[lll]][i+1, 1] <- 
          1/rgamma(1, ncluster/2 + a1, 
                   b1 + as.numeric(t(matrix(res.vec.i, ncol = 1))%*%
                                     (Dw.t.slop[[lll]] - parm.beta.lll[2]*W.t.slop[[lll]])%*%
                                     matrix(res.vec.i, ncol = 1))/2)
        
        # paste0("rho.", beta.rand.names[lll]) (M-H)
        inter = as.numeric( t(matrix(res.vec.i, ncol = 1))%*%
                              W.t.slop[[lll]]%*%matrix(res.vec.i, ncol = 1) )
        ll.rho = ll.rho.slop[[lll]] + rho.prior.val/(2*parm.beta.list[[lll]][i+1, 1]^2)*inter
        parm.beta.list[[lll]][i+1, 2] <- sample(x = rho.prior.val, size = 1, 
                                                prob = exp(ll.rho - max(ll.rho)))
        
      } 
      ## end 2.4
    } # end loop over time-varying coefficients 
    

    # 4. update theta (spatial random effects) 
    ## posterior mean and variance of e ##
    rand.eff.i <- get.randcovariate.list(X.design.rand.slop = 
                                           X.design.rand.slop,
                                         beta.rand.list = beta.rand.list,
                                         curr.index = i+1)
    rand.eff.i <- Reduce("+", rand.eff.i)
    if(spatial.str == "exch"){
      pre.s.0 <- Diagonal(x = rep((1/parm.mat[i, "tau2.s"]), nlocs))
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff.0  - rand.eff.i - offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff.0 - rand.eff.i
      }
      mu.s.post <- Sigma.s.post%*%
        (t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 1/rgamma(1, nlocs/2 + a1, 
                                          b1 + as.numeric(sum(theta.mat[i+1, ]^2))/2)
      
    }else if(spatial.str == "CAR"){
      
      pre.s.0 <- (1/parm.mat[i, "tau2.s"])*(Dw.s - parm.mat[i, "rho.s"]*W.s)
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff.0  - rand.eff.i - offset.vec
      }else{
        res.vec.i <-  z.vec - fix.eff.0  - rand.eff.i
      }
      mu.s.post <- Sigma.s.post%*%(t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 
        1/rgamma(1, nlocs/2 + a1, 
                 b1 + as.numeric(t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                                   (Dw.s - parm.mat[i, "rho.s"]*W.s)%*%
                                   matrix(theta.mat[i+1, ], ncol = 1))/2)
      
      
      # 5. update rho.t (M-H)
      inter = as.numeric( t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                            W.s%*%matrix(theta.mat[i+1, ], ncol = 1) )
      ll.rho = ll.rho.s.1 + rho1.prior.val/(2*parm.mat[i+1, "tau2.s"]^2)*inter
      parm.mat[i+1, "rho.s"] <- sample(x = rho1.prior.val, size = 1, 
                                       prob = exp(ll.rho - max(ll.rho)))
      
    }
    
    
    # update xi (over dispersion) 
    if(!is.null(offset.vec)){
      eta.vec = offset.vec + fix.eff.0 + rand.eff.i + 
        as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
      eta.vec <- as.numeric(eta.vec)
    }else{
      eta.vec = fix.eff.0 + rand.eff.i + 
        as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
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
  parm.post.comb <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                          do.call(cbind, beta.rand.list)[-c(1:burn_in), ],
                          theta.mat[-c(1:burn_in), ])
  X.design.comb <- cbind(X.design.fix.0, 
                         do.call(cbind, X.design.rand.slop),
                         X.design.spatial)
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
  
  if(all(covariate.AD %in% names(X.design.rand.slop))){
    index.rand <- !(names(X.design.rand.slop) %in% covariate.AD)
    parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                            do.call(cbind, beta.rand.list[index.rand])[-c(1:burn_in),],
                            theta.mat[-c(1:burn_in), ])
    X.design.null <- cbind(X.design.fix.0, 
                           do.call(cbind, X.design.rand.slop[index.rand]),
                           X.design.spatial)
    X.design.null <- as.matrix(X.design.null)
    pred.val.null.mat <- get_predvals(parm_post = parm.post.null,
                                      X_all = X.design.null,
                                      ntotal = ntotal)
  }else if(any(covariate.AD %in% names(X.design.rand.slop))){
    index.curr <- covariate.AD %in% names(X.design.rand.slop)
    index.fix <- covariate.AD[!index.curr]
    index.rand <- !names(X.design.rand.slop)%in%covariate.AD
    X.design.fix.null <- X.design.fix.0
    X.design.fix.null[,index.fix] <- 0
    beta.rand.comb <- beta.rand.list[index.rand]
    if(length(beta.rand.comb) == 0){
      parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                              theta.mat[-c(1:burn_in), ])
      X.design.null <- cbind(X.design.fix.0,
                             X.design.spatial)
    }else{
      parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                              do.call(cbind, beta.rand.comb)[-c(1:burn_in), ],
                              theta.mat[-c(1:burn_in), ])
      X.design.null <- cbind(X.design.fix.0,
                             do.call(cbind, 
                                     X.design.rand.slop[index.rand]),
                             X.design.spatial)
    }
    X.design.null <- as.matrix(X.design.null)
    pred.val.null.mat <- get_predvals(parm_post = parm.post.null,
                                      X_all = X.design.null,
                                      ntotal = ntotal)
  }else if(!any(covariate.AD %in% names(X.design.rand.slop))){
    X.design.fix.null <- X.design.fix.0
    X.design.fix.null[,covariate.AD] <- 0
    X.design.null <- cbind(X.design.fix.null, 
                           do.call(cbind, X.design.rand.slop),
                           X.design.spatial)
    parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                            do.call(cbind, beta.rand.list)[-c(1:burn_in), ],
                            theta.mat[-c(1:burn_in), ])
    X.design.null <- as.matrix(X.design.null)
    pred.val.null.mat <- get_predvals(parm_post = parm.post.null,
                                      X_all = X.design.null,
                                      ntotal = ntotal)
  }
  
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
               parm.beta.tvarying = parm.beta.list, 
               pred.counts = pred.counts.nb.mat,
               pred.null.counts = pred.counts.null.nb.mat,
               AD = AD.mat, 
               WAIC = -2*(lppd - pWAIC), 
               lppd = lppd, pWAIC = pWAIC)
  }else{
    re <- list(beta = beta.mat, 
               parm = parm.mat, 
               rand.int.s = theta.mat,
               beta.tvarying = beta.rand.list,
               parm.beta.tvarying = parm.beta.list,
               pred.counts = pred.counts.nb.mat,
               pred.null.counts = pred.counts.null.nb.mat,
               AD = AD.mat,
               accept.rate = mean(accept.mat[-c(1:burn_in), 1]), 
               WAIC = -2*(lppd - pWAIC), 
               lppd = lppd, pWAIC = pWAIC)
  }
  return(re)
  
  
}

### space-time additive model
### no random intercept (baseline temporal trend) and constant slopes ###
fit.NB.st.add.s4 <- function(outcome.vec,
                             X.design.fix, 
                             offset.vec = NULL,
                             niter, burn_in,
                             covariate.AD,
                             countsonly,
                             X.design.spatial,
                             spatial.str, W.s = NULL,
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
  beta.mat[1, ] <- rep(0, nbeta.fix)
  
  parm.mat <- as.data.frame(matrix(NA, ncol = 1, nrow = niter + 1))
  colnames(parm.mat) <- c("xi")
  parm.mat[1, ] <- c(round(mean(outcome.vec)/2))
  
  nlocs = ncol(X.design.spatial)
  theta.mat <- matrix(NA, ncol = nlocs, nrow = niter + 1)
  theta.mat[1, ] <- rnorm(nlocs, 0, 1)
  if(spatial.str == "exch"){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
  }else if(spatial.str != "exch" & !is.null(W.s)){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
    parm.mat[["rho.s"]] <- c(0.5, rep(NA, niter))
    Dw.s <- Diagonal(x = apply(W.s, 1, sum))
    
    lambda.s = eigen(solve(Dw.s)%*%W.s, only.values = TRUE)$values
    rho1.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho.s.1 = sapply(rho1.prior.val, function(x) 0.5*sum(log(1-x*lambda.s)), simplify = TRUE)
    
  }
  
  ## initial values ##
  accept.mat <- matrix(NA, ncol = 1, nrow = niter + 1)
  
  if(!is.null(offset.vec)){
    pg.parm1 <- outcome.vec + parm.mat[1, "xi"]
    offset.vec <- as.numeric(offset.vec)
    pg.parm2 <- offset.vec + X.design.fix%*%matrix(beta.mat[1, ], ncol = 1) + 
      X.design.spatial %*% matrix(theta.mat[1, ], ncol = 1)
    pg.parm2 <- as.numeric(pg.parm2)
    ntotal = num.obs
    omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  }else{
    pg.parm1 <- outcome.vec + parm.mat[1, "xi"]
    pg.parm2 <- X.design.fix%*%matrix(beta.mat[1, ], ncol = 1) +
      X.design.spatial %*% matrix(theta.mat[1, ], ncol = 1)
    pg.parm2 <- as.numeric(pg.parm2)
    ntotal = num.obs
    omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
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
    theta.vec.i <- as.numeric(X.design.spatial %*% 
                                matrix(theta.mat[i, ], ncol = 1))
    if(!is.null(offset.vec)){
      res.vec.i <- z.vec  - theta.vec.i - offset.vec
    }else{
      res.vec.i <- z.vec  - theta.vec.i
    }
    m <- M%*%(C.beta0.pre%*%c.beta0 + 
                t(sqrt(omega.vec)*X.design.fix)%*%
                (sqrt(omega.vec)*matrix(res.vec.i, ncol = 1)))
    beta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )
    
    fix.eff <- X.design.fix %*% matrix(beta.mat[i+1, ], ncol = 1)
    # 4. update theta (spatial random effects) 
    ## posterior mean and variance of e ##
    if(spatial.str == "exch"){
      pre.s.0 <- Diagonal(x = rep((1/parm.mat[i, "tau2.s"]), nlocs))
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff  - offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff 
      }
      mu.s.post <- Sigma.s.post%*%
        (t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 1/rgamma(1, nlocs/2 + a1, 
                                          b1 + as.numeric(sum(theta.mat[i+1, ]^2))/2)
      
    }else if(spatial.str == "CAR"){
      
      pre.s.0 <- (1/parm.mat[i, "tau2.s"])*(Dw.s - parm.mat[i, "rho.s"]*W.s)
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff - offset.vec
      }else{
        res.vec.i <-  z.vec - fix.eff 
      }
      mu.s.post <- Sigma.s.post%*%(t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 
        1/rgamma(1, nlocs/2 + a1, 
                 b1 + as.numeric(t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                                   (Dw.s - parm.mat[i, "rho.s"]*W.s)%*%
                                   matrix(theta.mat[i+1, ], ncol = 1))/2)
      
      
      # 5. update rho.t (M-H)
      inter = as.numeric( t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                            W.s%*%matrix(theta.mat[i+1, ], ncol = 1) )
      ll.rho = ll.rho.s.1 + rho1.prior.val/(2*parm.mat[i+1, "tau2.s"]^2)*inter
      parm.mat[i+1, "rho.s"] <- sample(x = rho1.prior.val, size = 1, 
                                       prob = exp(ll.rho - max(ll.rho)))
      
    }
    
    
    # 9. update xi (over dispersion) 
    if(!is.null(offset.vec)){
      eta.vec = offset.vec + X.design.fix%*%matrix(beta.mat[i+1, ], ncol = 1) + 
        as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
      eta.vec <- as.numeric(eta.vec)
    }else{
      eta.vec = X.design.fix%*%matrix(beta.mat[i+1, ], ncol = 1) + 
        as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
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
                          theta.mat[-c(1:burn_in), ])
  X.design.comb <- cbind(X.design.fix, 
                         X.design.spatial)
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
                         X.design.spatial)
  X.design.null <- as.matrix(X.design.null)
  pred.val.null.mat <- get_predvals(parm_post = parm.post.comb,
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
               rand.int.s = theta.mat,
               pred.counts = pred.counts.nb.mat,
               pred.null.counts = pred.counts.null.nb.mat,
               AD = AD.mat, 
               accept.rate = mean(accept.mat[-c(1:burn_in), 1]), 
               WAIC = -2*(lppd - pWAIC), 
               lppd = lppd, pWAIC = pWAIC)
  }
  
  return(re)
  
}






