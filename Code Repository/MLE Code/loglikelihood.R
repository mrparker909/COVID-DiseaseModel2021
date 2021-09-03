# Disease Analytics (COVID-19) Log Likelihood Function

#' @title loglik
#' @description   Used to calculate the negative of the log likelihood.
#' @param par     Vector with six elements (if no covariates): log(lambda), log(gamma), log(omega), logis(pdet), logis(pmort), logis(prec). If there are covariates, include a starting value for each covariate. Order is: gamma, gamma site covariates, gamma time covariates, omega, omega site covariates, omega time covariates, pdet, pdet site covariates, pdet time covariates, pmort, pmort site covariates, pmort time covariates, prec, prec site covariates, prec time covariates.
#' @param nit     R by T matrix of observed counts with R sites/rows and T sampling occassions/columns.
#' @param Dit     Deaths are assumed fully observed. R by T matrix of death counts with R sites/rows and T sampling occassions/columns.
#' @param rit     R by T matrix of observed recovery counts with R sites/rows and T sampling occassions/columns.
#' @param l_s_c   list of lambda (mean initial infections) site covariates (list of vectors of length R (number of sites), each vector is a covariate)
#' @param g_s_c   list of gamma (mean imported infections) site covariates (list of vectors of length R (number of sites), each vector is a covariate)
#' @param g_t_c   list of gamma (mean imported infections) time covariates (list of vectors of length T (number of sampling occasions), each vector is a covariate)
#' @param o_s_c   list of omega (mean infections per domestic infected) site covariates (list of vectors of length R (number of sites), each vector is a covariate)
#' @param o_t_c   list of omega (mean infections per domestic infected) time covariates (list of vectors of length T (number of sampling occasions), each vector is a covariate)
#' @param p_s_c   list of pdet (probability of detection) site covariates (list of vectors of length R (number of sites), each vector is a covariate)
#' @param p_t_c   list of pdet (probability of detection) time covariates (list of vectors of length T (number of sampling occasions), each vector is a covariate)
#' @param pm_s_c  list of pm (probability of death) site covariates (list of vectors of length R (number of sites), each vector is a covariate)
#' @param pm_t_c  list of pm (probability of death) time covariates (list of vectors of length T (number of sampling occasions), each vector is a covariate)
#' @param pr_s_c  list of pr (probability of recovery) site covariates (list of vectors of length R (number of sites), each vector is a covariate)
#' @param pr_t_c  list of pr (probability of recovery) time covariates (list of vectors of length T (number of sampling occasions), each vector is a covariate)
#' @param K       Upper bound on summations (population upper bound). Currently only a vector of K values (one entry for each site) is possible, eg: K=rep(100,times=R).
#' @param constLambda If not NULL, then lambda is not estimated (par[1] is ignored), instead constLambda contains the true value of lambda (must be an integer value between 0 and K). This is useful for setting the initial population size to zero with: constLambda = 0.
#' @param constGamma If not NULL, then gamma is not estimated (par[2] is ignored), instead constGamma contains the true value of gamma (must be a non-negative real number). This is useful for setting the population importation rate to zero with: constGamma = 0.
#' @param HPDfile File name to output the posterior density to. Posterior density is stored as an array (i,t,N-1), where each entry is the posterior likelihood of true population N at site i and time t.
#' @param outfile File name to output the current parameter values and negative log likelihood value.
#' @param VERBOSE If TRUE, prints the log likelihood to console.
#' @details Note that this function is adapted from the negative log likelihood function from the R package unmarked (Fiske and Chandler 2019), and uses the recursive method of computation described in Web Appendix A of Dail and Madsen 2011: Models for Estimating Abundance from Repeated Counts of an Open Metapopulation, published in Biometrics volume 67, issue 2.
#' @references Fiske, I., Chandler, R., Miller, D., Royle, A., Kery, M., Hostetler, J., Hutchinson, R., Smith, A., & Kellner, K. (2019). unmarked: Models for Data from Unmarked Animals (Version 0.13-1) [Computer software]. https://CRAN.R-project.org/package=unmarked
#' @export
loglik <- function(par, nit, K, Dit=NULL, rit=NULL,
                   l_s_c=NULL, g_s_c=NULL, g_t_c=NULL, 
                   o_s_c=NULL, o_t_c=NULL, p_s_c=NULL, 
                   p_t_c=NULL, pm_s_c=NULL, pm_t_c=NULL, 
                   pr_s_c=NULL, pr_t_c=NULL, constLambda=NULL,
                   constGamma=NULL,
                   HPDfile=NULL, outfile = NULL, VERBOSE=FALSE) {
  T <- ncol(nit)
  R <- nrow(nit)

  if(is.null(Dit)) { Dit = matrix(0, nrow=R, ncol=T) }
  if(is.null(rit)) { rit = matrix(0, nrow=R, ncol=T) }
  
  if(!is.null(HPDfile)) { HPDMat = array(0, dim = c(R,T,K[1,1]+1)) }
  
  # NOTE: pdet = probability of detection
  #       pmor = probability of mortality/death
  #       prec = probability of recovery
  
  param_length = 1
  
  # extract lambda estimates from par, setup lambda covariate matrix lamb, and covariate vector B_l
  lamb <- matrix(rep(1,times=R),ncol=1) # site covariates for lambda
  B_l <- par[param_length] # coefficients for lambda
  if(!is.null(l_s_c)) {
    lamb <- cbind(lamb, do.call(cbind, l_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
    B_l <- sapply(X = param_length - 1 + 1:(length(l_s_c)+1), FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_l)
  
  # extract gamma estimates from par, setup gamma covariate matrix gamm, and covariate vector B_g
  gamm <- matrix(rep(1,times=R),ncol=1) # site covariate coefficients for gamma
  B_g <- par[param_length] # covariates for gamma
  if(!is.null(g_s_c)) {
    gamm <- cbind(gamm, do.call(cbind, g_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
    B_g <- sapply(X = param_length - 1 + 1:(length(g_s_c)+1), FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_g)
  
  # extract gamma estimates from par, setup gamma covariate matrix gamt, and covariate vector B_gt
  gamt <- matrix(rep(0,times=T),ncol=1) # time covariate coefficients for gamma
  B_gt <- NULL # covariates for gamma
  if(!is.null(g_t_c)) {
    gamt <- do.call(cbind, g_t_c) # rows are times, cols are covariate values, here we are creating the design matrix
    B_gt <- sapply(X = param_length + 1:(length(g_t_c)) - 1, FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_gt)
  
  # extract omega estimates from par, setup omega covariate matrix omeg, and covariate vector B_o
  omeg <- matrix(rep(1,times=R),ncol=1) # site covariate coefficients for omega
  B_o <- par[param_length] # covariates for omega
  if(!is.null(o_s_c)) {
    omeg <- cbind(omeg, do.call(cbind, o_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
    B_o <- sapply(X = param_length - 1 + 1:(length(o_s_c)+1), FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_o)
  
  # extract omega estimates from par, setup omega covariate matrix omet, and covariate vector B_ot
  omet <- matrix(rep(0,times=T),ncol=1) # time covariate coefficients for omega
  B_ot <- NULL # covariates for omega
  if(!is.null(o_t_c)) {
    omet <- do.call(cbind, o_t_c) # rows are times, cols are covariate values, here we are creating the design matrix
    B_ot <- sapply(X = param_length + 1:(length(o_t_c)) - 1, FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_ot)
  
  # extract pdet estimates from par, setup pdet covariate matrix pdet, and covariate vector B_p
  pdet <- matrix(rep(1,times=R),ncol=1) # site covariate coefficients for pdet
  B_p <- par[param_length] # covariates for pdet
  if(!is.null(p_s_c)) {
    pdet <- cbind(pdet, do.call(cbind, p_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
    B_p <- sapply(X = param_length - 1 + 1:(length(p_s_c)+1), FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_p)
  
  # extract pdet estimates from par, setup pdet covariate matrix pt, and covariate vector B_pt
  pt   <- matrix(rep(0,times=T),ncol=1) # time covariate coefficients for pdet
  B_pt <- NULL # covariates for pdet
  if(!is.null(p_t_c)) {
    pt   <- do.call(cbind, p_t_c) # rows are times, cols are covariate values, here we are creating the design matrix
    B_pt <- sapply(X = param_length + 1:(length(p_t_c)) - 1, FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_pt)
  
  # extract pmor estimates from par, setup pmor covariate matrix pmor, and covariate vector B_pm
  pmor <- matrix(rep(1,times=R),ncol=1) # site covariate coefficients for pmor
  B_pm <- par[param_length] # covariates for pmor
  if(!is.null(pm_s_c)) {
    pmor <- cbind(pmor, do.call(cbind, pm_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
    B_pm <- sapply(X = param_length - 1 + 1:(length(pm_s_c)+1), FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_pm)
  
  # extract pmor estimates from par, setup pmor covariate matrix pmt, and covariate vector B_pmt
  pmt   <- matrix(rep(0,times=T),ncol=1) # time covariate coefficients for pdet
  B_pmt <- NULL # covariates for pmor
  if(!is.null(pm_t_c)) {
    pmt   <- do.call(cbind, pm_t_c) # rows are times, cols are covariate values, here we are creating the design matrix
    B_pmt <- sapply(X = param_length + 1:length(pm_t_c) - 1, FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_pmt)
  
  # extract prec estimates from par, setup prec covariate matrix prec, and covariate vector B_pr
  prec <- matrix(rep(1,times=R),ncol=1) # site covariate coefficients for prec
  B_pr <- par[param_length] # covariates for prec
  if(!is.null(pr_s_c)) {
    prec <- cbind(prec, do.call(cbind, pr_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
    B_pr <- sapply(X = param_length - 1 + 1:(length(pr_s_c)+1), FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_pr)
  
  # extract prec estimates from par
  prt   <- matrix(rep(0,times=T),ncol=1) # time covariate coefficients
  B_prt <- NULL # covariates
  if(!is.null(pr_t_c)) {
    prt   <- do.call(cbind, pr_t_c) # rows are times, cols are covariate values, here we are creating the design matrix
    B_prt <- sapply(X = param_length + 1:length(pr_t_c) - 1, FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_prt)
  
  Y         <- nit
  mode(Y)   <- "numeric"
  g1_t_star <- list()
  g1_t      <- list()
  g1        <- list()
  g2        <- list()
  g_star    <- list()
  g3        <- vector(length = R, mode = "list")
  g3_temp   <- list()

  
  if(is.null(g_t_c)) {
    t_gamt=matrix(0,nrow=T,ncol=1)
    t_B_gt=0
  } else {
    t_gamt=gamt
    t_B_gt=B_gt
  }
  if(is.null(o_t_c)) {
    t_omet=matrix(0,nrow=T,ncol=1)
    t_B_ot=0
  } else {
    t_omet=omet
    t_B_ot=B_ot
  }
  
  for(i in 1:R) {
    # check if tpMAT is time dependent
    time_dependent = TRUE
    if(is.null(g_t_c) && is.null(o_t_c) && is.null(p_t_c) && is.null(pr_t_c) && is.null(pm_t_c)) {
      time_dependent = FALSE
    }
    
    # tempMat holds transition probabilities in log space
    tempMat <- matrix(-Inf, nrow = K[i]+1, ncol = K[i]+1)
    for(t in 1:T) {
      if(t == 1 | time_dependent) {
        
        t_pmor <- plogis(sum(pmor[i,] * B_pm)+sum(pmt[t,] * B_pmt))
        t_prec <- plogis(sum(prec[i,] * B_pr)+sum(prt[t,] * B_prt))
        
        t_gamm <- gamm[i,]
        if(!is.null(constGamma)) { t_gamm = constGamma }
        
        g3_temp <- log_tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=t_omet[t,], B_ot=t_B_ot, gamm=t_gamm, B_g=B_g, gamt=t_gamt[t,], B_gt=t_B_gt)

        g3[[i]][[t]] <- g3_temp
      } else {
        g3[[i]][[t]] <- g3[[i]][[1]]
      }
    }
  }
  
  for(i in 1:R) {
    g1_t_star[[i]] <- rep(0, times=K[i]+1)
    g1_t[[i]]      <- numeric(K[i]+1)
    g1[[i]]        <- numeric(K[i]+1)
    g2[[i]]        <- numeric(K[i]+1)
    g_star[[i]]    <- matrix(0, nrow=K[i]+1, ncol=1)
  }
  
  # apply over sites (1 to R), TODO: this is a prime candidate for parallel computing since each site i is independent
  ll_i  <- vapply(X = 1:R, FUN = function(i, K, T, Y, Dit, rit, pdet, B_p, pt, B_pt, g3, g1_t_star, g1_t,g1,g2, g_star) {
    g3        <- g3[[i]]
    g1_t_star <- g1_t_star[[i]]
    g1_t      <- g1_t[[i]]
    g1        <- g1[[i]]
    g2        <- g2[[i]]
    g_star    <- g_star[[i]]
    K         <- K[i]
    
    lli = 0
    for(t in (T-1):1) {
      t_pdet <- plogis(sum(pdet[i,] * B_p)+sum(pt[t,] * B_pt))
      t_pmor <- plogis(sum(pmor[i,] * B_pm)+sum(pmt[t,] * B_pmt))
      t_prec <- plogis(sum(prec[i,] * B_pr)+sum(prt[t,] * B_prt))
      
      if(t_pmor + t_prec > 1) { return(-Inf) }
      
      start = rit[i,t+1]
      Mult  <- matrix(-Inf, nrow=1, ncol=(K+1))
     
      # for each possible active population size A from 0 up to K
      for(A in 0:K) {
        end = A-Dit[i,t+1]
        # sum over all possible values of Rit from rit up to A
        if(start > end) { Mult[A+1] = -Inf } else {
          for(Rit in start:end) {
            # Note: Nit = A_{it+1} + D_{it+1} + R_{it+1}
            Nit = A+Dit[i,t+1]+Rit
            Mult2 = dmultinomial(Nit,Dit[i,t+1],Rit,t_pmor,t_prec,log=T)
            Mult[A+1] = logSumExp(c(Mult[A+1], Mult2))
          }
        }
      }

      lli = lli + Mult
    }
    
    # Note: size 0:K = N[t] - (observed active cases) = N[t]-a[t-1]+r[t-1]+D[t-1] = N[t] - (a[t] - n[t])
    a = NULL
    for(t in 1:T) {
      a = c(a, sum(nit[i,1:t])-sum(rit[i,1:t])-sum(Dit[i,1:t]))
    }
    
    # loop backwards over times t, stopping at t==2
    for(t in T:2) {
      t_pdet <- plogis(sum(pdet[i,] * B_p)+sum(pt[t,] * B_pt))
      t_pmor <- plogis(sum(pmor[i,] * B_pm)+sum(pmt[t,] * B_pmt))
      t_prec <- plogis(sum(prec[i,] * B_pr)+sum(prt[t,] * B_prt))
      
      RecovProb = dmultinomial(a[t-1], Dit[i,t], rit[i,t], t_pmor, t_prec, log = T)
      
      # size takes possible value of N (0 to K) at time t (for t > 1)
      g1_t <- dbinom(x=Y[i,t], size = pmax(0,(0:K)-(a[t]-nit[i,t])), prob = t_pdet, log = T)
      
      g1_t_star = g1_t + t(g_star) + RecovProb

      if(!is.null(HPDfile)) {
        HPDMat[i,t,] <<- g1_t_star
      }
      
      # update g_star for t-1
      #g_star = g3[[t]] %*% g1_t_star
      g_star = Ax_log(g3[[t]], g1_t_star)
    }

    # t == 1:
    t_pdet <- plogis(sum(pdet[i,] * B_p)+sum(pt[1,] * B_pt))
    
    # Here we can set N==constLambda at t==1
    g1 <- log(numeric(K+1))
    if(!is.null(constLambda)) {
      g1[constLambda+1] <- dbinom(x = nit[i,1], size = constLambda, prob = t_pdet, log = T)
    } else {
      # apply binomial thinning at t==1
      g1 <- dpois(x = 0:K, lambda = exp(sum(lamb[i,] * B_l)), log = T) +
       dbinom(x = nit[i,1], size = 0:K, prob = t_pdet, log = T) # Note: size = N[1]-a[0]+r[0]+D[0] = N[1]
      
      if(nit[i,1] > 0) {
        g1[1:(nit[i,1])] <- -Inf
      }
    }

    g1 = g1 + t(lli)

    if(!is.null(HPDfile)) {
      HPDMat[i,1,] <<- g1
    }

    # apply recursive definition of likelihood
    return( logSumExp(g1 + g_star) )
  }, FUN.VALUE = numeric(1), K=K, T=T, Y=Y, Dit=Dit, rit=rit, pdet=pdet, B_p=B_p, pt=pt, B_pt=B_pt, g3=g3, g1_t_star=g1_t_star, g1_t=g1_t,g1=g1,g2=g2,g_star=g_star)
  
  ll <- sum(unlist(ll_i))
  if(is.nan(ll)) { ll <- -Inf }
  
  if(!is.null(HPDfile)) { 
    saveRDS(HPDMat, file = HPDfile)
  }
  
  if(!is.null(outfile)) { 
    datcsv = list(par=par, nll=-1*ll)
    saveRDS(datcsv, file = outfile)
  }
  
  if(VERBOSE) {
    if(typeof(ll)=="double") {
      try(print(paste0("log likelihood: ", ll)))
      try(print(paste0("parameters: ", par)))
    }
  }

  return(-1*ll)
}

log_tp_MAT <- function(M, omeg, B_o, omet, B_ot, gamm, B_g, gamt, B_gt) {
  if(is.null(B_o )) { B_o  = 1 }
  if(is.null(B_ot)) { B_ot = 0 }
  if(is.null(B_g )) { B_g  = 1 }
  if(is.null(B_gt)) { B_gt = 0 }
  
  omega = exp(sum(omeg * B_o)+sum(B_ot*omet))
  if(is.null(omega)) {stop("omega is NULL in tp_MAT")}
  gamma = exp(sum(B_g*gamm)+sum(B_gt*gamt))
  if(is.null(gamma)) {stop("gamma is NULL in tp_MAT")}
  
  # col-1 of M is cumulative total cases (at t)
  # row-1 of M is cumulative total cases (at t-1)
  # m = col - row is new cases at t
  K = nrow(M)
  
  # Note that we can't pass Rcpp functions through %dopar%
  M = foreach(row = 1:K, .combine="rbind", .export = c("conv_log_FFT", "conv_log", "logSumExp")) %dopar% {
    Mrow = M[row,]
    if(row-1==0) { # all new cases are imported cases
      Mrow = dpois(0:(K-1), gamma, log = TRUE)
    } else {   # cases may be imported or domestic
      romeg = max(0, omega*(row-1))
      # slow convolve in log space
      conv = conv_log(dpois(0:(K-1),gamma,log=TRUE), dpois(0:(K-1),romeg,log=TRUE))
      Mrow = conv[1:K]
    }
    return(Mrow)
  }
  return(M)
}

# multinomial density function
dmultinomial <- function(A, B, C, pB, pC, log = FALSE) 
{
  if(1-pB-pC <= 0) {
    ifelse(log==T, return(-Inf), return(0))
  }
  # note that there are three groups, A-B-C, B, C, with probabilities 1-pB-pC, pB, pC
  MC = lgamma(A+1) - lgamma(B+1) - lgamma(C+1) - lgamma(A-B-C+1)
  dmult = MC + B*log(pB) + C*log(pC) + (A-B-C)*log(1-pB-pC);
  if(!log) { dmult = exp(dmult) }
  
  return(dmult);
}

# compute convolution in log space
conv_log <- function(log_x,log_y) {
  
  nx = length(log_x)
  ny = length(log_y)
  if(nx != ny) {stop("length x must equal length y")}
  log_y = rev(log_y)
  n = nx
  m = 2*n - 1
  ans = numeric(m)
  
  # first half
  for(k in 1:n) {
    v = numeric(k)
    for(j in 1:k) {
      v[j] = log_x[j] + log_y[n-(k-j)]
    }

    if(max(v) == -Inf) { 
      ans[k] = -Inf 
    } else {
      ans[k] = logSumExp(v)
    }
  }
  
  # second half
  for(k in 1:(n-1)) {
    v = numeric(k)
    for(j in 1:k) {
      v[j] = log_x[n-(k-j)] + log_y[j]
    }

    if(max(v) == -Inf) { 
      ans[m-k+1] = -Inf 
    } else {
      ans[m-k+1] = logSumExp(v)
    }
  }
  
  return(ans)
}

logSumExp <- function(x) {
  if(all(is.infinite(x))) { return(x[1]) }
  x = x[which(is.finite(x))]
  ans = x[1]
  for(i in seq_along(x)[-1]) {
    ma = max(ans,x[i])
    mi = min(ans,x[i])
    ans = ma + log1p(exp(mi-ma))
  }
  return(ans)
}


# Multiply a matrix times a vector y=Ax in log space
Ax_log <- function(logA,logx) {
  xl = length(logx)
  if(is.null(xl)) { xl = nrow(logx) }
  Ac = ncol(logA)
  if(Ac != xl) {stop("columns of A must match length of x")}
  
  y = matrix(-Inf, nrow=nrow(logA), ncol=1)
  for(r in 1:nrow(logA)) {
    temp = logA[r,] + logx
    y[r,1] = logSumExp(temp)
  }
  return(y)
}
