### function needed
library(Matrix)
library(ks)
library(pbapply)
library(parallel)
library(MASS)
library(tictoc)
library(Surrogate)
#### function to compute the ICA using the PC algorithm
ICA.ContCont.MultS.PC = function (M = 1000, N, Sigma, Seed = 123, Show.Progress = FALSE) 
{
  is.PD = function(X, tol = 1e-08) {
    X[is.na(X)] = 0
    min.lambda = min(eigen(X, only.values = T, symmetric = T)$values)
    return(min.lambda > tol)
  }
  r.random = function(j, k, R) {
    d = ncol(R)
    r1 = R[j, (j + 1):(j + k - 1)]
    r3 = R[(j + 1):(j + k - 1), j + k]
    R2.inv = R[(j + 1):(j + k - 1), (j + 1):(j + k - 1)]
    R2 = solve(R2.inv)
    #r.c = extraDistr::rnsbeta(1, parm.b, parm.b, -1, 1)
    r.c = extraDistr::rnsbeta(1, 1 + 0.5 * (d - 1 - k), 1 + 
                                0.5 * (d - 1 - k), -1, 1)
    D2 = (1 - tcrossprod(crossprod(r1, R2), r1)) * (1 - tcrossprod(crossprod(r3, 
                                                                             R2), r3))
    if (D2 < 0 & D2 > -1e-08) {
      D2 = 0
    }
    r = tcrossprod(crossprod(r1, R2), r3) + r.c * sqrt(D2)
    return(r)
  }
  Correlation.matrix.PC = function(R, Range = c(-1, 1)) {
    Rf = R
    Rf[is.na(Rf)] = 0
    d = ncol(R)
    j.ind = do.call("c", lapply(1:(d - 1), function(x) {
      x:1
    }))
    k.ind = do.call("c", lapply(1:(d - 1), function(x) {
      1:x
    }))
    for (i in 1:(length(j.ind))) {
      j = j.ind[i]
      k = k.ind[i]
      if (is.na(R[j, j + k])) {
        if (k == 1) {
          R[j, j + k] = R[j + k, j] = extraDistr::rnsbeta(1, d/2, d/2, -1, 1)
        }
        else {
          R[j, j + k] = R[j + k, j] = r.random(j, k, 
                                               R)
        }
        if (is.nan(R[j, j + k])) {
          stop("error")
        }
      }
    }
    return(R)
  }
  MutivarICA.fun = function(R, Sigma, N) {
    d = nrow(R)
    sdMat = diag(sqrt(Sigma))
    rtn = sdMat %*% R %*% t(sdMat)
    var_diff <- function(cov_mat) {
      cov_val <- cov_mat[1, 1] + cov_mat[2, 2] - (2 * cov_mat[1, 
                                                              2])
      fit <- c(cov_val)
      fit
    }
    cov_2_diffs <- function(cov_mat) {
      cov_val <- (cov_mat[2, 2] - cov_mat[1, 2]) - (cov_mat[2, 
                                                            1] - cov_mat[1, 1])
      fit <- c(cov_val)
      fit
    }
    A <- matrix(var_diff(cov_mat = rtn[1:2, 1:2]), nrow = 1)
    B <- NULL
    Aantal <- (dim(R)[1] - 2)/2
    rtn_part <- rtn[c(3:dim(rtn)[1]), c(1, 2)]
    for (z in 1:Aantal) {
      cov_mat_here <- rtn_part[c((z * 2) - 1, z * 2), c(1:2)]
      B <- rbind(B, cov_2_diffs(cov_mat_here))
    }
    Dim <- dim(R)[1]
    Sub_mat_var <- rtn[c(3:Dim), c(3:Dim)]
    C <- matrix(NA, nrow = Aantal, ncol = Aantal)
    for (l in 1:Aantal) {
      for (k in 1:Aantal) {
        Sub_mat_var_hier <- Sub_mat_var[c((k * 2) - 1, 
                                          k * 2), c((l * 2) - 1, l * 2)]
        C[k, l] <- cov_2_diffs(cov_mat = Sub_mat_var_hier)
      }
    }
    Delta <- cbind(rbind(A, B), rbind(t(B), C))
    ICA <- (t(B) %*% solve(C) %*% B)/A
    Adj.ICA <- 1 - (1 - ICA) * ((N - 1)/(N - Aantal - 1))
    return(c(ICA, Adj.ICA))
  }
  set.seed(Seed)
  d = nrow(Sigma)[1]
  Vars = diag(Sigma)
  IND = ks::vec(matrix(1:d, ncol = 2, byrow = T), byrow = F)
  Sigma = Sigma[IND, IND]
  R = cov2cor(Sigma)
  if (!is.PD(R)) {
    alpha = uniroot(function(alpha, R, Rfixed, tol = 1e-04) {
      if (anyNA(R)) {
        R[is.na(R)] = 0
      }
      f = alpha * R + (1 - alpha) * Rfixed
      min(eigen(f)$values) - tol
    }, c(0, 1), R = R, Rfixed = diag(d), tol = 1e-08)$root
    R = R * alpha + (1 - alpha) * diag(d)
    warning(paste("The initial correlation matrix is not PD. TThe matrix was shrunk by a factor alpha=", 
                  alpha, " for correction", sep = ""))
  }
  Results = pbmapply(function(x) {
    IND.2 =  c(sample(1:(d/2)),sample((d/2+1):d))
    R.test = R[IND.2,IND.2]
    R.random = tryCatch(Correlation.matrix.PC(R.test,Range = c(-1, 1)),error=function(e){NULL})
    if(is.null(R.random)){
      R.random = matrix(NA,d,d)
      return(c(NA,NA, R.random[lower.tri(R.random)]))
    }
    R.random = R.random[order(IND.2),order(IND.2)]
    IND = ks::vec(matrix(1:d, ncol = 2), byrow = TRUE)
    R.random = R.random[IND, IND]
    ICA = MutivarICA.fun(R.random, Vars, N)
    if (Show.Progress == TRUE) {
      cat((x/M) * 100, "% done... ", sep = "")
    }
    return(c(ICA, R.random[lower.tri(R.random)]))
  }, x = 1:M)
  R2_H = Results[1, ]
  Corr.R2_H = Results[2, ]
  Outcome = list(R2_H = R2_H, Corr.R2_H = Corr.R2_H, Lower.Dig.Corrs.All = t(Results[-c(1:2), 
  ]))
  class(Outcome) <- "ICA.ContCont.MultS"
  return(Outcome)
}
#### function to compute the ICA using the modified PC algorithm (using three cores in parallel)
ICA.ContCont.MultS.MPC = function (M = 1000, N, Sigma,prob=NULL, Seed = 123, Save.Corr=F, nCores=3,
                                   Show.Progress = FALSE) 
{
  # M: number of simulations
  # Sigma: incomplete matrix of variance
  # prob: probabilities for each combination of r surrogates
  # Seed: seed
  # save.corr: save the generated correlation matrices
  # Show.Progress: show a progress bar
  is.PD = function(X, tol = 1e-08) {
    X[is.na(X)] = 0
    min.lambda = min(eigen(X, only.values = T, symmetric = T)$values)
    return(min.lambda > tol)
  }
  r.random = function(j, k, R) {
    d = ncol(R)
    r1 = R[j, (j + 1):(j + k - 1)]
    r3 = R[(j + 1):(j + k - 1), j + k]
    R2.inv = R[(j + 1):(j + k - 1), (j + 1):(j + k - 1)]
    R2 = solve(R2.inv)
    #r.c = extraDistr::rnsbeta(1, parm.b, parm.b, -1, 1)
    r.c = extraDistr::rnsbeta(1, 1 + 0.5 * (d - 1 - k), 1 + 
                                0.5 * (d - 1 - k), -1, 1)
    D2 = (1 - tcrossprod(crossprod(r1, R2), r1)) * (1 - tcrossprod(crossprod(r3, 
                                                                             R2), r3))
    if (D2 < 0 & D2 > -1e-08) {
      D2 = 0
    }
    r = tcrossprod(crossprod(r1, R2), r3) + r.c * sqrt(D2)
    return(r)
  }
  Correlation.matrix.PC = function(R, Range = c(-1, 1),parm.a) {
    Rf = R
    Rf[is.na(Rf)] = 0
    d = ncol(R)
    j.ind = do.call("c", lapply(1:(d - 1), function(x) {
      x:1
    }))
    k.ind = do.call("c", lapply(1:(d - 1), function(x) {
      1:x
    }))
    for (i in 1:(length(j.ind))) {
      j = j.ind[i]
      k = k.ind[i]
      if (is.na(R[j, j + k])) {
        if (k == 1) {
          R[j, j + k] = R[j + k, j] = extraDistr::rnsbeta(1, parm.a, parm.a, -1, 1)
        }
        else {
          R[j, j + k] = R[j + k, j] = r.random(j, k, R)
        }
        if (is.nan(R[j, j + k])) {
          stop("error")
        }
      }
    }
    return(R)
  }
  MultivarICA.fun = function(R, Sigma, N) {
    d = nrow(R)
    sdMat = diag(sqrt(Sigma))
    rtn = sdMat %*% R %*% t(sdMat)
    var_diff <- function(cov_mat) {
      cov_val <- cov_mat[1, 1] + cov_mat[2, 2] - (2 * cov_mat[1, 
                                                              2])
      fit <- c(cov_val)
      fit
    }
    cov_2_diffs <- function(cov_mat) {
      cov_val <- (cov_mat[2, 2] - cov_mat[1, 2]) - (cov_mat[2, 
                                                            1] - cov_mat[1, 1])
      fit <- c(cov_val)
      fit
    }
    A <- matrix(var_diff(cov_mat = rtn[1:2, 1:2]), nrow = 1)
    B <- NULL
    Aantal <- (dim(R)[1] - 2)/2
    rtn_part <- rtn[c(3:dim(rtn)[1]), c(1, 2)]
    for (z in 1:Aantal) {
      cov_mat_here <- rtn_part[c((z * 2) - 1, z * 2), c(1:2)]
      B <- rbind(B, cov_2_diffs(cov_mat_here))
    }
    Dim <- dim(R)[1]
    Sub_mat_var <- rtn[c(3:Dim), c(3:Dim)]
    C <- matrix(NA, nrow = Aantal, ncol = Aantal)
    for (l in 1:Aantal) {
      for (k in 1:Aantal) {
        Sub_mat_var_hier <- Sub_mat_var[c((k * 2) - 1, 
                                          k * 2), c((l * 2) - 1, l * 2)]
        C[k, l] <- cov_2_diffs(cov_mat = Sub_mat_var_hier)
      }
    }
    Delta <- cbind(rbind(A, B), rbind(t(B), C))
    ICA <- (t(B) %*% solve(C) %*% B)/A
    Adj.ICA <- 1 - (1 - ICA) * ((N - 1)/(N - Aantal - 1))
    return(c(ICA, Adj.ICA))
  }
  d = nrow(Sigma)[1]
  Vars = diag(Sigma)
  IND = ks::vec(matrix(1:d, ncol = 2, byrow = T), byrow = F)
  Sigma = Sigma[IND, IND]
  R = cov2cor(Sigma)
  if (!is.PD(R)) {
    alpha = uniroot(function(alpha, R, Rfixed, tol = 1e-04) {
      if (anyNA(R)) {
        R[is.na(R)] = 0
      }
      f = alpha * R + (1 - alpha) * Rfixed
      min(eigen(f)$values) - tol
    }, c(0, 1), R = R, Rfixed = diag(d), tol = 1e-08)$root
    R = R * alpha + (1 - alpha) * diag(d)
    warning(paste("The initial correlation matrix is not PD. TThe matrix was shrunk by a factor alpha=", 
                  alpha, " for correction", sep = ""))
  }
  p = (d-2)/2
  
  if(is.null(prob)){
    prob = choose(p,1:p)/sum(choose(p,1:p))
  }
  
  cl <- makeCluster(getOption("cl.cores", nCores))
  clusterExport(cl=cl, varlist=c('is.PD','r.random','Correlation.matrix.PC','MultivarICA.fun','p','R','Vars','prob',
                                 'N','M','d','Save.Corr'), envir=environment())
  
  clusterSetRNGStream(cl,Seed)
  
  if(Show.Progress==F){
    opb <- pboptions(type="none")
    on.exit(pboptions(opb))
  }
  
  Results = pblapply(X=1:M,function(X) {
    colnames(R) = rownames(R)=c('T0',paste('S',1:p,0,sep=''),'T1',paste('S',1:p,1,sep=''))
    sum.r = 0
    while(sum.r==0){
      r.num = sample(1:p,1,prob=prob)
      r = sample(c(rep(1,r.num),rep(0,p-r.num))) 
      sum.r=sum(r)
    }
    
    if(sum(r) < p){
      var.p = (1:p)[r == 1] 
      fixZero.p = (1:p)[r == 0]
      colr = 2+p + (1:p)[r==0]
      rowr = (1:p)[r==0] + 1
      R[1:(p+1),colr] = R[colr,1:(p+1)] = 0
      R[rowr,(2+p):d] = R[(2+p):d,rowr] = 0
      
      if(length(var.p)==1){
        if(length(fixZero.p)==1){
          p.index = c(var.p,fixZero.p)
        }else{
          p.index = c(var.p,sample(fixZero.p))        
        }
      }else{
        if(length(fixZero.p)==1){
          p.index = c(sample(var.p),fixZero.p)
        }else{
          p.index = c(sample(var.p),sample(fixZero.p))    
        }
      }
    }else{
      p.index=1:p
    }
    
    IND.2 =  c(1,1+p.index,p+2+p.index[p:1],p+2)
    R.test = R[IND.2,IND.2]
    parm.a = sum(r)+1
    R.random = tryCatch(Correlation.matrix.PC(R.test,Range = c(-1, 1),parm.a),error=function(e){NULL})
    if(is.null(R.random)){
      R.random = matrix(NA,d,d)
      return(c(NA,NA, r,R.random[lower.tri(R.random)]))
    }
    R.random = R.random[order(IND.2),order(IND.2)]
    IND = ks::vec(matrix(1:d, ncol = 2), byrow = TRUE)
    R.random = R.random[IND,IND]
    ICA = MultivarICA.fun(R.random, Vars, N)
    if(Save.Corr){
      return(c(ICA,r, R.random[lower.tri(R.random)]))      
    }else{
      return(c(ICA,r))      
    }
    
  }, cl=cl)
  stopCluster(cl)
  Results = do.call('rbind',Results)
  R2_H = Results[,1]
  Corr.R2_H = Results[,2 ]
  r = Results[,3:(p+2)]
  if(Save.Corr){
    Lower.Dig.Corrs.All = Results[,-c(1:(p+2))]
    Outcome = list(R2_H = R2_H, Corr.R2_H = Corr.R2_H, Lower.Dig.Corrs.All = Lower.Dig.Corrs.All,surr.eval.r=r)    
  }else{
    Outcome = list(R2_H = R2_H, Corr.R2_H = Corr.R2_H,surr.eval.r=r)    
  }
  return(Outcome)
}

ICA.ContCont.MultS = function (M = 500, N, Sigma, G = seq(from = -1, to = 1, by = 1e-05), 
                               Seed = c(123), Show.Progress = FALSE) 
{
  SDs <- sqrt(diag(Sigma))
  mu <- rep(0, times = length(SDs))
  results_here <- results <- all_delta_S_T <- all_delta_S_T_here <- Lower.Dig.Corrs.All <- NULL
  found <- 0
  set.seed(Seed)
  for (i in 1:M) {
    Sigma_c <- Sigma_c <- cov2cor(Sigma)
    num_elem <- dim(Sigma_c)[1]^2
    Sigma_c_orig <- Sigma_c
    size <- row_num <- col_num_ind <- 3
    total_size <- dim(Sigma_c)[1]
    here <- Sigma_c
    while (size <= total_size) {
      here <- Sigma_c[(1:size), (1:size)]
      here[is.na(here)] <- sample(x = G, size = length(here[is.na(here)]), 
                                  replace = TRUE)
      here[upper.tri(here)] = t(here)[upper.tri(here)]
      while (det(here) < 0) {
        here <- Sigma_c[(1:size), (1:size)]
        here[is.na(here)] <- sample(x = G, size = length(here[is.na(here)]), 
                                    replace = TRUE)
        here[upper.tri(here)] = t(here)[upper.tri(here)]
      }
      Sigma_c[1:row_num, 1:col_num_ind] <- here
      row_num <- row_num + 1
      col_num_ind <- col_num_ind + 1
      size <- size + 1
    }
    Sigma_c <- here
    Min.Eigen.Sigma <- try(min(eigen(Sigma_c)$values), TRUE)
    if ((class(Min.Eigen.Sigma) != "try-error") & (Min.Eigen.Sigma >= 
                                                   1e-11)) {
      found <- found + 1
      if (Show.Progress == TRUE) {
        prog = (i/M)*100
        if(as.integer(prog) == prog){
          cat((i/M) * 100, "% done... ", sep = "")          
        }
      }
      corMat <- Sigma_c
      varVec <- SDs^2
      n = nrow(corMat)
      sdMat = diag(sqrt(varVec))
      rtn = sdMat %*% corMat %*% t(sdMat)
      var_diff <- function(cov_mat) {
        cov_val <- cov_mat[1, 1] + cov_mat[2, 2] - (2 * 
                                                      cov_mat[1, 2])
        fit <- c(cov_val)
        fit
      }
      cov_2_diffs <- function(cov_mat) {
        cov_val <- (cov_mat[2, 2] - cov_mat[1, 2]) - 
          (cov_mat[2, 1] - cov_mat[1, 1])
        fit <- c(cov_val)
        fit
      }
      A <- matrix(var_diff(cov_mat = rtn[1:2, 1:2]), nrow = 1)
      B <- NULL
      Aantal <- (dim(Sigma_c)[1] - 2)/2
      rtn_part <- rtn[c(3:dim(rtn)[1]), c(1, 2)]
      for (z in 1:Aantal) {
        cov_mat_here <- rtn_part[c((z * 2) - 1, z * 2), 
                                 c(1:2)]
        B <- rbind(B, cov_2_diffs(cov_mat_here))
      }
      Dim <- dim(Sigma_c)[1]
      Sub_mat_var <- rtn[c(3:Dim), c(3:Dim)]
      C <- matrix(NA, nrow = Aantal, ncol = Aantal)
      for (l in 1:Aantal) {
        for (k in 1:Aantal) {
          Sub_mat_var_hier <- Sub_mat_var[c((k * 2) - 
                                              1, k * 2), c((l * 2) - 1, l * 2)]
          C[k, l] <- cov_2_diffs(cov_mat = Sub_mat_var_hier)
        }
      }
      Delta <- cbind(rbind(A, B), rbind(t(B), C))
      ICA <- (t(B) %*% solve(C) %*% B)/A
      Adj.ICA <- 1 - (1 - ICA) * ((N - 1)/(N - Aantal - 
                                             1))
      Lower.Dig.Corrs.Here <- Sigma_c[lower.tri(Sigma_c)]
      results_here <- (c(ICA, Adj.ICA))
      results <- rbind(results, results_here)
      Lower.Dig.Corrs.All <- data.frame(rbind(Lower.Dig.Corrs.All, 
                                              Lower.Dig.Corrs.Here))
    }
    Sigma_c <- Sigma_c_orig
  }
  row.names(results) <- NULL
  row.names(Lower.Dig.Corrs.All) <- NULL
  results <- data.frame(results, row.names = NULL)
  names(results) <- c("ICA", "Adj.ICA")
  fit <- list(R2_H = (as.numeric(results$ICA)), Corr.R2_H = (as.numeric(results$Adj.ICA)), 
              Lower.Dig.Corrs.Sigma = Lower.Dig.Corrs.All, Call = match.call())
  class(fit) <- "ICA.ContCont.MultS"
  fit
}