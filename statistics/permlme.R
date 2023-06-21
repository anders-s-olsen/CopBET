permlme <- function(lme0, lme1, data = NULL, seed = NULL, 
                    statistic = "Wald", perm.X = FALSE, return.index.perm = FALSE,
                    nperm = 1000, cpus = 1, trace = TRUE){ 
  
  statistic <- match.arg(statistic, c("LRT","Wald"), several.ok = TRUE)
  if(is.null(data)){data <- nlme::getData(lme1)}
  data <- as.data.frame(data)
  
  ## ** 0- normal test
  name.coef <- setdiff(names(fixef(lme1)),names(fixef(lme0)))
  if(length(name.coef)!=1 && "Wald" %in% statistic){
    stop("There should be exactly one parameter more in lme1 than in lme0. \n")
  }
  lme0.ML <- update(lme0, data = data, method = "ML")
  lme1.ML <- update(lme1, data = data, method = "ML")
  out <- list(call = match.call(),
              LRT = anova(lme0.ML, lme1.ML),
              Wald = summary(update(lme1, data = data, method = "REML"))$tTable[name.coef,,drop=FALSE],
              stat.perm = NULL,
              n.perm = c(LRT = NA, Wald = NA))
  
  ## ** 1- extract key quantities from input
  n.obs <- NROW(data) ## number of observations
  
  cluster <- getGroups(lme0) ## to which cluster (patient) each observation belongs
  U.cluster <- levels(cluster)
  n.cluster <- length(U.cluster)
  index.cluster <- lapply(U.cluster, function(iCluster){which(cluster==iCluster)}) ## used to restaure proper ordering after tapply
  
  Y <- getResponse(lme1) ## response
  name.Y <- all.vars(formula(lme1))[[1]] ## name of the response variable
  X0 <- model.matrix(formula(lme0),data) ## design matrix
  
  Omega0 <- getVarCov(lme0, type = "marginal", individuals = levels(cluster)) ## residual variance-covariance matrix
  beta0 <- fixef(lme0) ## regression coefficients
  
  ## ** 2- compute residuals
  Xbeta0 <- X0 %*% beta0
  residuals0 <- as.double(Y - X0 %*% beta0)
  
  ## ** 3- compute normalized residuals
  sqrtOmega0 <- lapply(Omega0,function(iOmega0){t(chol(iOmega0))})
  sqrtOmega0M1 <- lapply(sqrtOmega0,solve)
  
  residuals0N <- vector(length=n.obs, mode = "numeric")
  for(iCluster in 1:n.cluster){ ## iCluster <- 1
    residuals0N[index.cluster[[iCluster]]] <- sqrtOmega0M1[[iCluster]] %*% residuals0[index.cluster[[iCluster]]]
  }
  
  ## ** 4- estimate the distribution of the test statistics under the null
  warper <- function(iPerm){
    data.perm <- data
    
    ## permute residuals and fixed effects
    index.perm <- sample(1:n.obs)
    residuals0N.perm <- residuals0N[index.perm]
    
    ## rescale residuals
    for(iCluster in 1:n.cluster){ ## iCluster <- 1
      data.perm[[name.Y]][index.cluster[[iCluster]]] <- sqrtOmega0[[iCluster]] %*% residuals0N.perm[index.cluster[[iCluster]]]
    }
    ## add permuted fixed effects
    if(perm.X){
      data.perm[[name.Y]] <- data.perm[[name.Y]] + Xbeta0[index.perm,,drop=FALSE]
    }
    
    if("LRT" %in% statistic){
      lme0.permML <- try(update(lme0, data = data.perm, method = "ML"), silent = TRUE)
      lme1.permML <- try(update(lme1, data = data.perm, method = "ML"), silent = TRUE)
      if(inherits(lme0.permML,"try-error")||inherits(lme0.permML,"try-error")){
        LRT.stat <- NA
      }else{
        LRT.stat <- as.double(2*(logLik(lme1.permML)-logLik(lme0.permML)))
      }
    }else{
      LRT.stat <- NA
    }
    
    if("Wald" %in% statistic){
      lme1.permREML <- try(update(lme1, data = data.perm, method = "REML"), silent = TRUE)
      if(inherits(lme1.permREML,"try-error")){
        Wald.stat <- NA
      }else{
        Wald.stat <- summary(lme1.permREML)$tTable[name.coef,"t-value"]
      }
    }else{
      Wald.stat <- NA
    }
    iOut <- c(LRT = LRT.stat, Wald = Wald.stat)
    if(return.index.perm){
      attr(iOut,"index.perm") <- index.perm
    }
    return(iOut)
  }
  
  
  if(cpus==1 || is.null(cpus)){
    if(!is.null(seed)){set.seed(seed)}
    
    if(trace){
      requireNamespace("pbapply")
      ls.perm <- pbapply::pblapply(1:nperm, function(i){warper(i)})
    }else{
      ls.perm <- lapply(1:nperm, function(i){warper(i)})
    }
  }else{
    cl <- snow::makeSOCKcluster(cpus)
    doSNOW::registerDoSNOW(cl)
    
    pb <- txtProgressBar(max = nperm, style=3)
    opts <- list(progress = function(n) setTxtProgressBar(pb, n))
    if(is.null(seed)){
      ls.perm <- foreach::`%dopar%`(
        foreach::foreach(i=1:nperm, .options.snow=opts, .packages = "nlme"), {
          warper(i)
        })
    }else{
      requireNamespace(doRNG)
      set.seed(seed)
      ls.perm <- doRNG::`%dorng%`(
        foreach::foreach(i=1:nperm, .options.snow=opts, .packages = "nlme"), {
          warper(i)
        })
    }
    parallel::stopCluster(cl)
  }
  out$stat.perm <- do.call(rbind,ls.perm)
  if(return.index.perm){
    out$index.perm <- do.call(cbind,lapply(ls.perm,attr,"index.perm"))
  }
  
  ## ** 5- compare the observed statistic with its null distribution
  n.permLRT <- sum(!is.na(out$stat.perm[,"LRT"]))
  if(n.permLRT>0){
    out$LRT[["p.value.perm"]] <- c(NA,(sum(pmax(out$stat.perm[,"LRT"],0) >= out$LRT$L.Ratio[2], na.rm = TRUE) + 1)/(n.permLRT+1))
  }
  out$n.perm["LRT"]  <- n.permLRT
  
  n.permWald <- sum(!is.na(out$stat.perm[,"Wald"]))
  if(n.permWald>0){
    out$Wald <- data.frame(out$Wald, "p.value.perm" = (sum(abs(out$stat.perm[,"Wald"]) >= abs(out$Wald[,"t-value"]), na.rm = TRUE) + 1)/(n.permWald+1))
  }
  out$n.perm["Wald"]  <- n.permWald
  
  
  ## ** export
  class(out) <- append("permlme", class(out))
  return(out)
}

## * print.permlme
print.permlme <- function(x,...){
  if("p.value.perm" %in% names(x$LRT)){
    cat("Likelihood ratio test: \n")
    print(x$LRT)
  }
  if("p.value.perm" %in% names(x$Wald)){
    cat("Wald test: \n")
    print(x$Wald)
  }
  return(NULL)
}