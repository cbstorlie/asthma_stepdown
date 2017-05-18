
try(source("/data5/bsi/clip/s201827.clip_hpc/ASD/source_code/weibull_probit_DPM_4.R"), silent=TRUE)
try(source("~/mayo/s112216.BedsidePatientRescue/source_code/weibull_probit_DPM_4.R"), silent=TRUE)

library("sas7bdat")
library("lme4")
library("rjags")
library("truncnorm")


## Treat y1=ED/IP on a fine increment dt (e.g., day or week)
## Treat y2=Rescue Inhaler and y3=OCS and on a by round resolution
##
## Med level and delta are time varying on dt scale, account for it in "integral"
##
## First pass, move every medication level increase to the last ED/IP visit (if there is one), otherwise to last OB/OP visit, otherwise leave at beginning of round.
##

## Need to input y1, y2, y3, N, M, U, J, L=3, X, unit, robs, dobs, W, dfT, step, visit, aa.rho, bb.rho, aa.pp, bb.pp, cuts
## N - # obs on a fine inc scale, i.e., length(y1)
## M - # obs on a by round scale, i.e., length(y2) or length(y3)
## J - # fixed effect covariates
## unit - "individual" label of the n observations
## robs - which obs on fine inc scale correspond to obs at round scale
## dobs - which obs on fine inc scale correspond to "individual" u
## step - length M; the final_step variable on 1-6 scale, not 0-5 (this is the "rounded, average" step level for the round)
## visit - length N giving 1,2,3 for no visit, OB/OP visit, or ED/IP visit, respectively
## aa.rho, bb.rho - vectors of length 3, the dbeta params for prob of med change at visit type k
## aa.Z - vector of length 6, the ddirichlet param for prob of med level = m

model.str <- "model{

  for(i in 1:N){
    y1[i] ~ dpois(lambda[i,1])
  }
  for(i in 1:M){
    y2[i] ~ dpois(rlambda[i,2])
    y3[i] ~ dpois(rlambda[i,3])
  }

  for(i in 1:N){
    for(l in 1:L){
      log(lambda[i,l]) <- A[unit[i],l] - B[l]*Z[i] + X[i,]%*%C[,l] + delta[i,l]
    }
  }
  for(i in 1:M){
    for(l in 1:L){
      rlambda[i,l] <- sum(lambda[robs[i,1]:robs[i,2],l])
    }
  }

 ## Model for A and B
  for(u in 1:U){
    for(l in 1:L){
      A[u,l] <- A.raw[u,l]
    }
    A.raw[u,1:L] ~ dmnorm(mu.a.raw, Tau.a.raw[,])
  }
  for(l in 1:L){
    mu.a.raw[l] <- mu.a[l]
    mu.a[l] ~ dnorm(0,.0001)
    B.raw[l] ~ dlnorm(-1,.1)
    B[l] <- B.raw[l] + offset[l]
  }
  Tau.a.raw ~ dwish(WA[,], dfT)

  for(j in 1:J){
    for(l in 1:L){
      C[j,l] ~ dnorm(0,.01)
    }
  }

  for(u in 1:U){
    mud[dobs[u,1],1:L] <- rep(0,L)
    for(i in (dobs[u,1]+1):dobs[u,2]){
      mud[i,1:L] <- phi*delta[i-1,1:L]
    }
    for(i in dobs[u,1]:dobs[u,2]){
      delta[i,1:L] ~ dmnorm(mud[i,], Tau.d[,])
    }
    Del.now[u,1:L] <- delta[dobs[u,2],1:L]
  }

  phi ~ dbeta(aa.phi,bb.phi)
  Tau.d ~ dwish(WD[,], dfT)

 ## distribution for Zc
  for(i in 1:N){
    for(k in 1:K){
      Zc.mat[i,k] ~ dnorm(muZ[i,k], precZ[i])
    }
    Zc[i] <- Zc.mat[i,BS[i]]
    Z[i] <- trunc(Zc[i]+.5)*(Zc[i]>=0.5)*(Zc[i]<4.5) + 5*(Zc[i]>=4.5)
    Zcg0[i] <- Zc[i]*(Zc[i]>0)
  }
  for(k in 1:K){
    p0[k] <- 1/K
  }
  for(u in 1:U){
    for(k in 1:K){
      muZ[dobs[u,1],k] <- muZ0
    }
    precZ[dobs[u,1]] <- precZ0
    BS[dobs[u,1]] ~ dcat(p0)
    for(i in (dobs[u,1]+1):dobs[u,2]){
      for(k in 1:K){
        p.BS[i,k] <- rho[visit[i-1]]/(K-1)*(BS[i-1]!=k) + (1-rho[visit[i-1]])*(BS[i-1]==k)
      }
      BS[i] ~ dcat(p.BS[i,1:K])
      for(k in 1:K){
        muZ[i,k] <- Zc.mat[i-1,k]
      }
      precZ[i] <- 200
    }
  }
  for(i in 1:M){
    meanZ[i] <- mean(Zcg0[robs[i,1]:robs[i,2]])
    step[i] ~ dinterval(meanZ[i],cuts)
  }

  muZ0 ~ dnorm(M.Z, P.Z)
  precZ0 <- .02

  for(k in 1:3){
    rho[k] ~ dbeta(aa.rho[k], bb.rho[k])
  }
  
 
 ## Summarize Tau.d as variances and correlations
  Sig.d <- inverse(Tau.d)
  for(l in 1:L){
    sig.d[l] <- Sig.d[l,l]
  }
  cor.d[1] <- Sig.d[1,2]/sqrt(sig.d[1]*sig.d[2])
  cor.d[2] <- Sig.d[1,3]/sqrt(sig.d[1]*sig.d[3])
  cor.d[3] <- Sig.d[2,3]/sqrt(sig.d[2]*sig.d[3])

}"












## Independent Outcomes model... Not recommended...  Need to input A.Tau, B.Tau  in place of WA, WD, dfT

model.ind.str <- "model{

  for(i in 1:N){
    y1[i] ~ dpois(lambda[i,1])
  }
  for(i in 1:M){
    y2[i] ~ dpois(rlambda[i,2])
    y3[i] ~ dpois(rlambda[i,3])
  }

  for(i in 1:N){
    for(l in 1:L){
      log(lambda[i,l]) <- A[unit[i],l] - B[l]*Z[i] + X[i,]%*%C[,l] + delta[i,l]
    }
  }
  for(i in 1:M){
    for(l in 1:L){
      rlambda[i,l] <- sum(lambda[robs[i,1]:robs[i,2],l])
    }
  }

 ## Model for A and B
  for(u in 1:U){
    for(l in 1:L){
      A[u,l] <- A.raw[u,l]
    }
    A.raw[u,1:L] ~ dmnorm(mu.a.raw, Tau.a.raw[,])
  }
  for(l in 1:L){
    mu.a.raw[l] <- mu.a[l]
    mu.a[l] ~ dnorm(0,.0001)
    B.raw[l] ~ dlnorm(-1,.1)
    B[l] <- B.raw[l] + offset[l]
    tau.a.vec[l] ~ dgamma(A.Tau,B.Tau)
    Tau.a.raw[l,l] <- tau.a.vec[l]
    tau.d.vec[l] ~ dgamma(A.Tau,B.Tau)
    Tau.d[l,l] <- tau.d.vec[l]
  }
  for(l in 1:(L-1)){
    for(m in (l+1):L){
      Tau.a.raw[l,m] <- 0
      Tau.a.raw[m,l] <- 0
      Tau.d[l,m] <- 0
      Tau.d[m,l] <- 0
    }
  }

  for(j in 1:J){
    for(l in 1:L){
      C[j,l] ~ dnorm(0,.01)
    }
  }

  for(u in 1:U){
    mud[dobs[u,1],1:L] <- rep(0,L)
    for(i in (dobs[u,1]+1):dobs[u,2]){
      mud[i,1:L] <- phi*delta[i-1,1:L]
    }
    for(i in dobs[u,1]:dobs[u,2]){
      delta[i,1:L] ~ dmnorm(mud[i,], Tau.d[,])
    }
    Del.now[u,1:L] <- delta[dobs[u,2],1:L]
  }

  phi ~ dbeta(aa.phi,bb.phi)

 ## distribution for Zc
  for(i in 1:N){
    for(k in 1:K){
      Zc.mat[i,k] ~ dnorm(muZ[i,k], precZ[i])
    }
    Zc[i] <- Zc.mat[i,BS[i]]
    Z[i] <- trunc(Zc[i]+.5)*(Zc[i]>=0.5)*(Zc[i]<4.5) + 5*(Zc[i]>=4.5)
    Zcg0[i] <- Zc[i]*(Zc[i]>0)
  }
  for(k in 1:K){
    p0[k] <- 1/K
  }
  for(u in 1:U){
    for(k in 1:K){
      muZ[dobs[u,1],k] <- muZ0
    }
    precZ[dobs[u,1]] <- precZ0
    BS[dobs[u,1]] ~ dcat(p0)
    for(i in (dobs[u,1]+1):dobs[u,2]){
      for(k in 1:K){
        p.BS[i,k] <- rho[visit[i-1]]/(K-1)*(BS[i-1]!=k) + (1-rho[visit[i-1]])*(BS[i-1]==k)
      }
      BS[i] ~ dcat(p.BS[i,1:K])
      for(k in 1:K){
        muZ[i,k] <- Zc.mat[i-1,k]
      }
      precZ[i] <- 200
    }
  }
  for(i in 1:M){
    meanZ[i] <- mean(Zcg0[robs[i,1]:robs[i,2]])
    step[i] ~ dinterval(meanZ[i],cuts)
  }

  muZ0 ~ dnorm(M.Z, P.Z)
  precZ0 <- .02

  for(k in 1:3){
    rho[k] ~ dbeta(aa.rho[k], bb.rho[k])
  }
  
 
 ## Summarize Tau.d as variances and correlations
  Sig.d <- inverse(Tau.d)
  for(l in 1:L){
    sig.d[l] <- Sig.d[l,l]
  }
  cor.d[1] <- Sig.d[1,2]/sqrt(sig.d[1]*sig.d[2])
  cor.d[2] <- Sig.d[1,3]/sqrt(sig.d[1]*sig.d[3])
  cor.d[3] <- Sig.d[2,3]/sqrt(sig.d[2]*sig.d[3])
}"





get.cred.bands3 <- function(Y, alpha=.05, alpha1=alpha/2, alpha2=alpha/2, tol=1E-6){

  sim.conf.l <- function(gamma, Y){
      band <- apply(Y,2,quantile,prob=c(gamma))
    n.outside <- 0
    for(i in 1:nrow(Y))
      n.outside <- n.outside + any(Y[i,]<band)
    return((n.outside/nrow(Y)-alpha1)^2)
  }
  sim.conf.u <- function(gamma, Y){
      band <- apply(Y,2,quantile,prob=c(1-gamma))
    n.outside <- 0
    for(i in 1:nrow(Y))
      n.outside <- n.outside + any(Y[i,]>band)
    return((n.outside/nrow(Y)-alpha2)^2)
  }
  gamma1 <- robust.optimize(fn=sim.conf.l,par.lim=c(0,alpha1),npar=3,rel.tol=tol,Y=Y)$par
  gamma2 <- robust.optimize(fn=sim.conf.u,par.lim=c(0,alpha2),npar=3,rel.tol=tol,Y=Y)$par
print(gamma1)
print(gamma2)
  band <- apply(Y,2,quantile,prob=c(gamma1,1-gamma2))
  mu <- colMeans(Y)
  return(list(mu=mu, band=band))
}





robust.optimize <- function(fn, par.lim, npar=5, rel.tol=1E-2, ...){

 ## Conduct grid search on par
  par.min <- par.lim[1]
  par.max <- par.lim[2]

  inc <- (par.max - par.min)/(npar-1)
  par.vec <- seq(par.min, par.max, inc)
  par.old <- par.vec
  pen.best <- Inf
  pen.prev <- Inf

  repeat{
    npar.now <- length(par.vec)
    for(j in 1:npar.now){
      par <- par.vec[j]
      pen <- fn(par,...)
      if(pen < pen.best){
        pen.best <- pen
        par.best <- par
      }
    }

    if(inc/(abs(par.best)+1E-6) <= rel.tol)
      break
    else
      pen.prev <- pen.best

   ## create par.vec for next pass
    par.min <- par.best - floor(npar/2)/2*inc
    par.min <- max(par.min, par.lim[1])  
    par.max <- par.best + floor(npar/2)/2*inc
    par.max <- min(par.max, par.lim[2])  
    inc <- inc/2
    par.vec <- seq(par.min, par.max, inc) 
    ind <- which.equal(par.vec, par.old)
    par.vec <- par.vec[!ind]
    if(length(par.vec)==0)
      break
    par.old <- c(par.old, par.vec)
  }

  return(list(par=par.best, pen=pen.best))
}


which.equal <- function(x, y){

  n <- length(x)
  ans <- rep(0,n)
  for(i in 1:n){
    ans[i] <- any(approx.equal(y,x[i]))
  }
  return(as.logical(ans))
}


approx.equal <- function(x, y, tol=1E-9){

  return(abs(x - y) < tol)
}

