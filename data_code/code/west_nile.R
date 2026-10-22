bundle<-readRDS('data/west_nile.rds'); cfg<-bundle$config; cal<-bundle$calibration
cc<-cfg$cc; bg<-cfg$bg; vars<-cfg$vars; spde<-cfg$spde
Ac<-inla.spde.make.A(spde$mesh,loc=as.matrix(cc[,c('Longitude','Latitude')]))
Ab<-inla.spde.make.A(spde$mesh,loc=as.matrix(bg[,c('Longitude','Latitude')]))
standardize <- function(X, w = rep(1, nrow(X))) {
  den <- colSums((!is.na(X))*w)
  mu <- colSums(X*w, na.rm = TRUE)/den
  sd <- sqrt(colSums(sweep(X, 2, mu)^2*w, na.rm = TRUE)/(den-1))
  stopifnot(all(is.finite(sd)), all(sd > 0))
  Z <- sweep(sweep(X, 2, mu), 2, sd, '/')
  Z[is.na(Z)] <- 0 # matches handling in the original recovery experiment
  list(X = Z, sd = sd, mean = mu)
}
fit_binomial <- function(y, trials, X, A, link = 'logit', offset = NULL) {
  vs <- colnames(X)
  form <- as.formula(paste('y ~ -1 + Intercept + f(s,model=spde) +', paste(vs,collapse=' + ')))
  sk <- inla.stack(data=list(y=y, Ntrials=trials), A=list(A,1),
    effects=list(s=seq_len(spde$n.spde), data.frame(Intercept=1,X,check.names=FALSE)),
    remove.unused=FALSE)
  dat <- inla.stack.data(sk, spde=spde)
  a <- list(formula=form, family='binomial', data=dat, Ntrials=dat$Ntrials,
    control.fixed=cfg$fixed, control.predictor=list(A=inla.stack.A(sk),compute=FALSE,link=1),
    control.compute=list(dic=FALSE,waic=FALSE,cpo=FALSE,config=FALSE,return.marginals=FALSE),
    control.inla=list(strategy='adaptive',cmin=0), num.threads='1:1',
    verbose=FALSE, inla.mode='experimental')
  if (link != 'logit') a$control.family <- list(link=link)
  if (!is.null(offset)) a$offset <- offset
  do.call(INLA::inla, a)
}
coefs <- function(fit, rep, cell, vs, sd = rep(1,length(vs))) {
  stopifnot(isTRUE(fit$ok),fit$mode$mode.status==0)
  z <- fit$summary.fixed[vs, c('mean','sd','0.025quant','0.975quant'), drop=FALSE]
  z <- as.data.frame(sweep(as.matrix(z),1,sd,'/'))
  names(z) <- c('mean','sd','lower95','upper95')
  stopifnot(all(is.finite(as.matrix(z))), all(z$sd > 0),
            all(z$lower95 < z$upper95))
  data.frame(replicate=rep,cell=cell,parameter=vs,z,row.names=NULL)
}
fit_presence <- function(k, rep, cell, background = 'archived', health = FALSE) {
  vs <- c(vars, if(health) 'health_travel_log')
  case <- which(k > 0)
  if(background == 'archived') { XB <- as.matrix(bg[,vs]); AB <- Ab } else {
    set.seed(360927L+rep)
    j <- sample.int(nrow(cc), nrow(bg), replace=TRUE)
    XB <- as.matrix(cc[j,vs]); AB <- Ac[j,,drop=FALSE]
  }
  X <- rbind(as.matrix(cc[case,vs]),XB)
  w <- c(k[case],rep(1,nrow(XB)))
  z <- standardize(X,w)
  fit <- fit_binomial(c(k[case],rep(0,nrow(XB))), w, z$X,
                      rbind(Ac[case,,drop=FALSE],AB))
  coefs(fit,rep,cell,vs,z$sd)
}
fit_county <- function(k, rep, cell, link='cloglog', offset=FALSE, health=FALSE) {
  vs <- c(vars, if(health) 'health_travel_log')
  z <- standardize(as.matrix(cc[,vs]))
  fit <- fit_binomial(k,rep(cfg$years,nrow(cc)),z$X,Ac,link=link,
    offset=if(offset) log(cc$population/1e5) else NULL)
  coefs(fit,rep,cell,vs,z$sd)
}
fit_incidence <- function(counts, rep, cell, equal_population=FALSE) {
  a <- cfg$incidence_args
  a$.parent.frame <- environment()
  a$data$y <- counts[cfg$incidence_county_order]
  j <- cfg$incidence_covariate_rows
  a$data$NumCases[j] <- a$data$y
  a$data$Outbreak[j] <- as.integer(a$data$y>0)
  if(equal_population) {
    a$data$population[j] <- mean(cc$population)
    a$data$logpop[j] <- log(mean(cc$population)/1e5)
  }
  a$inla.call <- inla.getOption('inla.call')
  a$num.threads <- '1:1'; a$blas.num.threads <- 1
  a$verbose <- FALSE
  a$control.predictor$compute <- FALSE
  a$control.compute <- list(dic=FALSE,waic=FALSE,cpo=FALSE,config=FALSE,return.marginals=FALSE)
  a$control.mode <- list(theta=cfg$incidence_theta,restart=TRUE)
  fit <- do.call(INLA::inla,a)
  coefs(fit,rep,cell,rownames(cfg$archived_incidence))
}
simulate_cases <- function(rep, design = 'observed') {
  n <- nrow(cc); T <- cfg$years
  if(design == 'joint') {
    set.seed(660927L+rep)
    mu <- cc$population*cal$joint_rate
    p0 <- exp(-cal$size*log1p(mu/cal$size))
    p <- -expm1(-cal$size*log1p(mu/cal$size))
  } else {
    set.seed(260927L+rep)
    mu <- if(design=='equal') rep(-log1p(-cfg$target_positive/(T*n)),n) else cc$population*cal$rate
    p0 <- exp(-mu)
    p <- -expm1(-mu)
  }
  k <- rbinom(n,T,p)
  counts <- numeric(n)
  for(i in which(k>0)) {
    u <- pmin(p0[i]+runif(k[i])*p[i],1-.Machine$double.eps)
    zz <- if(design=='joint') qnbinom(u,size=cal$size,mu=mu[i]) else qpois(u,lambda=mu[i])
    stopifnot(all(is.finite(zz)),all(zz>=1))
    counts[i] <- sum(zz)
  }
  stopifnot(all(counts>=k),all((counts>0)==(k>0)))
  data.frame(ADMcode=cc$ADMcode,positive_years=k,cases=counts)
}
fit_job <- function(job,module) {
  rep<-job$replicate; cell<-job$cell
  if(cell=='observed_data')return(fit_incidence(cc$NumCases,rep,cell))
  design<-if(grepl('^joint',cell))'joint' else if(grepl('equal',cell))'equal' else 'observed'
  sim<-simulate_cases(rep,design); k<-sim$positive_years
  switch(cell,
    incidence_observed=fit_incidence(sim$cases,rep,cell),
    incidence_equal=fit_incidence(sim$cases,rep,cell,equal_population=TRUE),
    equal_archived=fit_presence(k,rep,cell),
    observed_archived=fit_presence(k,rep,cell),
    observed_county=fit_presence(k,rep,cell,background='county'),
    equal_county=fit_presence(k,rep,cell,background='county'),
    county_logit=fit_county(k,rep,cell,link='logit'),
    county_cloglog=fit_county(k,rep,cell),
    county_cloglog_offset=fit_county(k,rep,cell,offset=TRUE),
    health_presence=fit_presence(k,rep,cell,health=TRUE),
    health_offset=fit_county(k,rep,cell,offset=TRUE,health=TRUE),
    joint_presence=fit_presence(k,rep,cell),
    joint_incidence=fit_incidence(sim$cases,rep,cell))
}
