bundle <- readRDS('data/synthetic.rds')
obj<-bundle$config; cfg<-obj$cfg; spde<-obj$spde; A<-obj$A; grid<-bundle$grid
paired_background <- function(population,n,seed,occurrence_probability=NULL) {
  stopifnot(all(is.finite(population)),all(population>=0),sum(population)>0)
  wlog <- log1p(population)
  wlog[wlog==0] <- median(wlog[wlog>0])
  prob <- list(log=wlog/sum(wlog),linear=population/sum(population))
  if(!is.null(occurrence_probability))
    prob$occurrence <- occurrence_probability/sum(occurrence_probability)
  set.seed(seed); u<-runif(n)
  draws<-lapply(prob,function(p) {
    cdf<-cumsum(p);cdf[length(cdf)]<-1
    tabulate(findInterval(u,c(0,cdf),rightmost.closed=TRUE),nbins=length(p))
  })
  stopifnot(all(vapply(draws,sum,numeric(1))==n))
  list(counts=draws,probability=prob)
}

case_process <- function(regime,environment_beta,replicate,scenario_id) {
  exposure <- grid$population*exp(environment_beta*grid$environment)
  report_rate <- -log1p(-cfg$max_occurrence[[regime]])/max(exposure)
  infection_rate <- report_rate/cfg$report_probability
  p <- -expm1(-exposure*report_rate)
  set.seed(cfg$seed_cases+10000L*scenario_id+replicate)
  k <- as.vector(rmultinom(1,cfg$n_cases,p))
  stopifnot(sum(k)==cfg$n_cases,all(p>=0),all(p<1))
  list(k=k,p=p,report_rate=report_rate,infection_rate=infection_rate)
}

fit_occurrence <- function(k,b,theta=NULL) {
  nt<-k+b; y<-k; y[nt==0]<-NA; nt[nt==0]<-1L
  sk<-inla.stack(data=list(y=y,Ntrials=nt),A=list(A,1),
    effects=list(s=seq_len(spde$n.spde),
                 data.frame(Intercept=1,grid[,cfg$predictors])),remove.unused=FALSE)
  dd<-inla.stack.data(sk,spde=spde)
  f<-inla(y~-1+Intercept+urban+travel+environment+f(s,model=spde),
    family='binomial',data=dd,Ntrials=dd$Ntrials,
    control.fixed=cfg$fixed_prior,
    control.predictor=list(A=inla.stack.A(sk),compute=FALSE,link=1),
    control.compute=list(dic=FALSE,waic=FALSE,cpo=FALSE,config=FALSE,return.marginals=FALSE),
    control.inla=list(strategy='adaptive',cmin=0),
    control.mode=if(is.null(theta))list() else list(theta=theta,restart=TRUE),
    num.threads='1:1',verbose=FALSE,inla.mode='experimental')
  stopifnot(isTRUE(f$ok),f$mode$mode.status==0)
  z<-as.data.frame(f$summary.fixed[,c('mean','sd','0.025quant','0.975quant')])
  names(z)<-c('mean','sd','lower95','upper95')
  stopifnot(all(is.finite(as.matrix(z))),all(z$sd>0))
  z<-data.frame(parameter=rownames(z),z,row.names=NULL)
  attr(z,'theta')<-f$mode$theta
  z
}

fit_job <- function(job,module) {
  scenarios <- expand.grid(regime=c('rare','saturating'),environment=c('null','positive'))
  sid <- which(scenarios$regime==job$regime & scenarios$environment==job$environment)
  beta <- cfg$environmental_effect[[job$environment]]
  cp <- case_process(job$regime,beta,job$replicate,sid)
  b <- paired_background(grid$population,cfg$n_cases*cfg$background_ratio,
    cfg$seed_background+10000L*sid+job$replicate,
    if(job$environment=='null')cp$p else NULL)
  fit_occurrence(cp$k,b$counts[[job$scheme]])
}
