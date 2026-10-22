bundle<-readRDS('data/west_nile_background.rds'); cfg<-bundle$config; spde<-cfg$spde
d<-bundle$observations; pres<-d[d$presence==1,]; bg<-d[d$presence==0,]
pop<-bundle$population; cc<-bundle$counties
Ap<-inla.spde.make.A(spde$mesh,loc=as.matrix(pres[,c('Longitude','Latitude')]))
Ab<-inla.spde.make.A(spde$mesh,loc=as.matrix(bg[,c('Longitude','Latitude')]))
Ac<-inla.spde.make.A(spde$mesh,loc=as.matrix(cc[,c('Longitude','Latitude')]))
# Retain the archived Hawaii occurrence and its zero spatial projection.
fit_one <- function(counts, spec, simulation=NULL, theta=NULL) {
  vs <- c(cfg$vars,if(spec=='health_extension') 'health_travel_log')
  j <- which(counts>0)
  if(is.null(simulation)) {
    Xp <- as.matrix(pres[,vs]); AA <- Ap
    y <- rep(1,nrow(pres)); nt <- y
  } else {
    jj <- match(cc$ADMcode, simulation$ADMcode)
    stopifnot(!anyNA(jj))
    k <- simulation$positive_years[jj]; ii <- which(k>0)
    Xp <- as.matrix(cc[ii,vs]); AA <- Ac[ii,,drop=FALSE]
    y <- k[ii]; nt <- k[ii]
  }
  X <- rbind(Xp,as.matrix(bg[j,vs]))
  form <- as.formula(paste('y ~ -1 + Intercept + f(s,model=spde) +',paste(vs,collapse=' + ')))
  sk <- inla.stack(data=list(y=c(y,rep(0,length(j))),Ntrials=c(nt,counts[j])),
      A=list(rbind(AA,Ab[j,,drop=FALSE]),1),
      effects=list(s=seq_len(spde$n.spde),data.frame(Intercept=1,X)),remove.unused=FALSE)
  dat <- inla.stack.data(sk,spde=spde)
  f <- inla(form,family='binomial',data=dat,Ntrials=dat$Ntrials,
      control.fixed=cfg$fixed,
      control.predictor=list(A=inla.stack.A(sk),compute=FALSE,link=1),
      control.compute=list(dic=FALSE,waic=FALSE,cpo=FALSE,config=FALSE,return.marginals=FALSE),
      control.inla=list(strategy='adaptive',cmin=0),
      control.mode=if(is.null(theta)) list() else list(theta=theta,restart=TRUE),
      num.threads='1:1',verbose=FALSE,inla.mode='experimental')
  stopifnot(isTRUE(f$ok),f$mode$mode.status==0)
  z <- as.data.frame(f$summary.fixed[,c('mean','sd','0.025quant','0.975quant')])
  names(z) <- c('mean','sd','lower95','upper95')
  stopifnot(all(is.finite(as.matrix(z))),all(z$sd>0))
  data.frame(parameter=rownames(z),z,row.names=NULL)
}
draw_pair <- function(seed) {
  set.seed(seed); u <- runif(nrow(bg))
  sampler <- function(p) {
    cdf <- cumsum(p); cdf[length(cdf)] <- 1
    tabulate(findInterval(u,c(0,cdf),rightmost.closed=TRUE),nbins=length(p))
  }
  list(log=sampler(pop$log_probability),linear=sampler(pop$linear_probability))
}
fit_job <- function(job,module) {
  counts<-if(job$outcome=='baseline')rep(1L,nrow(bg)) else draw_pair(880928L+job$replicate)[[job$scheme]]
  sim<-if(job$outcome=='joint_null')bundle$null_cases[[job$replicate]] else NULL
  fit_one(counts,job$spec,sim)
}
