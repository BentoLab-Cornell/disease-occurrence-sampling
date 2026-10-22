d <- readRDS('data/dengue.rds')
jobs <- read.csv('data/dengue_jobs.csv')
prepare_job <- function(job) {
  vars <- d$model$vars
  cases <- d$reconstructed_cases
  if(job$experiment=='reconstruction') {
    archived_cases <- d$archived[d$archived$presence==1,]
    archived_bg <- d$archived[d$archived$presence==0,]
    bg <- d$reconstructed_archived_backgrounds
    if(job$scenario=='archived_all') { cases<-archived_cases; bg<-archived_bg }
    if(startsWith(job$scenario,'restore_')) {
      restore <- if(job$scenario=='restore_all_unresolved') d$unresolved else sub('^restore_','',job$scenario)
      cases[,restore] <- archived_cases[,restore,drop=FALSE]
      bg[,restore] <- archived_bg[,restore,drop=FALSE]
    }
  } else {
    key <- paste(job$seed,job$arm,sep='_')
    bg <- d$backgrounds[[key]]
    if(job$experiment=='omission') vars <- setdiff(vars,d$unresolved)
    if(job$experiment=='country') {
      flag <- paste0('omit_',job$country)
      cases <- cases[!d$country_cases[[flag]],]
      bg <- bg[!d$country_backgrounds[[key]][[flag]],]
    }
  }
  raw <- rbind(cases[,c('Longitude','Latitude',vars)],bg[,c('Longitude','Latitude',vars)])
  sc <- d$model$scaling[match(vars,d$model$scaling$parameter),]
  X <- sweep(sweep(as.matrix(raw[,vars]),2,sc$center,'-'),2,sc$sd,'/')
  list(y=c(rep(1,nrow(cases)),rep(0,nrow(bg))),X=X,
       loc=as.matrix(raw[,c('Longitude','Latitude')]),vars=vars)
}
fit_job <- function(job,module) {
  mesh <- d$model$mesh
  v <- data.frame(fixed_sd=1,range_threshold=100,sigma_threshold=10,mesh_factor=1)
  if(module=='prior') {
    v <- d$model$variants[d$model$variants$variant==job$variant,]
    a <- list(y=d$archived$presence,X=d$prior_model$X,loc=d$prior_model$loc,vars=d$model$vars)
    A <- d$prior_model$A
    if(v$mesh_factor!=1) {
      edge <- 3.8*v$mesh_factor
      mesh <- inla.mesh.2d(boundary=d$model$boundary,loc.domain=a$loc,
        max.edge=c(1,5)*edge,cutoff=edge/5,offset=c(edge,diff(range(a$loc[,1]))/10))
      A <- inla.spde.make.A(mesh,loc=a$loc)
    }
  } else {
    a <- prepare_job(jobs[jobs$id==job$id,,drop=FALSE])
    A <- inla.spde.make.A(mesh,loc=a$loc)
  }
  spde <- inla.spde2.pcmatern(mesh,prior.range=c(v$range_threshold,.1),
                            prior.sigma=c(v$sigma_threshold,.01))
  sk <- inla.stack(data=list(y=a$y),A=list(A,1,1),
    effects=list(s=seq_len(mesh$n),data.frame(Intercept=1,a$X),iidx=seq_along(a$y)),
    remove.unused=FALSE,tag='est')
  form <- as.formula(paste('y ~ -1 + Intercept + f(s,model=spde) +',paste(a$vars,collapse=' + ')))
  f <- inla(form,data=inla.stack.data(sk,spde=spde),family='binomial',
    control.fixed=list(mean.intercept=0,prec.intercept=1,mean=0,
                       prec=list(Intercept=1,default=1/v$fixed_sd^2)),
    control.predictor=list(A=inla.stack.A(sk),compute=FALSE,link=1),
    control.compute=list(cpo=FALSE,waic=FALSE,dic=FALSE,config=FALSE,return.marginals=FALSE),
    control.inla=list(strategy='adaptive',cmin=0),inla.mode='experimental',
    verbose=FALSE,num.threads='1:1')
  fixed_coefficients(f)
}
