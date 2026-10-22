modules <- c('dengue','prior','synthetic','west_nile','background')
result_files <- c(dengue='dengue',prior='prior_mesh',synthetic='synthetic',
                 west_nile='west_nile',background='west_nile_background')
read_reference <- function(module) {
  z <- read.csv(file.path('data/results',paste0(result_files[[module]],'.csv')))
  z$id <- switch(module,dengue=z$id,prior=z$variant,
    synthetic=paste(z$regime,z$environment,z$replicate,z$scheme,sep='__'),
    west_nile=paste(z$cell,z$replicate,sep='__'),
    background=paste(z$outcome,z$spec,z$replicate,z$scheme,sep='__'))
  z
}
load_inla <- function() {
  if(!requireNamespace('INLA',quietly=TRUE)) stop('Install INLA: Rscript code/install.R refit')
  if(packageVersion('INLA')!=numeric_version('23.03.26'))
    stop('The supplied models use INLA 23.03.26. See README.md for installation.')
  suppressPackageStartupMessages(library(INLA))
  suppressPackageStartupMessages(library(Matrix))
  inla.setOption(num.threads='1:1')
  Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1)
}
validate_fit <- function(z,expected) {
  stopifnot(nrow(z)==nrow(expected),!anyDuplicated(z$parameter),
    setequal(z$parameter,expected$parameter),all(z$id==expected$id[1]),
    all(is.finite(as.matrix(z[,c('mean','sd','lower95','upper95')]))),
    all(z$sd>0),all(z$lower95<z$upper95))
  invisible(z)
}
atomic_csv <- function(z,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- paste0(path,'.tmp'); write.csv(z,tmp,row.names=FALSE)
  if(!file.rename(tmp,path)) stop('Cannot save ',path)
}
fixed_coefficients <- function(f) {
  stopifnot(isTRUE(f$ok),f$mode$mode.status==0)
  z <- as.data.frame(f$summary.fixed[,c('mean','sd','0.025quant','0.975quant')])
  names(z)<-c('mean','sd','lower95','upper95')
  data.frame(parameter=rownames(z),z,row.names=NULL)
}
