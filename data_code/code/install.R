# Run in R 4.3.3 for the versions used in the analysis.
options(repos=c(CRAN='https://cloud.r-project.org'))
if(!requireNamespace('ggplot2',quietly=TRUE)) install.packages('ggplot2')
if('refit' %in% commandArgs(TRUE)) {
  for(p in c('Matrix','foreach','sp','lifecycle'))
    if(!requireNamespace(p,quietly=TRUE)) install.packages(p)
  if(!requireNamespace('INLA',quietly=TRUE) || packageVersion('INLA')!=numeric_version('23.03.26'))
    install.packages('https://inla.r-inla-download.org/R/testing/src/contrib/INLA_23.03.26.tar.gz',
                     repos=NULL,type='source')
}
