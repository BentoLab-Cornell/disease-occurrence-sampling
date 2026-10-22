# Run from the unpacked package directory.
args <- commandArgs(TRUE)
source('code/common.R')
action <- if(length(args)) args[1] else 'results'
if(action == 'results') {
  source('code/summarize.R')
  summarize_results(if(length(args)>1) args[2] else 'saved')
} else if(action %in% c('refit','jobs')) {
  module <- if(length(args)>1) args[2] else 'all'
  if(!module %in% c(modules,'all')) stop('Unknown analysis: ',module)
  selected <- if(module=='all') modules else module
  for(m in selected) {
    reference <- read_reference(m)
    ids <- unique(reference$id)
    if(length(args)>2) ids <- ids[ids==args[3]]
    if(!length(ids)) stop('Unknown job; use Rscript run.R jobs ',m)
    if(action=='jobs') { cat(paste(m,ids),sep='\n'); next }
    load_inla()
    e <- new.env(parent=globalenv())
    sys.source(file.path('code',paste0(if(m=='prior')'dengue' else m,'.R')),envir=e)
    for(id in ids) {
      expected <- reference[reference$id==id,]
      path <- file.path('output/refits',m,paste0(id,'.csv'))
      if(file.exists(path)) { validate_fit(read.csv(path),expected); next }
      cat('Fitting',m,id,'\n'); flush.console()
      z <- e$fit_job(expected[1,,drop=FALSE],m)
      z$id <- id
      validate_fit(z,expected)
      atomic_csv(z,path)
      cat('Saved',path,'\n'); flush.console()
      gc()
    }
  }
} else stop('Use: results [saved|refitted], jobs [analysis], or refit [analysis] [job-id]')
