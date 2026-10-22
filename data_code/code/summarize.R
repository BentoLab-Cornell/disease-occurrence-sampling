# All changes are paired differences of posterior means. Across-draw spread
# is reported separately from the marginal credible intervals of each fit.
by_group <- function(z,keys,fun) {
  groups <- split(seq_len(nrow(z)),interaction(z[,keys,drop=FALSE],drop=TRUE,lex.order=TRUE))
  out <- do.call(rbind,lapply(groups,function(i)cbind(z[i[1],keys,drop=FALSE],fun(z[i,,drop=FALSE]))))
  rownames(out)<-NULL; out
}
change_summary <- function(z,keys) by_group(z,keys,function(x)data.frame(
  n_pairs=nrow(x),mean_change=mean(x$change),median_change=median(x$change),
  min_change=min(x$change),max_change=max(x$change),
  q25_change=unname(quantile(x$change,.25)),q75_change=unname(quantile(x$change,.75)),
  draw_sd=sd(x$change),mcse=sd(x$change)/sqrt(nrow(x)),
  negative_changes=sum(x$change<0),positive_changes=sum(x$change>0),
  sign_reversals=sum(x$mean_from*x$mean_to<0)))
coefficient_summary <- function(z,keys) by_group(z,keys,function(x)data.frame(
  n_fits=nrow(x),mean=mean(x$mean),median=median(x$mean),
  min=min(x$mean),max=max(x$mean),
  positive_intervals=sum(x$lower95>0),negative_intervals=sum(x$upper95<0),
  intervals_include_zero=sum(x$lower95<=0 & x$upper95>=0)))
contrast_spec <- data.frame(
  contrast=c('empty_cells_B_minus_A','weighting_C_minus_B','area_D_minus_C','total_D_minus_A'),
  from_arm=c('legacy_area_log_floor','area_log_positive','area_linear_positive','legacy_area_log_floor'),
  to_arm=c('area_log_positive','area_linear_positive','linear_positive','linear_positive'),
  label=c('Exclude empty cells\nB - A','Log to linear\nC - B','Remove area factor\nD - C','Full recipe\nD - A'))
contrasts <- function(z,spec=contrast_spec) do.call(rbind,lapply(seq_len(nrow(spec)),function(i) {
  a<-z[z$arm==spec$from_arm[i],c('seed','parameter','mean','lower95','upper95')]
  b<-z[z$arm==spec$to_arm[i],c('seed','parameter','mean','lower95','upper95')]
  x<-merge(a,b,by=c('seed','parameter'),suffixes=c('_from','_to'))
  stopifnot(nrow(x)>0,!anyDuplicated(x[,c('seed','parameter')]))
  x$contrast<-spec$contrast[i];x$change<-x$mean_to-x$mean_from;x
}))
read_results <- function(module,mode) {
  z<-read_reference(module)
  if(mode=='saved') return(z)
  for(id in unique(z$id)) {
    path<-file.path('output/refits',module,paste0(id,'.csv'))
    if(!file.exists(path)) stop('Refits are incomplete: ',path)
    new<-read.csv(path); old<-z[z$id==id,]; validate_fit(new,old)
    cols<-c('mean','sd','lower95','upper95')
    z[z$id==id,cols]<-new[match(old$parameter,new$parameter),cols]
  }
  z
}
summarize_results <- function(mode='saved') {
  if(!mode %in% c('saved','refitted'))stop('Choose saved or refitted.')
  if(!requireNamespace('ggplot2',quietly=TRUE))stop('Install ggplot2: Rscript code/install.R')
  tables_dir<-file.path('output',mode,'tables');figure_dir<-file.path('output',mode,'figures')
  dir.create(tables_dir,recursive=TRUE,showWarnings=FALSE)
  dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
  save_table<-function(z,name)atomic_csv(z,file.path(tables_dir,paste0(name,'.csv')))
  r<-setNames(lapply(modules,read_results,mode=mode),modules)
  for(m in modules)save_table(r[[m]],paste0(m,'_coefficients'))
  d<-r$dengue; s<-r$synthetic; v<-r$prior
  original_seeds<-c(980019L,980027L,980035L,980043L);new_seeds<-980001:980016
  seeds<-original_seeds;arms<-c(A='legacy_area_log_floor',B='area_log_positive',C='area_linear_positive',D='linear_positive')
  sampling<-d[d$experiment=='sampling',];paired<-contrasts(sampling)
  summary<-change_summary(paired,c('parameter','contrast'))
  save_table(paired,'sampling_paired_changes');save_table(summary,'sampling_summary')
  save_table(coefficient_summary(sampling,c('arm','parameter')),'sampling_coefficients_summary')
  location<-d[d$experiment=='reconstruction',]
  baseline<-merge(location[location$scenario=='archived_all',c('parameter','mean','sd')],
    location[location$scenario=='reconstructed_all',c('parameter','mean','sd')],
    by='parameter',suffixes=c('_archived','_reconstructed'))
  baseline$change<-baseline$mean_reconstructed-baseline$mean_archived
  save_table(baseline,'reconstruction_baseline')
  spec<-contrast_spec[c(2,4),]
  full<-contrasts(sampling,spec);full$model<-'full_model'
  omit<-contrasts(d[d$experiment=='omission',],spec);omit$model<-'omit_five'
  omission_pairs<-rbind(full,omit)
  omission_summary<-change_summary(omission_pairs,c('model','contrast','parameter'))
  save_table(omission_pairs,'omission_paired_changes');save_table(omission_summary,'omission_summary')
  country_order<-c('All countries','Brazil','India','Philippines')
  country_names<-c(BRA='Brazil',IND='India',PHL='Philippines')
  country_pairs<-do.call(rbind,lapply(c('ALL',names(country_names)),function(code){
    z<-if(code=='ALL')sampling else d[d$experiment=='country' & d$country==code,]
    x<-contrasts(z,spec);x$country<-if(code=='ALL')'All countries' else country_names[[code]]
    x$contrast[x$contrast=='weighting_C_minus_B']<-'controlled_C_minus_B';x
  }))
  country_summary<-change_summary(country_pairs,c('country','contrast','parameter'))
  save_table(country_pairs,'country_paired_changes');save_table(country_summary,'country_summary')
  original<-contrasts(sampling,contrast_spec[4,]);original$cohort<-'original_4'
  additional<-contrasts(d[d$experiment=='additional',],contrast_spec[4,]);additional$cohort<-'new_16'
  draw_pairs<-rbind(original,additional); combined<-draw_pairs;combined$cohort<-'combined_20'
  draw_summary<-change_summary(rbind(draw_pairs,combined),c('cohort','parameter'))
  save_table(draw_pairs,'background_draws_paired_changes');save_table(draw_summary,'background_draws_summary')
  draw_coefs<-d[d$experiment %in% c('sampling','additional') & d$arm %in% arms[c('A','D')],]
  draw_coefs$cohort<-ifelse(draw_coefs$experiment=='sampling','original_4','new_16')
  both<-draw_coefs;both$cohort<-'combined_20'
  save_table(coefficient_summary(rbind(draw_coefs,both),c('cohort','arm','parameter')),'background_draws_coefficients_summary')
  save_table(coefficient_summary(s,c('regime','environment','scheme','parameter')),'synthetic_summary')
  syn_pairs<-merge(s[s$scheme=='log',],s[s$scheme=='linear',],
    by=c('regime','environment','replicate','parameter'),suffixes=c('_from','_to'))
  syn_pairs$change<-syn_pairs$mean_to-syn_pairs$mean_from
  save_table(change_summary(syn_pairs,c('regime','environment','parameter')),'synthetic_paired_summary')
  save_table(coefficient_summary(r$west_nile,c('cell','parameter')),'west_nile_summary')
  save_table(coefficient_summary(r$background,c('outcome','spec','scheme','parameter')),'west_nile_background_summary')
  prior_baseline<-v[v$variant=='original',]
  v$change_from_original<-v$mean-prior_baseline$mean[match(v$parameter,prior_baseline$parameter)]
  v$change_in_original_sd<-v$change_from_original/prior_baseline$sd[match(v$parameter,prior_baseline$parameter)]
  v$contains_zero<-v$lower95<=0 & v$upper95>=0
  save_table(v,'prior_coefficients')
  inputs<-readRDS('data/dengue.rds')
  save_table(inputs$model$scaling,'dengue_scaling')
  case_comparison<-do.call(rbind,lapply(inputs$model$vars,function(nm){
    a<-inputs$archived[inputs$archived$presence==1,nm];b<-inputs$reconstructed_cases[,nm]
    ok<-is.finite(a)&is.finite(b);rmse<-sqrt(mean((a[ok]-b[ok])^2))
    data.frame(parameter=nm,n=length(a),archived_finite=sum(is.finite(a)),
      reconstructed_finite=sum(is.finite(b)),n_both_finite=sum(ok),
      missingness_disagreements=sum(is.na(a)!=is.na(b)),correlation=cor(a[ok],b[ok]),
      rmse=rmse,rmse_per_archived_sd=rmse/inputs$model$scaling$sd[match(nm,inputs$model$scaling$parameter)])
  }))
  save_table(case_comparison,'case_value_comparison')
  counts<-do.call(rbind,lapply(names(inputs$country_backgrounds),function(key){
    do.call(rbind,lapply(names(country_names),function(code)data.frame(
      seed=as.integer(sub('_.*','',key)),arm=sub('^[0-9]+_','',key),country=country_names[[code]],
      n_cases=sum(!inputs$country_cases[[paste0('omit_',code)]]),
      n_backgrounds=sum(!inputs$country_backgrounds[[key]][[paste0('omit_',code)]]))))
  }))
  save_table(counts,'country_sample_counts')
  rm(inputs)
  source('code/figures.R',local=TRUE)
  plot_initial(s,sampling[sampling$arm %in% arms[c('A','D')],],v)
  plot_sampling(paired,summary,seeds,arms)
  plot_reconstruction(location,omission_pairs,omission_summary,seeds)
  plot_countries(country_pairs,country_summary,seeds,country_order)
  plot_draws(draw_pairs,draw_summary,original_seeds,new_seeds)
  writeLines(capture.output(sessionInfo()),file.path('output',mode,'session.txt'))
  cat('Tables and nine figures saved in',file.path('output',mode),'\n')
}
