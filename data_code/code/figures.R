suppressPackageStartupMessages(library(ggplot2))
pal<-c('Log weighting'='#3B9AB2','Population proportional'='#F21A00',
       'Known occurrence probability'='#E1AF00') # Wes Anderson Zissou1
th<-theme_minimal(base_size=10,base_family='sans')+theme(
 text=element_text(color='black'),axis.text=element_text(color='black'),
 panel.grid.minor=element_blank(),panel.grid.major.x=element_line(color='grey92',linewidth=.25),
 panel.grid.major.y=element_blank(),strip.text=element_text(face='plain',color='black'),
 legend.title=element_blank(),legend.position='bottom',legend.text=element_text(size=9),
 plot.margin=margin(8,10,8,8))

save_figure<-function(p,name,width,height){
  ggsave(file.path(figure_dir,paste0(name,'.pdf')),p,width=width,height=height,
         device=if(capabilities('cairo'))cairo_pdf else 'pdf')
  ggsave(file.path(figure_dir,paste0(name,'.png')),p,width=width,height=height,dpi=300,bg='white')
}

plot_initial<-function(s,e,v){
# Synthetic null: paired posterior means, not individual observations or intervals.
z<-subset(s,environment=='null'&parameter%in%c('urban','travel'))
stopifnot(nrow(z)==600,!anyDuplicated(z[,c('regime','parameter','replicate','scheme')]))
scheme_labels<-c(log='Log weighting',linear='Population proportional',occurrence='Known occurrence probability')
z$weight<-factor(scheme_labels[z$scheme],levels=names(pal))
z$predictor<-factor(ifelse(z$parameter=='urban','Urban cover','Healthcare travel time'),levels=c('Urban cover','Healthcare travel time'))
z$regime_label<-factor(ifelse(z$regime=='rare','Rare occurrence','Saturating occurrence'),levels=c('Rare occurrence','Saturating occurrence'))
p1<-ggplot(z,aes(weight,mean,group=replicate))+
 geom_hline(yintercept=0,color='grey45',linewidth=.35)+
 geom_line(color='grey70',alpha=.35,linewidth=.25)+
 geom_point(aes(color=weight),size=.75,alpha=.6)+
 stat_summary(aes(group=weight),fun=median,geom='point',shape=23,size=2.5,fill='white',color='black')+
 facet_grid(regime_label~predictor)+scale_color_manual(values=pal)+
 scale_x_discrete(labels=c('Log','Population','Known p'))+
 labs(x=NULL,y='Posterior mean coefficient')+th+theme(panel.grid.major.x=element_blank())
save_figure(p1,'Figure_1_synthetic_null',7.1,3.8)

# Empirical reconstruction: preserve all 4 paired draws and all marginal intervals.
stopifnot(nrow(e)==104,!anyDuplicated(e[,c('seed','arm','parameter')]))
seeds<-sort(unique(e$seed));stopifnot(length(seeds)==4)
e$draw<-match(e$seed,seeds)
e$weight<-factor(ifelse(e$arm=='legacy_area_log_floor','Log weighting','Population proportional'),levels=names(pal)[1:2])
f<-subset(e,parameter%in%c('urban_cover','health_travel_log'))
f$predictor<-factor(ifelse(f$parameter=='urban_cover','Urban cover','Healthcare travel time'),levels=c('Urban cover','Healthcare travel time'))
f$y<-5-f$draw+ifelse(f$arm=='legacy_area_log_floor',.12,-.12)
ln<-merge(subset(f,arm=='legacy_area_log_floor'),subset(f,arm=='linear_positive'),by=c('seed','predictor'),suffixes=c('_log','_pop'))
p2<-ggplot(f,aes(mean,y,color=weight))+
 geom_vline(xintercept=0,color='grey45',linewidth=.35)+
 geom_segment(data=ln,aes(x=mean_log,xend=mean_pop,y=y_log,yend=y_pop),inherit.aes=FALSE,color='grey78',linewidth=.4)+
 geom_segment(aes(x=lower95,xend=upper95,yend=y),linewidth=.65)+
 geom_point(size=2)+facet_wrap(~predictor,scales='free_x')+
 scale_y_continuous(breaks=1:4,labels=4:1,limits=c(.6,4.4))+
 scale_color_manual(values=pal[1:2])+labs(x='Coefficient per archived predictor SD',y='Paired background draw')+th
save_figure(p2,'Figure_2_empirical_focal',7.1,3.5)

# Differences between posterior means. Across-draw ranges are not posterior CIs.
l<-subset(e,arm=='legacy_area_log_floor');r<-subset(e,arm=='linear_positive')
d<-merge(l,r,by=c('seed','parameter'),suffixes=c('_log','_pop'))
d$change<-d$mean_pop-d$mean_log
d$magnitude_reduction_pct<-100*(abs(d$mean_log)-abs(d$mean_pop))/abs(d$mean_log)
d$magnitude_reduction_pct[abs(d$mean_log)<.1]<-NA
write.csv(d,file.path(tables_dir,'empirical_paired_changes.csv'),row.names=FALSE)
sm<-do.call(rbind,lapply(split(d,d$parameter),function(x)data.frame(parameter=x$parameter[1],
 log_median=median(x$mean_log),population_median=median(x$mean_pop),
 paired_median=median(x$change),paired_min=min(x$change),paired_max=max(x$change),
 magnitude_reduction_pct=median(x$magnitude_reduction_pct),
 sign_reversals=sum(x$mean_log*x$mean_pop<0),
 log_intervals_include_zero=sum(x$lower95_log<=0&x$upper95_log>=0),
 population_intervals_include_zero=sum(x$lower95_pop<=0&x$upper95_pop>=0))))
write.csv(sm,file.path(tables_dir,'empirical_summary.csv'),row.names=FALSE)
labs_all<-c(urban_cover='Urban cover',health_travel_log='Healthcare travel time',
 biodiv_intact='Biodiversity intactness',crop_expansion='Cropland expansion',
 evi_dissimilarity='Landscape heterogeneity',forest_cover='Forest cover',forest_loss='Forest loss',
 precip_change='Precipitation change',protected_areas='Protected and conserved areas',
 social_vulnerability='Social vulnerability',tmean_change='Temperature change',urban_expansion='Urban expansion')
q<-subset(sm,parameter!='Intercept');q$label<-factor(labs_all[q$parameter],levels=rev(labs_all))
p3<-ggplot(q,aes(paired_median,label))+
 geom_vline(xintercept=0,color='grey45',linewidth=.35)+
 geom_segment(aes(x=paired_min,xend=paired_max,yend=label),color=pal[1],linewidth=.75)+
 geom_point(color=pal[2],size=2)+labs(x='Paired coefficient change',y=NULL)+th
save_figure(p3,'Figure_8_all_empirical_slopes',7.1,4.0)

# Earlier prior/mesh checks are separate fits using archived extracted covariates.
v<-subset(v,parameter%in%c('urban_cover','health_travel_log'))
variant_names<-c(original='Original settings',fixed_sd_half='Slope prior SD halved',fixed_sd_double='Slope prior SD doubled',
 range_half='Spatial range threshold halved',range_double='Spatial range threshold doubled',
 spatial_sd_half='Spatial SD threshold halved',spatial_sd_double='Spatial SD threshold doubled',finer_mesh='Finer spatial mesh')
stopifnot(nrow(v)==16,all(v$variant%in%names(variant_names)))
v$label<-factor(variant_names[v$variant],levels=rev(variant_names))
v$predictor<-factor(ifelse(v$parameter=='urban_cover','Urban cover','Healthcare travel time'),levels=c('Urban cover','Healthcare travel time'))
v$setting<-ifelse(v$variant=='original','Original','Alternative')
p4<-ggplot(v,aes(mean,label,color=setting))+
 geom_segment(aes(x=lower95,xend=upper95,yend=label),linewidth=.55)+geom_point(size=1.7)+
 facet_wrap(~predictor,scales='free_x')+scale_color_manual(values=c(Original=unname(pal[1]),Alternative=unname(pal[2])))+
 labs(x='Coefficient per archived predictor SD',y=NULL)+th+theme(legend.position='none')
save_figure(p4,'Figure_9_prior_mesh_controls',7.1,3.3)
}

plot_sampling<-function(paired,summaries,seeds,arms){
focal_parameters<-c("urban_cover","health_travel_log")
palette <- c(A='#3B9AB2', B='#78B7C5', C='#E1AF00', D='#F21A00')
predictor_labels <- c(urban_cover='Urban cover', health_travel_log='Healthcare travel time')
plot_theme <- theme_minimal(base_size=11, base_family='sans') + theme(
  text=element_text(color='black'), axis.text=element_text(color='black'),
  panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
  strip.text=element_text(face='plain', color='black'),
  legend.title=element_text(size=10), legend.position='bottom',
  plot.margin=margin(9,12,8,8))
focal_pairs <- paired[paired$parameter %in% focal_parameters, ]
focal_summaries <- summaries[summaries$parameter %in% focal_parameters, ]
for (name in c('focal_pairs','focal_summaries')) {
  z <- get(name)
  z$predictor <- factor(predictor_labels[z$parameter], levels=unname(predictor_labels))
  z$contrast_label <- factor(contrast_spec$label[match(z$contrast, contrast_spec$contrast)],
                            levels=rev(contrast_spec$label))
  z$y <- as.numeric(z$contrast_label)
  assign(name,z)
}
focal_pairs$draw <- factor(match(focal_pairs$seed,seeds), levels=1:4)
focal_pairs$y <- focal_pairs$y + (as.numeric(focal_pairs$draw)-2.5)*.075
contrast_colors <- setNames(c(palette['B'],palette['C'],palette['D'],palette['A']),
                           contrast_spec$contrast)
p2 <- ggplot() +
  geom_vline(xintercept=0, color='grey60', linewidth=.35) +
  geom_segment(data=focal_summaries,
    aes(x=min_change,xend=max_change,y=y,yend=y),color='grey55',linewidth=.7) +
  geom_point(data=focal_pairs, aes(change,y,color=contrast,shape=draw),size=2) +
  geom_point(data=focal_summaries, aes(median_change,y),shape=23,fill='white',color='black',size=2.6) +
  scale_color_manual(values=contrast_colors, guide='none') +
  scale_shape_manual(values=c(16,17,15,18), name='Background draw') +
  scale_y_continuous(breaks=1:4,labels=rev(contrast_spec$label),limits=c(.6,4.4)) +
  facet_wrap(~predictor,scales='free_x') +
  labs(x='Change in posterior mean coefficient',y=NULL) + plot_theme +
  theme(panel.grid.major.y=element_blank(),panel.grid.major.x=element_line(color='grey93',linewidth=.25))
save_figure(p2,'Figure_3_sampling_choices',8.1,4.5)
}

plot_reconstruction<-function(location,paired_contrasts,paired_summary,seeds){
palette<-c(full_model='#3B9AB2',omit_five='#F21A00')
labels<-c(urban_cover='Urban cover',health_travel_log='Healthcare travel time')
plot_theme<-theme_minimal(base_size=11,base_family='sans')+theme(
  text=element_text(color='black'),axis.text=element_text(color='black'),
  panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),
  strip.text=element_text(face='plain',color='black'),legend.position='bottom',
  plot.margin=margin(9,12,8,8))
scenario_labels<-c(archived_all='Archived values',reconstructed_all='Reconstructed values',
  restore_forest_loss='Restore forest loss',restore_precip_change='Restore precipitation change',
  restore_protected_areas='Restore protected areas',restore_social_vulnerability='Restore social vulnerability',
  restore_tmean_change='Restore temperature change',restore_all_unresolved='Restore all five')
focal_locations<-location[location$parameter %in% names(labels),]
focal_locations$predictor<-factor(labels[focal_locations$parameter],levels=unname(labels))
focal_locations$scenario_label<-factor(scenario_labels[focal_locations$scenario],levels=rev(unname(scenario_labels)))
focal_locations$group<-ifelse(focal_locations$scenario=='archived_all','archived',
  ifelse(focal_locations$scenario=='reconstructed_all','reconstructed','restored'))
p1<-ggplot(focal_locations,aes(mean,scenario_label))+
  geom_segment(aes(x=lower95,xend=upper95,yend=scenario_label,color=group),linewidth=.6)+
  geom_point(aes(color=group),size=2.3)+facet_wrap(~predictor,scales='free_x')+
  scale_color_manual(values=c(archived='#3B9AB2',reconstructed='#F21A00',restored='#78B7C5'),guide='none')+
  labs(x='Coefficient per archived predictor SD',y=NULL)+plot_theme
focal_pairs<-paired_contrasts[paired_contrasts$parameter %in% names(labels),]
focal_summary<-paired_summary[paired_summary$parameter %in% names(labels),]
row_levels<-c('Log to linear: full model','Log to linear: omit five','Full recipe: full model','Full recipe: omit five')
for(name in c('focal_pairs','focal_summary')) {
  z<-get(name);z$predictor<-factor(labels[z$parameter],levels=unname(labels))
  z$row_label<-paste0(ifelse(z$contrast=='weighting_C_minus_B','Log to linear: ','Full recipe: '),
    ifelse(z$model=='full_model','full model','omit five'))
  z$y<-as.numeric(factor(z$row_label,levels=rev(row_levels)));assign(name,z)
}
focal_pairs$draw<-factor(match(focal_pairs$seed,seeds),levels=1:4)
focal_pairs$y<-focal_pairs$y+(as.numeric(focal_pairs$draw)-2.5)*.075
p2<-ggplot()+geom_vline(xintercept=0,color='grey60',linewidth=.35)+
  geom_segment(data=focal_summary,aes(x=min_change,xend=max_change,y=y,yend=y),color='grey55',linewidth=.7)+
  geom_point(data=focal_pairs,aes(change,y,color=model,shape=draw),size=2)+
  geom_point(data=focal_summary,aes(median_change,y),shape=23,fill='white',color='black',size=2.6)+
  scale_color_manual(values=palette,guide='none')+scale_shape_manual(values=c(16,17,15,18),name='Background draw')+
  scale_y_continuous(breaks=1:4,labels=rev(row_levels),limits=c(.6,4.4))+
  facet_wrap(~predictor,scales='free_x')+labs(x='Change in posterior mean coefficient',y=NULL)+plot_theme
save_figure(p1,'Figure_4_reconstruction',8.3,5)
save_figure(p2,'Figure_5_omitted_predictors',8.6,4.7)
}

plot_countries<-function(paired_contrasts,contrast_summary,seeds,country_order){
labels<-c(urban_cover='Urban cover',health_travel_log='Healthcare travel time')
focal_pairs<-paired_contrasts[paired_contrasts$parameter%in%names(labels),]
focal_summary<-contrast_summary[contrast_summary$parameter%in%names(labels),]
display_countries<-c('All countries','Omit Brazil','Omit India','Omit Philippines')
for(name in c('focal_pairs','focal_summary')){
  z<-get(name);z$predictor<-factor(labels[z$parameter],levels=unname(labels))
  z$comparison<-factor(ifelse(z$contrast=='controlled_C_minus_B','Log to linear (C - B)','Full recipe (D - A)'),
    levels=c('Log to linear (C - B)','Full recipe (D - A)'))
  z$country<-factor(z$country,levels=country_order)
  z$y<-5-as.numeric(z$country);assign(name,z)
}
focal_pairs$draw<-factor(match(focal_pairs$seed,seeds),levels=1:4)
focal_pairs$y<-focal_pairs$y+(as.numeric(focal_pairs$draw)-2.5)*.075
palette<-setNames(c('#3B9AB2','#78B7C5','#F21A00','#E1AF00'),country_order)
figure<-ggplot()+geom_vline(xintercept=0,color='grey60',linewidth=.35)+
  geom_segment(data=focal_summary,aes(x=min_change,xend=max_change,y=y,yend=y),color='grey55',linewidth=.65)+
  geom_point(data=focal_pairs,aes(change,y,color=country,shape=draw),size=2.1)+
  geom_point(data=focal_summary,aes(median_change,y),shape=23,fill='white',color='black',size=2.7)+
  scale_color_manual(values=palette,guide='none')+scale_shape_manual(values=c(16,17,15,18),name='Background draw')+
  scale_y_continuous(breaks=1:4,labels=rev(display_countries),limits=c(.6,4.4))+
  facet_grid(comparison~predictor,scales='free_x')+labs(x='Change in posterior mean coefficient',y=NULL)+
  theme_minimal(base_size=11,base_family='sans')+theme(text=element_text(color='black'),axis.text=element_text(color='black'),
    panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),strip.text=element_text(face='plain',color='black'),
    strip.text.y=element_text(angle=0),legend.position='bottom',plot.margin=margin(8,10,8,8))
save_figure(figure,'Figure_6_country_omissions',9.3,6.3)
}

plot_draws<-function(paired_contrasts,contrast_summary,original_seeds,new_seeds){
labels<-c(urban_cover='Urban cover',health_travel_log='Healthcare travel time')
focal_pairs<-paired_contrasts[paired_contrasts$parameter%in%names(labels),]
focal_summary<-contrast_summary[contrast_summary$cohort!='combined_20' & contrast_summary$parameter%in%names(labels),]
for(name in c('focal_pairs','focal_summary')){
  z<-get(name);z$predictor<-factor(labels[z$parameter],levels=unname(labels))
  z$cohort<-factor(z$cohort,levels=c('new_16','original_4'))
  z$y<-3-as.numeric(z$cohort);assign(name,z)
}
# Deterministic vertical offsets separate draw points without adding randomness.
focal_pairs$draw_index<-ifelse(focal_pairs$cohort=='new_16',match(focal_pairs$seed,new_seeds),match(focal_pairs$seed,original_seeds))
focal_pairs$cohort_n<-ifelse(focal_pairs$cohort=='new_16',length(new_seeds),length(original_seeds))
focal_pairs$y<-focal_pairs$y+(focal_pairs$draw_index-(focal_pairs$cohort_n+1)/2)*(.35/pmax(1,focal_pairs$cohort_n-1))
palette<-c(new_16='#3B9AB2',original_4='#F21A00')
figure<-ggplot()+geom_vline(xintercept=0,color='grey65',linewidth=.35)+
  geom_segment(data=focal_summary,aes(x=min_change,xend=max_change,y=y,yend=y),color='grey60',linewidth=.6)+
  geom_segment(data=focal_summary,aes(x=q25_change,xend=q75_change,y=y,yend=y),color='grey35',linewidth=1.3)+
  geom_point(data=focal_pairs,aes(change,y,color=cohort),size=2.1,alpha=.85)+
  geom_point(data=focal_summary,aes(mean_change,y),shape=23,fill='white',color='black',size=3)+
  scale_color_manual(values=palette,guide='none')+
  scale_y_continuous(breaks=c(1,2),labels=c('Original 4 pairs','New 16 pairs'),limits=c(.65,2.35))+
  facet_wrap(~predictor,scales='free_x',nrow=1)+labs(x='Change in posterior mean coefficient (D - A)',y=NULL)+
  theme_minimal(base_size=11,base_family='sans')+theme(text=element_text(color='black'),axis.text=element_text(color='black'),
    panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),strip.text=element_text(face='plain',color='black'),
    plot.margin=margin(8,10,8,8))
save_figure(figure,'Figure_7_additional_draws',8.7,4.2)
}
