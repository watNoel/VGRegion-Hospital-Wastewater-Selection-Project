
# Author: Noel Waters Feb 2026
# Purpose: Negative binomial regression models for modelling resistance rates in e.coli in
# municipal wastewater treatment plant influents vs hospitals and biofilms vs free flwoing water 
# Compares resistance rates between biofilm and wastewater at both hospital and municipal influent sites,
# Also compares rates betweem biofilms across site types, and between watewaters across site types.

{
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(DHARMa)
  
}


# the merged_counts_ecoli can be found as a tab in the Source Data in the submitted article, named "Merged Counts Ecoli"
# Filter out the  Östra biofilm, since it was not sampled correctly.
df<-readxl::read_xlsx("merged_counts_ecoli.xlsx") |> filter(!(Site=="Östra" & Sample_Type=="Biofilm"))

# make Site Type effluent influent or hospital
df$Site_Type<-ifelse(df$Site_Type=="Hospital","Hospital",
                     ifelse(grepl("influent| incoming",ignore.case = TRUE,x = df$Site),"Influent",
                            ifelse(grepl("effluent",ignore.case = TRUE,x = df$Site),"Effluent","Fix")
                     )
)

# Calculate estimators for use in statistical modelling

estimates<-df |> group_by(Site,Cleaned_SampleName,plate_type,Site_Type,Sample_Type) |> summarise(
                                                                                                 avg_count=mean(Count),
                                                                                                 median_count=median(Count)
)


# put the estimates on long format to facilitate choice of estimator to use ( median or mean)

estimators_long<-estimates|> 
  pivot_longer(col=-c(Site,Cleaned_SampleName,plate_type,Sample_Type,Site_Type),names_to="estimator",values_to = "CFU_ml") |> 
  group_by(Site,Cleaned_SampleName,Sample_Type,Site_Type,estimator) |> mutate(reference=CFU_ml[plate_type=="No anti"],
                                                                              rate=CFU_ml/reference,
                                                                              logref=log(reference),
                                                                              lograte=log(rate)) |> 
  mutate(Sample_Type=as.factor(Sample_Type),
         Sample_Type=relevel(Sample_Type,ref="Water"),
         Site_Type=as.factor(Site_Type),
         Site_Type=relevel(Site_Type,ref="Influent"))
estimators_long



estimators_long |> 
  filter(Site_Type=="Hospital" & Sample_Type=="Biofilm") |> 
  filter(estimator=="median_count") |> ungroup() |> 
  dplyr::select(Site,plate_type,CFU_ml) |> pivot_wider(id_cols = Site, names_from = plate_type,values_from = CFU_ml)



# For the manuscript, the median of technical replicates is used. 

# plotting the counts by site types and plate type. 
ggplot(estimators_long |> filter(estimator=="median_count"),
       aes(y=CFU_ml,
           x=Site,
           color=Sample_Type, group=Site))+
  geom_point(position=position_dodge(width=0.4)) + theme(axis.text.x = element_text(angle=45, vjust = 0.5)) + facet_wrap(~Site_Type*plate_type, scales="free") 




# No effluent samples are to be used in the statistical models. 
# We do not need the no antibiotic plates in the model, since
#the count on no antibiotic plates is represented in the "reference" column.


for_nb_simple_model<-estimators_long |> filter(plate_type!="No anti") |> filter(estimator=="median_count") |> filter(Site_Type !="Effluent")


# plotting the rates by site type..

# here we see that biofilm resistance rates are always higher than wastewater resistance rates, for both antibiotics

ggplot(for_nb_simple_model ,
       aes(y=rate,x=Site,color=Sample_Type, group=Site))+
  geom_point(position=position_dodge(width=0.8)) +
  theme(axis.text.x = element_text(angle=90)) +
  facet_wrap(~Site_Type*plate_type, scales="free") 

# Just looking at no anti plates

ggplot(for_nb_simple_model ,
       aes(y=logref,x=Site,color=Sample_Type))+
  geom_point(position=position_dodge(width=0.4)) + theme(axis.text.x = element_text(angle=45, vjust=0.5)
                                                         ) + facet_wrap(~Site_Type, scales="free")


# look at the effluents, they are really low! maybe this is sriving the statistics a bit much on the municipal end?
# plot average hospital, influent,effluent levels as a cascade?
ggplot(for_nb_simple_model |> filter(Sample_Type=="Water") ,
       aes(y=log(rate),x=Site,color=plate_type, group=Site))+
  geom_point(position=position_dodge(width=0.8)) +
  theme(axis.text.x = element_text(angle=90)) +
  facet_wrap(~plate_type*Site_Type, scales="free_x") 






# individual models for each antibiotic, but still both sample and site types within with the interaction..
large_nb_model_fitter_by_antibiotic<-function(df,antibiotic="CTX", reference_site_type="Influent"){
  
    indat<-df |> filter( plate_type %in% antibiotic)
  
    fit<-glmmTMB(formula = CFU_ml ~ Sample_Type*Site_Type + offset(logref)+ (1 |Site/Cleaned_SampleName), data=indat, family = nbinom2)
  
    # Does not work!
    #fit2<-glmmTMB(formula = CFU_ml ~ Sample_Type*Site_Type + offset(logref)+ (Sample_Type|Site), data=indat, family = nbinom2)
    
    emms_link <- emmeans(fit,
                       ~  Site_Type* Sample_Type,
                       type = "link",
                       offset = 0)
  
    # Resistant fraction = AB / No Anti (with multiplicity adjustment if need be, here we are only performing one comparison therefore set adjust=none)
  contrasts_link_by_site_type <- contrast(emms_link , method = "trt.vs.ctrl", 
                                          ref="Water",
                                          by = c("Site_Type"), 
                                          adjust = "none",
                                          type="link",
                                          infer=TRUE)
  
  contrasts_link_by_site_type
  # to compare hospital vs reference (influent)
  contrasts_link_by_sample_type <- contrast(emms_link ,
                                            method = "trt.vs.ctrl",
                                            ref=reference_site_type,
                                            by = c("Sample_Type"),
                                            adjust = "none",
                                            type="link", infer=TRUE)
  
  # use this if you want to make stuff plottable
  emms_link |> data.frame() |> mutate(ratio=exp(emmean),
                                      lower=exp(asymp.LCL),
                                      upper=exp(asymp.UCL)) 
  
  
  # also dong the fit with no interaction term to compare
  fit_no_int<-glmmTMB(formula = CFU_ml ~ Sample_Type+Site_Type + offset(logref)+ (1 |Site/Cleaned_SampleName), data=indat, family = nbinom2)
  
  
  
  return(list(fit=fit,
              diffs_df=emms_link |> data.frame() |> mutate(plate_type=antibiotic),
              key_stats_by_site=contrasts_link_by_site_type |> data.frame() |> mutate(plate_type=antibiotic),
              key_stats_by_type=contrasts_link_by_sample_type |> data.frame() |> mutate(plate_type=antibiotic),
              plate_type=antibiotic,
              fit_no_int=fit_no_int)
  )
}


fits_by_ab<-purrr::map(c("CTX","CIP"),.f = \(x) {large_nb_model_fitter_by_antibiotic(for_nb_simple_model,x)})


fits_by_ab[[1]]$fit |> summary()
# Look at residual plots for these models 
# 1 is CTX
# 2 is CIP
DHARMa::simulateResiduals(fits_by_ab[[1]]$fit) |> plot()
#DHARMa::simulateResiduals(fits_by_ab[[1]]$fit_no_int) |> plot()


fits_by_ab[[2]]$fit |> summary()
DHARMa::simulateResiduals(fits_by_ab[[2]]$fit) |> plot()

#DHARMa::simulateResiduals(fits_by_ab[[2]]$fit_no_int) |> plot()


# Is the interaction significant? Compare models with and without the interaction term using anova
# The anova results are reported in the manuscript!

# for CTX
#AIC(fits_by_ab[[1]]$fit_no_int,fits_by_ab[[1]]$fit)

anova(fits_by_ab[[1]]$fit_no_int,
      fits_by_ab[[1]]$fit) # A significant interaction


# for CIP
AIC(fits_by_ab[[2]]$fit_no_int,fits_by_ab[[2]]$fit)
anova(fits_by_ab[[2]]$fit_no_int,
      fits_by_ab[[2]]$fit) # a very significant interaction!



key_stats_both_abs <- purrr::map_df(fits_by_ab,purrr::pluck("key_stats_by_site"))
key_stats_both_abs |> mutate(ratio=exp(estimate),
                             lower=exp(asymp.LCL),
                             upper=exp(asymp.UCL))




# Since there was a significant interaction between the effect of sample type (Water or Biofilm) and Site Type ( Hospital or Influent)
# We split into two separate models, one for each site type and antibiotic

simple_nb_model_fitter<-function(df=for_simple_nb_model, antibiotic, site_type){
  #df<-for_nb_simple_model
  #antibiotic="CTX"
  #site_type="Hospital"
  indat<- df |> filter(plate_type==antibiotic & Site_Type==site_type)
  
  fit<-glmmTMB(formula = CFU_ml ~ Sample_Type + offset(logref)+ (1 |Site/Cleaned_SampleName), data=indat, family = nbinom2)
  
  # This one does not work...
  #fit<-glmmTMB(formula = CFU_ml ~ Sample_Type + offset(logref)+ (Sample_Type|Site), data=indat, family = nbinom2)
  
  #
  emm_link <- emmeans(fit, ~ Sample_Type, type="link", offset=0)
  difference_over_sample_type<-contrast(emm_link, method = "trt.vs.ctrl", infer=TRUE)
  
  emm_rate<-emmeans(fit, ~ Sample_Type, type="response", offset=0)
  #rate_over_sample_type<-contrast(emm_rate, method = "trt.vs.ctrl", infer=TRUE)
  #difference_over_sample_type |> data.frame()
  
  return(list(fit=fit,
              Site_Type=site_type,
              plate_type=antibiotic,
              key_stats=difference_over_sample_type |> data.frame() |> mutate(Site_Type=site_type,plate_type=antibiotic),
              diffs_df=emm_link |> data.frame() |> mutate(Site_Type=site_type,plate_type=antibiotic,                                       ),
              rates_df=emm_rate |> data.frame()|> mutate(Site_Type=site_type,plate_type=antibiotic,
              )
  )
  )
  
}


# This is analogous to modelling by site type but instead by sample type. 
#So either comparing hospital biofilm to influent biofilm or comparing hospital ww to influent ww.
simple_nb_model_fitter_by_type<-function(df=for_simple_nb_model, antibiotic, sample_type){
  #df<-for_nb_simple_model
  #antibiotic="CTX"
  #sample_type="Water"
  indat<- df |> filter(plate_type==antibiotic & Sample_Type==sample_type)
  
  fit<-glmmTMB(formula = CFU_ml ~ Site_Type + offset(logref)+ ( 1|Site/Cleaned_SampleName), data=indat, family = nbinom2)
  

  # Estimate the marginal means on the link scale
  emm_link <- emmeans(fit, ~ Site_Type, type="link", offset=0)
  difference_over_site_type<-contrast(emm_link, method = "trt.vs.ctrl", infer=TRUE)
  
  
  emm_rate<-emmeans(fit, ~ Site_Type, type="response", offset=0)

  
  return(list(fit=fit,
              Sample_Type=sample_type,
              plate_type=antibiotic,
              key_stats=difference_over_site_type |> data.frame() |> mutate(Sample_Type=sample_type,plate_type=antibiotic),
              diffs_df=emm_link |> data.frame() |> mutate(Sample_Type=sample_type,plate_type=antibiotic,
              ),
              rates_df=emm_rate |> data.frame()|> mutate(Sample_Type=sample_type,plate_type=antibiotic,
              )
  )
  )
  
}


# Run the four individual models ( two antibiotics, two site types)

combinations<-expand.grid(unique(for_nb_simple_model$Site_Type),c("CTX","CIP"), stringsAsFactors = FALSE)
combinations

all_four_individual_nb<-purrr::map2(combinations$Var1,
                              combinations$Var2,
                              .f = function(.x,.y){simple_nb_model_fitter(df = for_nb_simple_model,
                                                                          site_type = .x,
                                                                          antibiotic = .y)})



# Looking at residual plots for all models:
# No strong deviations
simulateResiduals(all_four_individual_nb[[1]]$fit) |> plot()
simulateResiduals(all_four_individual_nb[[2]]$fit) |> plot()
simulateResiduals(all_four_individual_nb[[3]]$fit) |> plot()
simulateResiduals(all_four_individual_nb[[4]]$fit) |> plot() # This one has some quantile deviations detected, but severity is minor and qq plot is ok.




key_stats_individual_nb<-purrr::map_df(all_four_individual_nb,pluck,"key_stats")
key_stats_individual_nb

for_plot_individual_nb<-purrr::map_df(all_four_individual_nb,pluck,"diffs_df") |> mutate(ratio=exp(emmean),
                                                                                         lower=exp(asymp.LCL),
                                                                                         upper=exp(asymp.UCL)) |> 
  mutate(Sample_Type=relevel(Sample_Type,ref="Water"))




# For the submitted manuscript, the relevant outputs are- again: -----------------------------------------------------------


# models with both site types and sample types.

  # both these yield a highly significant interaction term!
  fits_by_ab[[1]]$fit |> summary()
  simulateResiduals(fits_by_ab[[1]]$fit) |> plot()
  
  fits_by_ab[[2]]$fit |> summary()
  simulateResiduals(fits_by_ab[[2]]$fit) |> plot()
  # Join test for the significance of interaction
  # for CTX
  a_ctx<-anova(fits_by_ab[[1]]$fit_no_int,
        fits_by_ab[[1]]$fit) # A significant interaction, seen with significant likelihood ratio test-
  
  # for CIP
  a_cip<-anova(fits_by_ab[[2]]$fit_no_int,
        fits_by_ab[[2]]$fit) # A significant interaction, seen with significant likelihood ratio test-
  
  anova_df<-rbind(a_ctx |> mutate(Antibiotic="CTX"),a_cip |> mutate(Antibiotic="CIP"))
  
  


# Also, looking at differences between waters across hospital vs municipal and differences between biofilms across hospital vs municipal

  # Comparing water vs water and biofilm vs biofilm based on contrasts calculated from the overall model
  key_stats_by_type<-purrr::map_df(fits_by_ab,purrr::pluck,"key_stats_by_type")
  key_stats_by_type



# Model checks for the four models fitted for each site type and a --------
  
  
  all_four_individual_nb
  # always check residuals with DHARMA..
  simulateResiduals(all_four_individual_nb[[1]]$fit) |> plot()
  simulateResiduals(all_four_individual_nb[[2]]$fit) |> plot()
  simulateResiduals(all_four_individual_nb[[3]]$fit) |> plot()
  simulateResiduals(all_four_individual_nb[[4]]$fit) |> plot()


# Output:

writexl::write_xlsx(x = list(ecoli_anova=anova_df,
                             ecoli_mean_diffs_by_site_type=for_plot_individual_nb,
                             ecoli_stats_by_site_type=key_stats_individual_nb,
                             ecoli_stats_by_sample_type=key_stats_by_type,
                             ecoli_median_data=for_nb_simple_model), path = "./output/ecoli_stats.xlsx")










# Below are additional visualizaitons of the data.


# A figure to visualize the model ( not in paper): ----------------------------
# Note that the bares represent confidence in the mean estimates, wherefore individual points can and are indeed likely to fall outside of the bars.

four_models_plt<-ggplot(for_plot_individual_nb , aes(y=ratio,x=Sample_Type,color=Site_Type))+
  geom_errorbar(aes(ymax=upper,ymin=lower), position=position_dodge(width=0.95))+
  geom_point(position=position_dodge(width=0.95))+
  facet_wrap(~plate_type)+ 
  ggtitle("Ecoli resistance rates From four models, one per antibiotic and Site type \n and model fits for these models")

four_models_plt+geom_point(data=for_nb_simple_model,inherit.aes = FALSE,
                           mapping = aes(x = Sample_Type,y=rate,color=Site_Type),
                           position=position_jitterdodge(jitter.width=0.1),pch=24,size=2,alpha=0.6)



