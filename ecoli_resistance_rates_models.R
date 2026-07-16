
# Author: Noel Waters Feb 2026
# Purpose: Negative binomial regression models for modelling resistance rates in e.coli in
# municipal wastewater treatment plant influents vs hospitals and biofilms vs free flowing water 
# Compares resistance rates between biofilm and wastewater at both hospital and municipal influent sites,
# Also compares rates between biofilms across site types, and between wastewaters across site types.

{
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(DHARMa)
  
}

# Read in data for processing this
# the merged_counts_ecoli can be found as a tab in the Source Data in the submitted article, named "Merged Counts Ecoli"
# In this data, the östra biofilm data have been filtered out as it was not sampled correctly.
df<-readxl::read_xlsx("../for_paper/input/merged_counts_ecoli.xlsx") 




{
# make Site Type effluent influent or hospital
df$Site_Type<-ifelse(df$Site_Type=="Hospital","Hospital",
                     ifelse(grepl("influent| incoming",ignore.case = TRUE,x = df$Site),"Influent",
                            ifelse(grepl("effluent",ignore.case = TRUE,x = df$Site),"Effluent","Fix")
                     )
)


estimates<-df |> group_by(Site,Cleaned_SampleName,plate_type,Site_Type,Sample_Type) |> summarise(
  avg_count=mean(Count),
  median_count=median(Count),
) 

# put the estimates on long format to facilitate choice of estimator to use ( median or mean)

estimators_long<-estimates |>
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

# For the modelling, choose the median of the techincal replicates as input data.
for_nb_simple_model<-estimators_long |> filter(plate_type!="No anti") |> filter(estimator=="median_count") |> filter(Site_Type !="Effluent")
}





estimators_long |> 
  filter(Site_Type=="Hospital" & Sample_Type=="Biofilm") |> 
  filter(estimator=="median_count") |> ungroup() |> 
  dplyr::select(Site,plate_type,CFU_ml) |> pivot_wider(id_cols = Site, names_from = plate_type,values_from = CFU_ml)

# No medan value is not an integer...
any(estimates$median_count %% 1 != 0, na.rm = TRUE)


# For the manuscript, the median of technical replicates is used. 

# plotting the counts by site types and plate type. 
ggplot(estimators_long |> filter(estimator=="median_count"),
       aes(y=CFU_ml,
           x=Site,
           color=Sample_Type, group=Site))+
  geom_point(position=position_dodge(width=0.4)) + theme(axis.text.x = element_text(angle=45, vjust = 0.5)) +
  facet_wrap(~Site_Type*plate_type, scales="free") 




# No effluent samples are to be used in the statistical models. 
# We do not need the no antibiotic plates in the model, since the count on no antibiotic plates is represented in the "reference" column.

# plotting the rates by site type we see biofilm rates higher than muncipal rates for many hospital samples.

ggplot(for_nb_simple_model,
       aes(y=rate,x=Site,color=Sample_Type, group=Site))+
  geom_point(position=position_dodge(width=0.8)) +
  theme(axis.text.x = element_text(angle=90)) +
  facet_wrap(~Site_Type*plate_type, scales="free") 



# individual models for each antibiotic, but still both sample and site types within with the interaction..
large_nb_model_fitter_by_antibiotic<-function(df,antibiotic="CTX", reference_site_type="Influent"){
  
    indat<-df |> filter( plate_type %in% antibiotic)
  
    fit<-glmmTMB(formula = CFU_ml ~ Sample_Type*Site_Type + offset(logref)+ (1 |Site/Cleaned_SampleName), data=indat, family = nbinom2)
    
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
                                            type="link", 
                                            infer=TRUE)
  

  # also dong the fit with no interaction term to compare
  fit_no_int<-glmmTMB(formula = CFU_ml ~ Sample_Type+Site_Type + offset(logref)+ (1 |Site/Cleaned_SampleName), data=indat, family = nbinom2)
  
  
  
  return(list(fit=fit,
              fit_no_int=fit_no_int,
              emms_df=emms_link |> data.frame() |> mutate(plate_type=antibiotic) |> mutate(resistance_rate=exp(emmean),
                                                                                            lower=exp(asymp.LCL),
                                                                                            upper=exp(asymp.UCL)),
              key_stats_by_site_type=contrasts_link_by_site_type |> data.frame() |>
                mutate(plate_type=antibiotic) |>
                mutate(p.adjusted=p.adjust(p.value,method="bonferroni")) |> 
                mutate(ratio=exp(estimate),
                          lower=exp(asymp.LCL),
                          upper=exp(asymp.UCL)),
              key_stats_by_sample_type=contrasts_link_by_sample_type |> data.frame() |>
                mutate(plate_type=antibiotic) |>
                mutate(p.adjusted=p.adjust(p.value,method="bonferroni")) |> 
                mutate(ratio=exp(estimate),
                       lower=exp(asymp.LCL),
                       upper=exp(asymp.UCL)),
              plate_type=antibiotic)
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
#AIC(fits_by_ab[[2]]$fit_no_int,fits_by_ab[[2]]$fit)
anova(fits_by_ab[[2]]$fit_no_int,
      fits_by_ab[[2]]$fit) # a very significant interaction!

# for CTX
a_ctx<-anova(fits_by_ab[[1]]$fit_no_int,
             fits_by_ab[[1]]$fit) # A significant interaction, seen with significant likelihood ratio test-

# for CIP
a_cip<-anova(fits_by_ab[[2]]$fit_no_int,
             fits_by_ab[[2]]$fit) # A significant interaction, seen with significant likelihood ratio test-

anova_df<-rbind(a_ctx |> mutate(Antibiotic="CTX"),a_cip |> mutate(Antibiotic="CIP"))


# Looking at contrasts across sample types in both hospital and municipal :
key_stats_by_site_type_both_abs <- purrr::map_df(fits_by_ab,purrr::pluck("key_stats_by_site_type"))  #|> mutate(p.adjusted=p.adjust(p.value,method="bonferroni"))
key_stats_by_site_type_both_abs

# Looking at contrasts across site types in both biofilms and water:

key_stats_by_sample_type_both_abs <- purrr::map_df(fits_by_ab,purrr::pluck("key_stats_by_sample_type")) 
# All of these are significant after adjustment.
key_stats_by_sample_type_both_abs

# Here we can probably stop for the paper: 

emmeans_by_sample_and_site_type<-purrr::map_df(fits_by_ab,purrr::pluck("emms_df"))
emmeans_by_sample_and_site_type

writexl::write_xlsx(x = list(ecoli_anova=anova_df,
                             ecoli_emmeans=emmeans_by_sample_and_site_type,
                             ecoli_stats_by_site_type=key_stats_by_site_type_both_abs,
                             ecoli_stats_by_sample_type=key_stats_by_sample_type_both_abs,
                             ecoli_median_data=for_nb_simple_model),
                    path = "./output/ecoli_stats.xlsx")



# Below are some visualization of the models, not included in the  --------


for_plot_model_with_interaction_term<-purrr::map_df(fits_by_ab,pluck,"emms_df") |> mutate(ratio=exp(emmean),
                                                                                                               lower=exp(asymp.LCL),
                                                                                                               upper=exp(asymp.UCL)) |> 
  mutate(Sample_Type=relevel(Sample_Type,ref="Water"))

# Visualizing the 
plt<-ggplot(for_plot_model_with_interaction_term , aes(y=ratio,x=Site_Type,color=Sample_Type))+
  geom_errorbar(aes(ymax=upper,ymin=lower), position=position_dodge(width=0.95))+
  geom_point(position=position_dodge(width=0.95))+
  facet_wrap(~plate_type)+ 
  ggtitle(" Estimated marginal mean ecoli resistance rates per Antibiotic, Site Type And Sample Type")+
  ylab("Resistance rate")
# Or on the log scale:
ggplot(for_plot_model_with_interaction_term , aes(y=emmean,x=Site_Type,color=Sample_Type))+
  geom_errorbar(aes(ymax=asymp.UCL,ymin=asymp.LCL), position=position_dodge(width=0.95))+
  geom_point(position=position_dodge(width=0.95))+
  facet_wrap(~plate_type)+ 
  ggtitle(" log(Estimated marginal mean ecoli resistance rates) per Antibiotic, Site Type And Sample Type")+
  ylab("Log(Resistance rate)")


# If we want to add all the individual data points
plt+geom_point(data=for_nb_simple_model,inherit.aes = FALSE,
             mapping = aes(x = Site_Type,y=rate,color=Sample_Type),
             position=position_jitterdodge(jitter.width=0.1),pch=24,size=2,alpha=0.6)


ggplot(for_nb_simple_model, aes(x=plate_type,y=rate,fill=Sample_Type))+
       #geom_point(data=for_nb_simple_model,inherit.aes = FALSE,
        #       mapping = aes(x = Site_Type,y=rate,color=Sample_Type),
  geom_point(
               position=position_jitterdodge(jitter.width=0.1),pch=24,size=2,alpha=0.6)+ facet_wrap(~Site,scales="free")+theme(axis.text.x = element_text(angle=45, vjust = 0.5))


