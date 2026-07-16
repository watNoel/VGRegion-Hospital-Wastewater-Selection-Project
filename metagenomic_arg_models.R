# Author: Noel Waters Feb 2026
# Purpose: generalized linear mixed effects models for modelling arg carriage rates in bacteria
# in municipal wastewater treatment plant influents vs hospitals and biofilms vs wastewaters

# Separate analyses are performed for the data with and without the mölndal biofilm sample included
# outputs include statistics used in the publication, such as anova and estimated differences 
# in arg carriages between waters/biofilms in hospitals and municipal influent sites.


{
  library(tidyverse)
  library(DHARMa)
  library(emmeans)
  library(glmmTMB)
}




# Read in arg counts and metadata --------------------------------------------------

# modify paths to the data and file reading function needed! Examples are provided for csv,tsv and xlsx

resfinder_counts<-readxl::read_xlsx("./input/ResFinder_DB_counts_by_group.xlsx")
#read_tsv("./input/ResFinder_DB_counts_by_group.tsv")

#read.csv("../../metagenomics/ResFinder_Counts_Final.csv",header = TRUE) 



# Only use the arg counts, not distinguished by class.
resfinder_counts<-resfinder_counts |> dplyr::select(Sample,Count)


metadata<-readxl::read_xlsx("../../metagenomics/metadata.xlsx") |> rename(Sample_Type=Type,
                                                                          Site_Type=Type_Site)
metadata$Sample_Type[metadata$Sample_Type=="Wastewater"]<-"Water"

#colnames(metadata)
# Should say
#"Sample"            "Total Reads_After" "Sample_ID"         "Sample_Type"       "Site_Type"         "Group"             "Site"              "Date"             
#"Organism"          "geo_loc_name"      "sample_title" 




# merge these two:
resfinder_counts<- resfinder_counts |> merge(metadata,by="Sample")
resfinder_counts<- resfinder_counts |> rename(Total_Reads=`Total Reads_After`) 
resfinder_counts$CPM<-resfinder_counts$Count/resfinder_counts$Total_Reads*1e6






# Check that names are correctly spelled.
resfinder_counts$Site<-str_replace_all(str_to_title(str_replace_all(resfinder_counts$Site,"_"," "))," ","_")
resfinder_counts$Site_Type<-ifelse(resfinder_counts$Site_Type=="Hospital","Hospital",
                                   ifelse(grepl("influent|incoming",ignore.case = TRUE,x = resfinder_counts$Site),"Influent",
                                          ifelse(grepl("effluent",ignore.case = TRUE,x = resfinder_counts$Site),"Effluent","Fix")
                                   )
)


# sum across ALL arg types since we care about the overall
resfinder_totals<-resfinder_counts |>  group_by(Sample_ID,Sample,Sample_Type,Site_Type,Site,Total_Reads) |> 
  summarise(Total_CPM=sum(CPM),
            Total_Count=sum(Count)) |>
  ungroup() |>
  mutate(Total_Millions=Total_Reads/1e6) |> 
  mutate(Site_Type=factor(Site_Type)) |> 
  mutate(Site_Type=relevel(Site_Type,ref="Influent")) |> 
  mutate(Sample_Type=relevel(factor(Sample_Type),ref="Water")) |> ungroup() |> 
  mutate(log2CPM=log(Total_CPM,base=2)) |> 
  mutate(logTotalMillions=log(Total_Millions))





# Getting an overview of the counts data by site type and water type
# We see that 1: mölndal biofilm stands out considerably and there are no effluent wastewater samples sequenced.
# We will, like for the e.coli rates, not 
ggplot(resfinder_totals,aes(x=Site_Type,y=log2CPM,
                              tal_CPM,fill=Sample_Type)) + geom_boxplot()+ geom_point(position=position_jitterdodge(jitter.width=0.1), pch=24,size=2,alpha=0.6)


# Here we see the distribution across site types
ggplot(resfinder_totals ,
       aes(y=log2CPM,x=Site,color=Sample_Type, group=Sample))+
  geom_point(position=position_dodge(width=0.4)) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~Site_Type, scales="free_x")

# We se consistently lower arg carraige in influent biofilms vs waters, and the opposite trend for hospitals.
# With effluents there are no water data so can no say.
# Now to check if differences are statistically significant, let's run a linear mixed effects model that can account 
# for the fact that some samples come from the same sites.




# For the modelling, we drop the effluent samples and do models with and without mölndal the outlying sample.

table(resfinder_totals$Site_Type)

resfinder_totals_no_effluent<-resfinder_totals |> filter(Site_Type!="Effluent")

# the outlying mölndal biofilm sample
mölndal<-resfinder_totals_no_effluent |> filter(Site=="Mölndal" & Sample_Type=="Biofilm") |> pull(Sample)
no_mölndal<- resfinder_totals_no_effluent |> filter(Sample!=mölndal)




# First, A Comparison of effluent and influent biofilms ----------------------------
municipal_biofilms<-resfinder_totals |> filter(Sample_Type=="Biofilm" & Site_Type !="Hospital")
influents<-municipal_biofilms |> filter(Site_Type=="Influent") |> pull(log2CPM)
effluents<-municipal_biofilms |> filter(Site_Type=="Effluent") |> pull(log2CPM)

t.test(influents,effluents, var.equal = FALSE)

rm(municipal_biofilms,influents,effluents)


#unpaired, comparing all to all is highly significant.


# Modelling of arg counts -------------------------------------------------
# Idea is to Model the Total Counts and use the million reads as offset, such that what
#is effectively modelled are the ARG counts per million reads

# This is anologous to the model used for ecoli where it was resistant ecoli per ecoli.
# Since we have counts data again,and some relatedness of samples (some come from the same site) 
# negative binomial regression with random effect for the sites comes in handy.

# 1: Fit one overall model for all site and sample types.
# 2: Check the interaction term for significance, and also compare the resistance rates in biofilm vs water in the two different site types.


nb_fitter<-function(df){
  
  fit1<-glmmTMB(
    Total_Count ~ Sample_Type * Site_Type + (1|Site) + offset(logTotalMillions) ,
    data = df,
    family = nbinom2)  
  
  fit1_no_interaction<-glmmTMB(
    Total_Count ~ Sample_Type + Site_Type + (1|Site) + offset(logTotalMillions), 
    data = df,
    family = nbinom2)  
  
  
  # For more complex random structure. Only use if it improved fit ( assess by AIC and DHARMA)
  fit2<-glmmTMB(
    Total_Count ~ Sample_Type * Site_Type + (1|Site/Sample) + offset(logTotalMillions) , 
    data = df,
    family = nbinom2)  
  
  
  fit2_no_interaction<-glmmTMB(
    Total_Count ~ Sample_Type + Site_Type + (1|Site/Sample) + offset(logTotalMillions), 
    data = df,
    family = nbinom2)  
  return(list(fit1=fit1,
              fit1_no_interaction=fit1_no_interaction,
              fit2=fit2,
              fit2_no_interaction=fit2_no_interaction))
}

nb_fit<-nb_fitter(df=resfinder_totals_no_effluent)
DHARMa::simulateResiduals(nb_fit$fit1) |> plot()
DHARMa::simulateResiduals(nb_fit$fit2) |> plot()
# Same analysis with mölndal biofilm sample dropped -------------------------------------------
nb_fit_no_mölndal<-nb_fitter(df=no_mölndal)
DHARMa::simulateResiduals(nb_fit_no_mölndal$fit1) |> plot()
DHARMa::simulateResiduals(nb_fit_no_mölndal$fit2) |> plot()


# fit1 and fit2 are quite similar, the KS test is significant for fit1.
# Refitting without mölndal outlier removes the deviation and also the single point which deviated on the right side plot.
# This indicates its not really a problem with the model specification, just this outlier is inherently hard to model.

get_aic_and_anova<-function(four_fits){
print(AIC(four_fits$fit1,four_fits$fit1_no_interaction,
    four_fits$fit2,four_fits$fit2_no_interaction))
  
print(anova(four_fits$fit1_no_interaction,four_fits$fit1))

print(anova(four_fits$fit2_no_interaction,four_fits$fit2))  
}

get_aic_and_anova(nb_fit)
get_aic_and_anova(nb_fit_no_mölndal)
# Again, we see the fit1 and fit2 be almost identical in AIC, an highly significant interaction term even when mölndal outlier is included.
# Since the fits are so similar, and AIC is lower for the simpler model, we can proceed with the simpler random effect structure.

get_contrasts_from_fit<-function(fit,reference_site_type="Influent"){
  
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
  
  return(list(fit=fit,
              emms_df=emms_link |> data.frame()  |> mutate(resistance_rate=exp(emmean),
                                                                                            lower=exp(asymp.LCL),
                                                                                            upper=exp(asymp.UCL)),
              key_stats_by_site_type=contrasts_link_by_site_type |> data.frame() |>
                mutate(p.adjusted=p.adjust(p.value,method="bonferroni")) |> 
                mutate(ratio=exp(estimate),
                       lower=exp(asymp.LCL),
                       upper=exp(asymp.UCL)),
              key_stats_by_sample_type=contrasts_link_by_sample_type |> data.frame() |>
                mutate(p.adjusted=p.adjust(p.value,method="bonferroni")) |> 
                mutate(ratio=exp(estimate),
                       lower=exp(asymp.LCL),
                       upper=exp(asymp.UCL))))
}



# Getting the stats to present in the paper
key_stats_full<-get_contrasts_from_fit(fit=nb_fit$fit1)
key_stats_no_mölndal<-get_contrasts_from_fit(fit=nb_fit_no_mölndal$fit1)


# We have significant differences between biofilm and water, regardless of whether the outlier is included in the model or not.

# Save these to an excel file
anova_res<-anova(nb_fit$fit1_no_interaction,nb_fit$fit1) 
anova_res_no_mölndal<-anova(nb_fit_no_mölndal$fit1_no_interaction,nb_fit_no_mölndal$fit1) 


writexl::write_xlsx(list(anova_for_interaction=anova_res,
                         metagenomics_stats=key_stats_full$key_stats_by_site_type,
                         metagenomics_means=key_stats_full$emms_df,
                         anova_for_inter_no_mölndal=anova_res_no_mölndal,
                         metagenomics_stats_no_mölndal=key_stats_no_mölndal$key_stats_by_site_type,
                         metagenomics_means_no_mölndal=key_stats_no_mölndal$emms_df,
                         metagenomics_data=resfinder_totals), 
                    path = "./output/metagenomics_stats.xlsx")





#
