# VGRegion-Hospital-Wastewater-Selection-Project
On selection toward antibiotic-resistant bacteria from hospital and municipal wastewater in the Västra Götaland Region of west Sweden


To run the scripts for the statistical analysis of resistance rates in e.coli and ARG carriage rates,*ecoli_resistance_rates_models.R* and *metagenomic_arg_models.R* , the following R packages are required along with an R installation ( R version 4.4.0 tested). The scripts require input data avaialble in the submitted manuscript's source data or supplementary information.

| Package   | Version   | Purpose                                         |
|-----------|-----------|-------------------------------------------------|
| tidyverse | 2.0.0     | Data wrangling                                  |
| glmmTMB   | 1.1.14    | Negative binomial mixed effects modelling       |
| emmeans   | 2.0.1    | Estimated marginal means and contrasts          |
| DHARMa    | 0.4.7    | Checking residuals and model fit          |

