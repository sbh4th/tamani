#  program:  tamani-bacon.R
#  task:     run bacon decomposition for TAMANI analyses
#  input:    hh_female_impact_man.dta
#  output:   none
#  project:  TAMANI
#  author:   sam harper \ 2021-10-21

## Assumes the following file structure: 

 ## - r-scripts are located in a folder named "code".
 ## - Derived datasets are in the "data-clean" folder.

## Requires the following packages: "here", "tidyverse", "kableExtra",
##   "stargazer", "ggplot2", "ggeffects"

#-----------------------------------------------------------------------
#  SECTION 0: Load packages and data
#  Notes: */
#-----------------------------------------------------------------------
library(here) 
library(tidyverse)
library(haven)
library(bacondecomp)
library(ggplot2)
library(fixest)
library(glue)
library(modelsummary)
library(TwoWayFEWeights)
library(kableExtra)
devtools::install_github("bcallaway11/did")
library(did)

d <- read_dta(here("data-clean", "hh_female_impact_rpt.dta")) %>%
  mutate(district = as_factor(hh1, levels = "labels"),
         time = as.numeric(delwave) + 1,
         sba_birth = as.numeric(sba_birth))

# restrict to variables for analysis
d_ind <- d %>% 
  select(district, sba_birth, txdel, time) %>%
  mutate(dist_id = as.numeric(district)) %>%
  drop_na() %>%
  group_by(district, dist_id) %>%
  mutate(group = min(if_else(txdel==1, time, 5)),
         g2 = if_else(group==2, 1, 0),
         g3 = if_else(group==3, 1, 0),
         g4 = if_else(group==4, 1, 0)) 

# aggregate dataset
d_sba <- d_ind %>%
  group_by(district, dist_id, time) %>%
  summarise(tsba = sum(sba_birth),
            tpop = n(),
            psba = tsba / tpop, # % SBA
            txdel = mean(txdel),
            group = mean(group))

#--------------------------- End of SECTION 0 --------------------------#


#-----------------------------------------------------------------------
#  SECTION 1: Run the usual two way fixed effects DD model
#  Notes: */
#-----------------------------------------------------------------------

# Get TWFE coefficient from individual-level dataset
twfe <- fixest::feols(sba_birth ~ txdel | district + time, 
                      data = d,
                      weights = NULL, 
                      cluster = ~district)

summary(twfe)

# Compare to pooled sample
twfep <- fixest::feols(psba ~ txdel | district + time, 
                      data = d_sba,
                      weights = ~tpop, 
                      cluster = ~district)

summary(twfep)

# Without weights
twfepu <- fixest::feols(psba ~ txdel | district + time, 
                       data = d_sba,
                       weights =NULL, 
                       cluster = ~district)

# summary of TWFE models
models <- list("Individual" = twfe, "Aggregate" = twfep, 
               "Unweighted" = twfepu)
modelsummary(models, gof_omit = 'DF|Deviance|R2|AIC|BIC|Log.Lik')

#--------------------------- End of SECTION 1 --------------------------#

# create group-time aggregate data
d_agg <- d_sba %>% group_by(group, time) %>%
  summarise(tsba = sum(tsba),
         tpop = sum(tpop),
         psba = tsba / tpop, # % SBA
         txdel = mean(txdel),
         g2 = mean(g2),
         g3 = mean(g3),
         g4 = mean(g4))
  
d_agg %>% select(group, time, psba) %>%
  pivot_wider(names_from = time, values_from = psba) 

d_agg %>% filter(time < 3) %>%
  group_by(g2, time) %>%
  summarise(tsba = sum(tsba),
            tpop = sum(tpop),
            psba = tsba / tpop) %>%
  select(g2, time, psba) %>%
  pivot_wider(names_from = time, values_from = psba) %>%
  mutate(D = `2` - `1`) %>%
  group_by() %>%
  mutate(DD = D - lag(D, default = NA))

d_agg %>% filter(time > 1 & time < 4) %>%
  group_by(g3, time) %>%
  summarise(tsba = sum(tsba),
            tpop = sum(tpop),
            psba = tsba / tpop) %>%
  select(g3, time, psba) %>%
  pivot_wider(names_from = time, values_from = psba) %>%
  mutate(D = `3` - `2`) %>%
  group_by() %>%
  mutate(DD = D - lag(D, default = NA))
  
d_agg %>% filter(time > 2 & time < 5) %>%
  group_by(g4, time) %>%
  summarise(tsba = sum(tsba),
            tpop = sum(tpop),
            psba = tsba / tpop) %>%
  select(g4, time, psba) %>%
  pivot_wider(names_from = time, values_from = psba) %>%
  mutate(D = `4` - `3`) %>%
  group_by() %>%
  mutate(DD = D - lag(D, default = NA))
  


#-----------------------------------------------------------------------
#  SECTION 2: Run the Bacon decomposition
#  Notes: */
#-----------------------------------------------------------------------

df_bacon <- bacon(psba ~ txdel,
                  data = d_sba,
                  id_var = "district",
                  time_var = "time")


# Get de Chaisemartin and D'Haultfoeuille Decomposition
dCDH_decomp <- twowayfeweights(
  df = d_sba, 
  Y = "psba", 
  G = "district",
  T = "time", 
  D ="txdel",
  #weights = "tpop",
  cmd_type =  "feTR"
)

# Weakly Positive weights
dCDH_positive <- sum(dCDH_decomp$weight[dCDH_decomp$weight>=0])

# Negative weights
dCDH_negative <- sum(dCDH_decomp$weight[dCDH_decomp$weight<0])

#--------------------------- End of SECTION 2 --------------------------#




# CS approach
# Use not-yet-treated as comparison group
atts_cs <- did::att_gt(yname = "psba", # name of the LHS variable
                       tname = "time", # name of the time variable
                       idname = "dist_id", # name of the id variable
                       gname = "group", # name of the first treatment period
                       data = d_sba, # name of the data
                       xformla = NULL,
                       weightsname = "tpop",
                       est_method = "reg", # estimation method.
                       control_group = "notyettreated", # set the control group
                       bstrap = TRUE, # if TRUE compute boostrapped SE
                       biters = 1000, # number of boostrap interations
                       print_details = FALSE, # if TRUE, print detailed results
                       panel = FALSE) # panel or repeated cross-sectional
summary(atts_cs)

agg_effects_cs <- aggte(atts_cs, type = "dynamic", min_e = -2, max_e = 2)
summary(agg_effects_cs)

