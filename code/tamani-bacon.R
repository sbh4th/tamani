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
  select(district, sba_birth, txdel, time, personid) %>%
  mutate(dist_id = as.numeric(district),
         p_id = row_number(personid)) %>% # unique ID for each individual
  drop_na() %>%
  group_by(district, dist_id) %>%
  mutate(group = min(if_else(txdel==1, time, 5))) 

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
                      cluster = NULL)

summary(twfep)

# Without weights
twfepu <- fixest::feols(psba ~ txdel | district + time, 
                       data = d_sba,
                       weights =NULL, 
                       cluster = NULL)

# summary of TWFE models
models <- list("Individual" = twfe, "Aggregate" = twfep, 
               "Unweighted" = twfepu)
modelsummary(models, gof_omit = 'DF|Deviance|R2|AIC|BIC|Log.Lik')

#--------------------------- End of SECTION 1 --------------------------#



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
                       panel = TRUE) # panel or repeated cross-sectional
summary(atts_cs)

agg_effects_cs <- aggte(atts_cs, type = "dynamic", min_e = -2, max_e = 2)
summary(agg_effects_cs)

# Simple aggregation by population size
agg.simple <- aggte(atts_cs, type = "simple")
summary(agg.simple)

# Dynamic effects
agg.es <- aggte(atts_cs, type = "dynamic")
summary(agg.es)
ggdid(agg.es)

# Group-specific effects
agg.gs <- aggte(atts_cs, type = "group")
summary(agg.gs)
ggdid(agg.gs)

# Calendar time effects
agg.ct <- aggte(atts_cs, type = "calendar")
summary(agg.ct)
ggdid(agg.ct)


## event study estimates
# put estimates into a table
es_df <- data.frame(
  type          = "dynamic",
  term = paste0('ATT(', agg_effects_cs$egt, ")"),
  event.time= agg_effects_cs$egt,
  estimate  = agg_effects_cs$att.egt,
  std.error = agg_effects_cs$se.egt,
  conf.low  = agg_effects_cs$att.egt - 
    agg_effects_cs$crit.val.egt * agg_effects_cs$se.egt,
  conf.high = agg_effects_cs$att.egt 
  + agg_effects_cs$crit.val.egt  * agg_effects_cs$se.egt,
  point.conf.low  = agg_effects_cs$att.egt - 
    stats::qnorm(1 - agg_effects_cs$DIDparams$alp/2) * 
    agg_effects_cs$se.egt,
  point.conf.high = agg_effects_cs$att.egt + 
    stats::qnorm(1 - agg_effects_cs$DIDparams$alp/2) * 
    agg_effects_cs$se.egt
)

# Graph it. 
p_es <- ggplot(data = es_df,
  mapping = aes(x = event.time, y = estimate)) +
  geom_vline(xintercept = 0-0.05, color = 'grey', 
             size = 1.2, linetype = "dotted") + 
  geom_ribbon(aes(ymin= point.conf.low, ymax=  point.conf.high), 
              alpha = 0.5, size = 1, fill = "steelblue") +
  geom_ribbon(aes(ymin=  conf.low, ymax =  conf.high), 
              alpha =  0.3, size = 1, fill = "steelblue") +
  geom_line(mapping = aes(x = event.time, y=estimate), 
            colour = "black", size = 0.6, linetype = "dashed") +
  geom_line(size = 1.2, alpha = 2, colour = "darkblue") +
  geom_hline(yintercept = 0, colour="black", size = 0.25,
             linetype = "dotted") +
  xlab('Event time') + ylab("Event-Study Estimate") +
  scale_x_continuous(breaks = c(-5:5)) +
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title = element_text(color="black",  size = 12))+
  theme(plot.title=ggtext::element_markdown(size = 12,
                                            #face = "bold",
                                            color="black",
                                            hjust=0,
                                            lineheight=1.2)) +
  annotate(geom="text", x=3, y=-0.01, 
  label="Average of post-treatment ES coef's:
           0.0758 (0.0124)",
           color="black")

p_es

# individual-level analysis
atts_csi <- did::att_gt(yname = "sba_birth", # name of the LHS variable
                       tname = "time", # name of the time variable
                       idname = "p_id", # name of the id variable
                       gname = "group", # name of the first treatment period
                       data = d_ind, # name of the data
                       xformla = NULL,
                       weightsname = NULL,
                       est_method = "reg", # estimation method.
                       control_group = "notyettreated", # set the control group
                       bstrap = TRUE, # if TRUE compute boostrapped SE
                       biters = 1000, # number of boostrap interations
                       print_details = FALSE, # if TRUE, print detailed results
                       panel = FALSE) # panel or repeated cross-sectional
summary(atts_csi)

agg.simplei <- aggte(atts_csi, type = "simple")
summary(agg.simplei)

agg.dynamici <- aggte(atts_csi, type = "dynamic")
summary(agg.dynamici)
ggdid(agg.dynamici)

agg.cti <- aggte(atts_csi, type = "calendar")
summary(agg.cti)
ggdid(agg.cti)

# individual-level analysis
# accounting for clustering by district
atts_csic <- did::att_gt(yname = "sba_birth", # name of the LHS variable
                        tname = "time", # name of the time variable
                        idname = "p_id", # name of the id variable
                        gname = "group", # name of the first treatment period
                        data = d_ind, # name of the data
                        xformla = NULL,
                        weightsname = NULL,
                        est_method = "reg", # estimation method.
                        control_group = "notyettreated", # set the control group
                        clustervars = "dist_id",
                        bstrap = TRUE, # if TRUE compute boostrapped SE
                        biters = 1000, # number of boostrap iterations
                        print_details = FALSE, # if TRUE, print detailed results
                        panel = FALSE) # panel or repeated cross-sectional
summary(atts_csic)

atts_csac <- did::att_gt(yname = "psba", # name of the LHS variable
                       tname = "time", # name of the time variable
                       idname = "dist_id", # name of the id variable
                       gname = "group", # name of the first treatment period
                       data = d_sba, # name of the data
                       xformla = NULL,
                       weightsname = "tpop",
                       est_method = "reg", # estimation method.
                       control_group = "notyettreated", # set the control group
                       bstrap = FALSE, # if TRUE compute boostrapped SE
                       biters = 1000, # number of boostrap interations
                       print_details = FALSE, # if TRUE, print detailed results
                       panel = TRUE) # panel or repeated cross-sectional
summary(atts_csac)


# Asymptotic inference at individual-level 
atts_csi <- did::att_gt(yname = "sba_birth", # name of the LHS variable
                        tname = "time", # name of the time variable
                        idname = "p_id", # name of the id variable
                        gname = "group", # name of the first treatment period
                        data = d_i, # name of the data
                        xformla = NULL,
                        weightsname = NULL,
                        est_method = "reg", # estimation method.
                        control_group = "notyettreated", # set the control group
                        bstrap = FALSE, # if TRUE compute boostrapped SE
                        biters = 1000, # number of boostrap interations
                        print_details = FALSE, # if TRUE, print detailed results
                        panel = FALSE) # panel or repeated cross-sectional
summary(atts_csi)

cs_est <- rbind(atts_cs$att, atts_csi$att)

# try and calculate the SE by hand for each ATT
# ATT(2,2)
di_22 <- d_ind %>% filter(time < 3) %>%
  mutate(g22 = if_else(group==2,1,0))

twfe_22 <- fixest::feols(sba_birth ~ txdel | g22 + time, 
                      data = di_22,
                      weights = NULL, 
                      cluster = ~district)
summary(twfe_22)

# ATT(2,3)
di_23 <- d_ind %>% filter((time==1 | time==3) & group!=3) %>%
  mutate(g23 = if_else(group==2,1,0))

twfe_23 <- fixest::feols(sba_birth ~ txdel | g23 + time, 
                      data = di_23,
                      weights = NULL, 
                      cluster = ~district)
summary(twfe_23)


# ATT(4,4)
di_44 <- d_ind %>% filter((time ==3 | time == 4) & group >= 4) %>%
  mutate(g44 = if_else(group==4,1,0)) 

di_44 %>%
  group_by(g44, time) %>%
  summarise(tsba = mean(sba_birth))

twfe_44 <- fixest::feols(sba_birth ~ txdel | g44 + time, 
                      data = di_44,
                      weights = NULL, 
                      cluster = NULL)
summary(twfe_44)

# Bayesian implementation for ATT(4,4)
# random effects for each cluster
b44 <-
  brm(data = di_44, 
      family = bernoulli(),
      sba_birth ~ 1 + (1 | dist_id) + txdel + g44 + time,
      prior = c(prior(normal(0, 10), class = Intercept),  # bar alpha
                prior(normal(0, 10), class = b),          # betas
                prior(exponential(1), class = sd)),        # sigma
      iter = 5000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes",
      seed = 13)

plot(b44)
print(b44)


b44f <-
  brm(data = di_44, 
      family = bernoulli(),
      sba_birth ~ 1 + txdel + g44 + time,
      prior = c(prior(normal(0, 10), class = Intercept),  # bar alpha
                prior(normal(0, 10), class = b)),          # betas
      iter = 5000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes",
      seed = 13)

print(b44f)

library(RStata)
options("RStata.StataVersion" = 16)
options("RStata.StataPath"= '/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp')

s_od1 <- '
prtest sba_birth if g44==1, by(time)
'
stata(s_od1, data.in=di_44)
