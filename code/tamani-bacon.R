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
library(tidybayes)
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
library(RStata)
options("RStata.StataVersion" = 16)
options("RStata.StataPath"= '/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp')

d <- read_dta(here("data-clean", "hh_female_impact_rpt.dta")) %>%
  mutate(district = as_factor(hh1, levels = "labels"),
         time = as.numeric(delwave) + 1,
         sba_birth = as.numeric(sba_birth))

# restrict to variables for analysis
d_ind <- d %>% 
  select(district, sba_birth, txdel, time, personid) %>%
  mutate(dist_id = as.numeric(district),
         p_id = row_number(personid)) %>% # unique ID for each individual
  select(-personid) %>%
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
# Use not-yet-treated as comparison group, aggregate data, panel structure, no clustering, no bootstrap
atts_cs <- did::att_gt(yname = "psba", # name of the LHS variable
                       tname = "time", # name of the time variable
                       idname = "dist_id", # name of the id variable
                       gname = "group", # name of the first treatment period
                       data = d_sba, # name of the data
                       xformla = NULL,
                       weightsname = "tpop",
                       est_method = "reg", # estimation method.
                       control_group = "notyettreated", # set the control group
                       bstrap = FALSE, # if TRUE compute boostrapped SE
                       biters = 1000, # number of boostrap iterations
                       print_details = FALSE, # if TRUE, print detailed results
                       panel = TRUE) # panel or repeated cross-sectional
summary(atts_cs)

# Use not-yet-treated as comparison group, aggregate data, repeated cross section structure, no clustering, no bootstrap
atts_cs_rc <- did::att_gt(yname = "psba", # name of the LHS variable
                       tname = "time", # name of the time variable
                       idname = "dist_id", # name of the id variable
                       gname = "group", # name of the first treatment period
                       data = d_sba, # name of the data
                       xformla = NULL,
                       weightsname = "tpop",
                       est_method = "reg", # estimation method.
                       control_group = "notyettreated", # set the control group
                       bstrap = FALSE, # if TRUE compute boostrapped SE
                       biters = 1000, # number of boostrap iterations
                       print_details = FALSE, # if TRUE, print detailed results
                       panel = FALSE) # panel or repeated cross-sectional
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
                       bstrap = TRUE, # if TRUE compute bootstrapped SE
                       biters = 1000, # number of bootstrap iterations
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

# check the comparison of proportions 'by hand' in Stata
# as well as regression estimates with different functional 
# form and SEs

s_cs22 <- '
prtest sba_birth if g22==1, by(time)
scalar d1 = r(P_diff)
scalar d1_se = r(se_diff)
prtest sba_birth if g22==0, by(time)
scalar d0 = r(P_diff)
scalar d0_se = r(se_diff)
display "ATT(2,2)= " %5.4f d1 - d0 "  SE = " %5.4f sqrt(d1_se^2 + d0_se^2)
reg sba_birth i.txdel i.g22 i.time, robust
reg sba_birth i.txdel i.g22 i.time, vce(cl dist_id)
xtreg sba_birth i.txdel i.g22 i.time, 
qui logit sba_birth i.txdel i.g22 i.time
margins(r.txdel)
qui logit sba_birth i.txdel i.g22 i.time, vce(cl dist_id)
margins(r.txdel)
'
stata(s_cs22, data.in=di_22)

## use the survey setup to adjust for clustering and test equality of proportions

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

twfe_23 <- fixest::feols(sba_birth ~ att | group + time, 
                      data = di_23,
                      weights = NULL, 
                      cluster = ~dist_id)
summary(twfe_23)

# Bayesian implementation for ATT(4,4)
# random effects for each cluster
library(brms)
library(bayestestR)
library(modelr)
library(cmdstanr)
b44_flat <-
  brm(data = di_44, 
      family = bernoulli(),
      sba_birth ~ 1 + (1 | dist_id) + txdel + g44 + time,
      prior = c(prior(normal(0, 10), class = Intercept),  # bar alpha
                prior(normal(0, 10), class = b),          # betas
                prior(exponential(1), class = sd)),        # sigma
      iter = 5000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes",
      seed = 13,
      file = "code/fits/b44_flat.rda")

# if fit already exists
# b44_flat <- readRDS("code/fits/b44_flat.rds")

plot(b44_flat)


# Same model with regularizing priors
b44_reg <-
  brm(data = di_44, 
      family = bernoulli(),
      sba_birth ~ 1 + (1 | dist_id) + txdel + g44 + time,
      prior = c(prior(normal(0, 2), class = Intercept),  # bar alpha
                prior(normal(0, 1), class = b),          # betas
                prior(exponential(1), class = sd)),        # sigma
      iter = 5000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes",
      seed = 24,
      file = "code/fits/b44_reg")

# if fit already exists
# b44_reg <- readRDS("code/fits/b44_reg.rds")

# Visualize prior distributions for baseline value of SBA
bind_rows(prior_draws(b44_flat),
          prior_draws(b44_reg)) %>% 
  mutate(p = inv_logit_scaled(Intercept),
         w = factor(rep(c(10, 2), each = n() / 2),
                    levels = c(10, 2))) %>% 
  
  # plot
  ggplot(aes(x = p, fill = w)) +
  geom_density(size = 0, alpha = 3/4, adjust = 0.1) +
  scale_fill_manual(expression(italic(prior)), values = c("red", "blue")) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(title = expression(alpha%~%Normal(0*", "*italic(prior))),
       x = "prior prob SBA") + theme_bw() +
  ggsave("output/b44_intercept_priors.png", width = 9, height = 6)

# Visualize prior distributions for ATT(4,4)
bind_rows(prior_draws(b44_flat),
          prior_draws(b44_reg)) %>% 
  mutate(p = inv_logit_scaled(Intercept),
         b = inv_logit_scaled(Intercept + b),
         w = factor(rep(c(10, 1), each = n() / 2),
                    levels = c(10, 1))) %>%
  mutate(diff = abs(b - p)) %>%
  
  # plot
  ggplot(aes(x = diff, fill = w)) +
  geom_density(size = 0, alpha = 3/4, adjust = 0.1) +
  scale_fill_manual(expression(italic(prior)), values = c("red", "blue")) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(title = expression(beta%~%Normal(0*", "*italic(prior))),
       x = "prior prob SBA") + theme_bw() +
  ggsave("output/b44_beta_priors.png", width = 9, height = 6)


# Compare models
b44_flat <- add_criterion(b44_flat, c("loo", "waic"))
b44_reg <- add_criterion(b44_reg, c("loo", "waic"))

loo_compare(b44_flat, b44_reg, criterion = "loo") %>% 
  print(simplify = F)
model_weights(b44_flat, b44_reg, weights = "loo") %>% 
  round(digits = 2)


mf_0 <- di_44 %>% mutate(txdel = 0)
mf_1 <- di_44 %>% mutate(txdel = 1)
mf_me <- bind_rows(mf_0, mf_1)

mf_mpp0 <- di_44 %>%
  data_grid(dist_id, txdel, g44, time) %>%
  add_epred_draws(b44_flat, allow_new_levels = TRUE) %>%
  group_by(txdel, .draw) %>%
  summarise(`Pr(y)` = mean(`.epred`)) %>%
  ungroup()

mf_try <- add_epred_draws(b44_flat, newdata=mf_me) %>%
  mutate(txdel = recode_factor(txdel, `0` = "Control", `1` = "Treated")) %>%
  group_by(txdel, .draw) %>%
  summarise(`Pr(y)` = mean(`.epred`))

mf_try %>% group_by(txdel) %>%
  median_qi(`Pr(y)`)

p1 <- mf_try %>% group_by(txdel) %>%
  mutate(q05 = quantile(`Pr(y)`, probs = 0.025),
         q50 = quantile(`Pr(y)`, probs = 0.50 ),
         q95 = quantile(`Pr(y)`, probs = 0.975 )) %>%
  ggplot(aes(x = `Pr(y)`, fill = `txdel`)) +
  geom_density(alpha = 0.25, aes(y = ..scaled..)) + 
  geom_segment(aes(x = q05, xend = q95, y = 0, yend = 0, color = `txdel`), 
               position = position_nudge(y = -0.025)) +
  geom_point(aes(x = q50, y = 0, color = `txdel`), 
             position = position_nudge(y = -0.025)) +
  annotate("text", x = 0.72, y = 1.1, label="Control") +
  annotate("text", x = 0.85, y = 1.1, label="Treated") +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous("Density") +
  scale_fill_manual(values = c('#1b9e77','#d95f02')) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y = element_blank(),
        axis.ticks = element_blank())


p2 <- mf_try %>% 
  pivot_wider(names_from = txdel, values_from = `Pr(y)`) %>%
  mutate(diff = `Treated` - `Control`) %>%
  mutate(q05 = quantile(`diff`, probs = 0.025),
         q50 = quantile(`diff`, probs = 0.50 ),
         q95 = quantile(`diff`, probs = 0.975 )) %>%
  ggplot(aes(x = `diff`)) +
  geom_density(alpha = 0.25, fill='#7570b3',  aes(y = ..scaled..)) + 
  geom_segment(aes(x = q05, xend = q95, y = 0, yend = 0), 
               position = position_nudge(y = -0.025)) +
  geom_point(aes(x = q50, y = 0), 
             position = position_nudge(y = -0.025)) +
  annotate("text", x = 0.5, y = 1.1, label="Difference") +
  scale_x_continuous(limits=c(-0.5,0.5)) +
  scale_y_continuous("") +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank())

f2 <- p1 + p2 + plot_layout(widths = c(2, 1))
f2

ggsave(here("output", "u2s-fig2.png"), plot = f2,
       width = 6.5, height = 3)


# estimate marginal effects
mf <- model.frame(di_44)
mf$txdel <- 0
yhat0b <- fitted(b44_flat,
                 mf,
                 scale = "response", summary = FALSE)
mf$txdel <- 1
yhat1b <- fitted(b44_flat,
                 mf,
                 scale = "response", summary = FALSE)

describe_posterior(rowMeans(yhat0b))
describe_posterior(rowMeans(yhat1b))
describe_posterior(rowMeans(yhat1b - yhat0b))

test2 <- tibble(param = "ATT(2,3)", est = test)
ggplot(test2, aes(y = param, x=est)) + stat_pointinterval()

di_23 %>%
  data_grid(group, time, att) %>%
  add_epred_draws(b23_reg) %>%
  ggplot(aes(x = .epred, y = att)) +
  stat_pointinterval(.width = c(.68, .95)) 

+
  coord_flip() +
  xlab("predicted probability") +
  scale_x_continuous(breaks = seq(0, 0.24, 0.02))


b44_reg %>%
  spread_draws(b_txdel) %>%
  head(10)

bind_rows(spread_draws(b44_flat),
          spread_draws(b44_reg)) %>%
  head(20)



s_cs44 <- '
prtest sba_birth if g44==1, by(time)
scalar d1 = r(P_diff)
scalar d1_se = r(se_diff)
prtest sba_birth if g44==0, by(time)
scalar d0 = r(P_diff)
scalar d0_se = r(se_diff)
display "ATT(4,4)= " d1 - d0
display "ATT(4,4) SE = " sqrt(d1_se^2 + d0_se^2)
reg sba_birth i.txdel i.g44 i.time, robust
reg sba_birth i.txdel i.g44 i.time, vce(cl dist_id)
qui logit sba_birth i.txdel i.g44 i.time
margins(r.txdel)
qui logit sba_birth i.txdel i.g44 i.time, vce(cl dist_id)
margins(r.txdel)
'
stata(s_cs44, data.in=di_44)

s_cs <- '
csdid psba [weight=tpop], ivar(dist_id) time(time) gvar(group) method(reg) notyet
'
stata(s_cs, data.in=d_sba)
