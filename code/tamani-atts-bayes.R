#  program:  tamani-atts-bayes.R
#  task:     Bayesian analysi of group-time ATTs
#  input:    hh_female_impact_man.dta
#  output:   none
#  project:  TAMANI
#  author:   sam harper \ 2021-11-19

## Assumes the following file structure: 

 ## - r-scripts are located in a folder named "code".
 ## - Derived datasets are in the "data-clean" folder.

#-----------------------------------------------------------------------
#  SECTION 0: Load packages and data
#  Notes: */
#-----------------------------------------------------------------------
library(here) 
library(tidyverse)
library(tidybayes)
library(haven)
library(ggplot2)
library(modelsummary)
library(kableExtra)
library(brms)
library(bayestestR)
library(modelr)
library(cmdstanr)


d <- read_dta(here("data-clean", "hh_female_impact_rpt.dta")) %>%
  mutate(district = as_factor(hh1, levels = "labels"),
         time = as.numeric(delwave) + 1,
         sba_birth = as.numeric(sba_birth))

# restrict to variables for analysis
d_i <- d %>% 
  select(district, sba_birth, txdel, time, personid) %>%
  mutate(dist_id = as.numeric(district),
         p_id = row_number(personid), # unique indiv ID
         tx = recode_factor(txdel, `0`= "Control",
                            `1` = "Treated")) %>% 
  select(-personid) %>%
  drop_na() %>%
  group_by(district, dist_id) %>%
  mutate(group = min(if_else(txdel==1, time, 5))) 

#--------------------------- End of SECTION 0--------------------------#


#-----------------------------------------------------------------------
#  SECTION 1: Create each of the datasets for the group-time ATTs
#  Notes: */
#-----------------------------------------------------------------------

# ATT(2,2)
di_22 <- d_i %>% filter(time < 3) %>% rename(gp = group) %>%
  mutate(group = if_else(gp==2,1,0),
         att = if_else(group == 1 & time == 2, 1, 0))

# ATT(2,3)
di_23 <- d_i %>% filter((time==1 | time==3) & group!=3) %>%
  rename(gp = group) %>%  
  mutate(group = if_else(gp==2,1,0),
         att = if_else(group == 1 & time == 3, 1, 0))

# ATT(2,4)
di_24 <- d_i %>% filter((time==1 | time==4) & 
  group!=3 & group!=4) %>%
  rename(gp = group) %>%  
  mutate(group = if_else(gp==2,1,0),
         att = if_else(group == 1 & time == 4, 1, 0))

# ATT(3,2)
di_32 <- d_i %>% filter(time < 3 & group != 2) %>%
  rename(gp = group) %>%  
  mutate(group = if_else(gp==3,1,0),
         att = if_else(group == 1 & time == 2, 1, 0))
  
# ATT(3,3)
di_33 <- d_i %>% filter(time >1 & time < 4 & group != 2) %>%
  rename(gp = group) %>%  
  mutate(group = if_else(gp==3,1,0),
         att = if_else(group == 1 & time == 3, 1, 0))

# ATT(3,4)
di_34 <- d_i %>% filter((time == 2 | time == 4) & 
  (group == 3 | group == 5)) %>%
  rename(gp = group) %>% 
  mutate(group = if_else(gp==3,1,0),
         att = if_else(group == 1 & time == 4, 1, 0))

# ATT(4,2)
di_42 <- d_i %>% filter(time < 3 & group >= 3) %>%
  rename(gp = group) %>% 
  mutate(group = if_else(gp==4,1,0),
         att = if_else(group == 1 & time == 2, 1, 0))

# ATT(4,3)
di_43 <- d_i %>% filter((time ==2 | time == 3) & group >= 4) %>%
  rename(gp = group) %>% 
  mutate(group = if_else(gp==4,1,0),
         att = if_else(group == 1 & time == 3, 1, 0))

# ATT(4,4)
di_44 <- d_i %>% filter((time ==3 | time == 4) & group >= 4) %>%
  rename(gp = group) %>% 
  mutate(group = if_else(gp==4,1,0),
         att = if_else(group == 1 & time == 4, 1, 0)) 

#--------------------------- End of SECTION 1--------------------------#


#-----------------------------------------------------------------------
#  SECTION 2: Estimate Bayesian models for Group 2
#  Notes: Flat priors first, then skeptical */
#-----------------------------------------------------------------------

# flat priors
pflat <- c(prior(normal(0, 10), class = Intercept), 
                prior(normal(0, 10), class = b),          
                prior(exponential(1), class = sd))

# regularizing priors 
preg <- c(prior(normal(0, 2), class = Intercept), 
                prior(normal(0, 1), class = b),          
                prior(exponential(1), class = sd))

# ATT(2,2)
b22_flat <- brm(data = di_22, family = bernoulli(),
      sba_birth ~ 1 + (1 | dist_id) + group + time + att,
      prior = pflat, iter = 5000, warmup = 1000,
      chains = 4, cores = 4, seed = 220,
      file = "code/fits/b22_flat")

b22_reg <- brm(data = di_22, family = bernoulli(),
      sba_birth ~ 1 + (1 | dist_id) + group  + time + att,
      prior = preg, iter = 5000, warmup = 1000,
      chains = 4, cores = 4, seed = 221,
      file = "code/fits/b22_reg")


# ATT(2,3)
b23_flat <- brm(data = di_23, family = bernoulli(),
      sba_birth ~ 1 + (1 | dist_id) + group + time + att,
      prior = pflat, iter = 5000, warmup = 1000,
      chains = 4, cores = 4, seed = 230,
      file = "code/fits/b23_flat")

b23_reg <- brm(data = di_23, family = bernoulli(),
      sba_birth ~ 1 + (1 | dist_id) + group + time + att,
      prior = preg, iter = 5000, warmup = 1000,
      chains = 4, cores = 4, seed = 231,
      file = "code/fits/b23_reg")


# ATT(2,4)
b24_flat <- brm(data = di_24, family = bernoulli(),
      sba_birth ~ 1 + (1 | dist_id) + group + time + att,
      prior = pflat, iter = 5000, warmup = 1000,
      chains = 4, cores = 4, seed = 240,
      file = "code/fits/b24_flat.rda")

b24_reg <- brm(data = di_24, family = bernoulli(),
      sba_birth ~ 1 + (1 | dist_id) + group + time + att,
      prior = preg, iter = 5000, warmup = 1000,
      chains = 4, cores = 4, seed = 241,
      file = "code/fits/b24_reg.rda")


# ATT(3,2)
b32_flat <- brm(data = di_32, family = bernoulli(),
      sba_birth ~ 1 + (1 | dist_id) + group + time + att,
      prior = pflat, iter = 5000, warmup = 1000,
      chains = 4, cores = 4, seed = 320,
      file = "code/fits/b32_flat")

b32_reg <- brm(data = di_32, family = bernoulli(),
      sba_birth ~ 1 + (1 | dist_id) + group + time + att,
      prior = preg, iter = 5000, warmup = 1000,
      chains = 4, cores = 4, seed = 321,
      file = "code/fits/b32_reg")


# ATT(3,3)
b33_flat <- brm(data = di_33, family = bernoulli(),
      sba_birth ~ 1 + (1 | dist_id) + group + time + att,
      prior = pflat, iter = 5000, warmup = 1000,
      chains = 4, cores = 4, seed = 330,
      file = "code/fits/b33_flat")

b33_reg <- brm(data = di_33, family = bernoulli(),
      sba_birth ~ 1 + (1 | dist_id) + group + time + att,
      prior = preg, iter = 5000, warmup = 1000,
      chains = 4, cores = 4, seed = 331,
      file = "code/fits/b33_reg")

#--------------------------- End of SECTION 2--------------------------#


#-----------------------------------------------------------------------
#  SECTION 3: Estimate marginal effects
#  Notes: Flat priors first, then skeptical */
#-----------------------------------------------------------------------


# marginal effects for ATT(2,2)
# first for flat priors
mf <- model.frame(di_22)
mf$att <- 0
yhat0b <- fitted(b22_flat, mf,
                 scale = "response", summary = FALSE)
mf$att <- 1
yhat1b <- fitted(b22_flat, mf,
                 scale = "response", summary = FALSE)

describe_posterior(rowMeans(yhat0b))
describe_posterior(rowMeans(yhat1b))
describe_posterior(rowMeans(yhat1b - yhat0b))

# put into a data frame
me_22f <- tibble(param = "ATT(2,2)",
                prior = "(0,10)",
                est = rowMeans(yhat1b - yhat0b))

mf <- model.frame(di_22)
mf$att <- 0
yhat0b <- fitted(b22_reg, mf,
                 scale = "response", summary = FALSE)
mf$att <- 1
yhat1b <- fitted(b22_reg, mf,
                 scale = "response", summary = FALSE)

describe_posterior(rowMeans(yhat0b))
describe_posterior(rowMeans(yhat1b))
describe_posterior(rowMeans(yhat1b - yhat0b))

me_22r <- tibble(param = "ATT(2,2)",
                prior = "(0,0.5)",
                est = rowMeans(yhat1b - yhat0b))


# marginal effects for ATT(2,3)
# first for flat priors
mf <- model.frame(di_23)
mf$att <- 0
yhat0b <- fitted(b23_flat, mf,
                 scale = "response", summary = FALSE)
mf$att <- 1
yhat1b <- fitted(b23_flat, mf,
                 scale = "response", summary = FALSE)

describe_posterior(rowMeans(yhat0b))
describe_posterior(rowMeans(yhat1b))
describe_posterior(rowMeans(yhat1b - yhat0b))

# put into a data frame
me_23f <- tibble(param = "ATT(2,3)",
                prior = "(0,10)",
                est = rowMeans(yhat1b - yhat0b))

mf <- model.frame(di_23)
mf$att <- 0
yhat0b <- fitted(b23_reg, mf,
                 scale = "response", summary = FALSE)
mf$att <- 1
yhat1b <- fitted(b23_reg, mf,
                 scale = "response", summary = FALSE)

describe_posterior(rowMeans(yhat0b))
describe_posterior(rowMeans(yhat1b))
describe_posterior(rowMeans(yhat1b - yhat0b))

me_23r <- tibble(param = "ATT(2,3)",
                prior = "(0,0.5)",
                est = rowMeans(yhat1b - yhat0b))


# plot of marginal effects
me_22f %>% bind_rows(me_22r) %>%
  ggplot(aes(x = est, y = param, color = prior)) +
    stat_pointinterval(position = "dodge")
