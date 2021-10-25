library(tidyverse)
library(tidybayes)
library(broom)
library(magrittr)
library(patchwork)
library(haven)
library(rethinking)
library(modelsummary)
library(kableExtra)
library(brms)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(cmdstanr)
library(bayestestR)
library(modelr)
library(ggplot2)

d <- read_dta("hh_female_impact_rpt.dta")

d1 <- d %>%
  # trim to analytic variables
  select("hh1","sba_birth","txdel","delwave") %>%
  
  # limit to complete cases
  drop_na() %>% 
  
  # create indicators for wave
  mutate(delwave = factor(delwave),
         txdel = as.integer(txdel),
         sba_birth = as.integer(sba_birth)) 
glimpse(d1)

m1 <-
  brm(data = d1, 
      family = bernoulli(),
      sba_birth ~ 1 + (1 | hh1) + delwave + txdel,
      prior = c(prior(normal(0, 1.5), class = Intercept),  # bar alpha
                prior(normal(0, 0.5), class = b),          # betas
                prior(exponential(1), class = sd)),        # sigma
      chains = 4, cores = 4,
      sample_prior = "yes",
      seed = 56, 
      file = "m1fit")

plot(m1)
print(m1)
mcmc_plot(m1, pars = c("^r_", "^b_", "^sd_")) +
  theme(axis.text.y = element_text(hjust = 0))

mf <- d1
mf$txdel <- 0
yhat0b <- fitted(m1,
                t_0,
                scale = "response", summary = FALSE)
mf$txdel <- 1
yhat1b <- fitted(m1,
                t_1,
                scale = "response", summary = FALSE)

describe_posterior(rowMeans(yhat0b))
describe_posterior(rowMeans(yhat1b))
describe_posterior(rowMeans(yhat1b - yhat0b))

# set up counterfactual populations
t_0 <- d1 %>% mutate(txdel = 0)
t_1 <- d1 %>% mutate(txdel = 1)
t_pop <- bind_rows(t_0, t_1)


mp <- add_epred_draws(m1, newdata=t_pop) 

# posterior predictions for new clusters
post_m1 <- as_draws_df(m1)
glimpse(post_m1)



f <-
  post_m1 %>% 
  pivot_longer(b_b_school_1No:b_b_school_1Yes) %>% 
  mutate(fitted = inv_logit_scaled(b_a_Intercept + value)) %>% 
  mutate(school_1 = factor(str_remove(name, "b_b_school_1"),
                            labels = labels)) %>% 
  select(name:school_1) %>%
  group_by(school_1) %>%
  # note we're using 80% intervals
  mean_qi(fitted, .width = .95)

f

mf_try <- add_epred_draws(m1, newdata=t_pop, draws=1000) %>%
    mutate(txdel = recode_factor(txdel, `0` = "Control", `1` = "Treated")) %>%
  group_by(txdel, .draw) %>%
  summarise(`Pr(y)` = mean(`.epred`))

mf_try %>% group_by(txdel) %>%
  median_qi(`Pr(y)`)

p1 <- mf_try %>% group_by(Tx) %>%
  mutate(q05 = quantile(`Pr(y)`, probs = 0.025),
    q50 = quantile(`Pr(y)`, probs = 0.50 ),
    q95 = quantile(`Pr(y)`, probs = 0.975 )) %>%
  ggplot(aes(x = `Pr(y)`, fill = `Tx`)) +
  geom_density(alpha = 0.25, aes(y = ..scaled..)) + 
  geom_segment(aes(x = q05, xend = q95, y = 0, yend = 0), 
    position = position_nudge(y = -0.025)) +
  geom_point(aes(x = q50, y = 0), 
    position = position_nudge(y = -0.025)) +
  annotate("text", x = 0.25, y = 1.1, label="Control") +
  annotate("text", x = 0.75, y = 1.1, label="Treated") +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous("Density") +
  scale_fill_manual(values = c('#1b9e77','#d95f02')) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y = element_blank(),
        axis.ticks = element_blank())


p2 <- mf_try %>% 
  pivot_wider(names_from = Tx, values_from = `Pr(y)`) %>%
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
  scale_x_continuous(limits=c(0,1)) +
    scale_y_continuous("") +
    theme_classic() + 
  theme(legend.position = "none", axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank())

f2 <- p1 + p2 + plot_layout(widths = c(2, 1))
f2

ggsave(here("output", "u2s-fig2.png"), plot = f2,
       width = 6.5, height = 3)


m2 <-
  brm(data = d1, 
      family = bernoulli(),
      sba_birth ~ 1 + (1 + txdel | hh1) + delwave + txdel,
      prior = c(prior(normal(0, 1.5), class = Intercept),  # bar alpha
                prior(normal(0, 0.5), class = b),          # betas
                prior(exponential(1), class = sd),
                prior(lkj(2), class = cor)),        
      chains = 4, cores = 4,
      sample_prior = "yes",
      seed = 2785, 
      file = "m2fit")

plot(m2)
print(m2)
mcmc_plot(m2, pars = c("^r_", "^b_", "^sd_")) +
  theme(axis.text.y = element_text(hjust = 0))

# compare RE and RC models
m1 <- add_criterion(m1, c("loo", "waic"))
m2 <- add_criterion(m2, c("loo", "waic"))

loo_compare(m1, m2, criterion = "loo") %>% print(simplify = F)

library(did)
# Use here package to facilitate relative paths
library(here)
# Use these for data manipulation, and plots
library(tidyverse)
library(ggplot2)
library(scales)
library(ggrepel)
library(dplyr)
library(bacondecomp) 
library(TwoWayFEWeights)
library(fixest)
library(glue)

# Get TWFE coefficient
twfe <- fixest::feols(sba_birth ~ txdel | hh1 + delwave, 
                      data = d1,
                      weights = NULL, 
                      cluster = ~hh1)

summary(twfe)

d1pool <- d1 %>%
  mutate(delwave = as.numeric(delwave)) %>%
  group_by(hh1) %>%
  mutate(rtx = min(if_else(txdel==1, delwave, 5))) %>%
  group_by(hh1, delwave) %>%
  summarise(tsba = sum(sba_birth),
         pop = n(),
         psba = tsba / pop,
         txdel = mean(txdel),
         rtx = mean(rtx))

df_bacon <- bacon(psba ~ txdel,
                  data = d1pool,
                  id_var = "hh1",
                  time_var = "delwave")

### Control groups: Not yet treated

d2 <- d1 %>%
  mutate(delwave = as.numeric(delwave)) %>%
  group_by(hh1) %>%
  mutate(rtx = min(if_else(txdel==1, delwave, 5)))

# Use not-yet-treated as comparison group
atts_ny <- did::att_gt(yname = "psba", # name of the LHS variable
                       tname = "delwave", # name of the time variable
                       idname = "hh1", # name of the id variable
                       gname = "rtx", # name of the first treatment period variable
                       data = d1pool, # name of the data
                       xformla = NULL,
                       weightsname = "pop",
                       est_method = "reg", # estimation method. "dr" means doubly robust
                       control_group = "notyettreated", # set the control group which is either "nevertreated" or "notyettreated" 
                       bstrap = TRUE, # if TRUE compute boostrapped SE
                       biters = 1000, # number of boostrap interations
                       print_details = TRUE, # if TRUE, print detailed results
                       panel = TRUE) # whether the data is panel or repeated cross-sectional
summary(atts_ny)
agg_effects_es_ny <- aggte(atts_ny, type = "dynamic", min_e = -5, max_e = 5)
summary(agg_effects_es_ny)


# Use not-yet-treated as comparison group (drop all never-treated)
atts_ny2 <- did::att_gt(yname = "psba", # name of the LHS variable
                        tname = "delwave", # name of the time variable
                        idname = "hh1", # name of the id variable
                        gname = "rtx", # name of the first treatment period variable
                        data = d1pool, # name of the data
                        xformla = NULL,
                        weightsname = "pop",
                        est_method = "reg", # estimation method. "dr" means doubly robust
                        control_group = "notyettreated", # set the control group which is either "nevertreated" or "notyettreated" 
                        bstrap = TRUE, # if TRUE compute boostrapped SE
                        biters = 1000, # number of boostrap interations
                        print_details = FALSE, # if TRUE, print detailed results
                        panel = TRUE) # whether the data is panel or repeated cross-sectional
summary(atts_ny2)
agg_effects_es_ny2 <- aggte(atts_ny2, type = "dynamic", min_e = -5, max_e = 5)
summary(agg_effects_es_ny2)


# Now we are ready to go! Let me put all the outputs into a table

event_study_diff_ny <-   data.frame(
  type          = "dynamic",
  term = paste0('ATT(', agg_effects_es_ny$egt, ")"),
  event.time= agg_effects_es_ny$egt,
  estimate  = agg_effects_es_ny$att.egt,
  std.error = agg_effects_es_ny$se.egt,
  conf.low  = agg_effects_es_ny$att.egt - agg_effects_es_ny$crit.val.egt * agg_effects_es_ny$se.egt,
  conf.high = agg_effects_es_ny$att.egt + agg_effects_es_ny$crit.val.egt  * agg_effects_es_ny$se.egt,
  point.conf.low  = agg_effects_es_ny$att.egt - stats::qnorm(1 - agg_effects_es_ny$DIDparams$alp/2) * agg_effects_es_ny$se.egt,
  point.conf.high = agg_effects_es_ny$att.egt + stats::qnorm(1 - agg_effects_es_ny$DIDparams$alp/2) * agg_effects_es_ny$se.egt
)

event_study_diff_ny2 <-   data.frame(
  type          = "dynamic",
  term = paste0('ATT(', agg_effects_es_ny2$egt, ")"),
  event.time= agg_effects_es_ny2$egt,
  estimate  = agg_effects_es_ny2$att.egt,
  std.error = agg_effects_es_ny2$se.egt,
  conf.low  = agg_effects_es_ny2$att.egt - agg_effects_es_ny2$crit.val.egt * agg_effects_es_ny2$se.egt,
  conf.high = agg_effects_es_ny2$att.egt + agg_effects_es_ny2$crit.val.egt  * agg_effects_es_ny2$se.egt,
  point.conf.low  = agg_effects_es_ny2$att.egt - stats::qnorm(1 - agg_effects_es_ny2$DIDparams$alp/2) * agg_effects_es_ny2$se.egt,
  point.conf.high = agg_effects_es_ny2$att.egt + stats::qnorm(1 - agg_effects_es_ny2$DIDparams$alp/2) * agg_effects_es_ny2$se.egt
)


# Set ggplot theme
theme_set(
  #theme_clean() + 
  theme_classic() +
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      legend.background = element_rect(color = "white"),
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.spacing = unit(10, "lines"))
)


# Plots with not-yet-treated
#---------------------------------------------------------------------------------------
# First option
p_es_ny1 <- ggplot(data = event_study_diff_ny,
                   mapping = aes(x = event.time, y = estimate)) +
  geom_vline(xintercept = 0-0.05, color = 'grey', size = 1.2, linetype = "dotted") + 
  geom_ribbon(aes(ymin= point.conf.low, ymax=  point.conf.high), alpha = 0.5, size = 1, fill = "steelblue")+
  geom_ribbon(aes(ymin=  conf.low, ymax =  conf.high), alpha =  0.3, size = 1, fill = "steelblue")+
  geom_line(mapping = aes(x = event.time, y=estimate), colour = "black", size = 0.6, linetype = "dashed") +
  geom_line(size = 1.2, alpha = 2, colour = "darkblue") +
  geom_hline(yintercept = 0, colour="black", size = 0.25, linetype = "dotted")+
  xlab('Event time') +
  ylab("Event-Study Estimate") +
  scale_x_continuous(breaks = c(-5:5)) +
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title = element_text(color="black",  size = 12))+
  theme(plot.title=ggtext::element_markdown(size = 12,
                                            #face = "bold",
                                            color="black",
                                            hjust=0,
                                            lineheight=1.2)
  )+
  annotate(geom="text", x=3, y=-0.01, label="Average of post-treatment ES coef's:
           0.0758 (0.0124)",
           color="black")

p_es_ny1

# This is another option

p_es_ny2<- ggplot(data = event_study_diff_ny,
                  mapping = aes(x = event.time, y = estimate)) +
  geom_line(size = 0.5, alpha = 2, colour = "black") +
  geom_vline(xintercept = 0-0.05, color = 'grey', size = 1.2, linetype = "dotted") + 
  geom_hline(yintercept = 0, colour="black",  linetype = "dotted")+
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), show.legend = FALSE, linetype= 1, size = 1.1,
                  color = "red")+
  geom_pointrange(aes(ymin = point.conf.low, ymax = point.conf.high), show.legend = FALSE, size = 1.1)+
  xlab("Event time") +
  ylab("Event Study Estimate") +
  scale_x_continuous(breaks = c(-5:5)) +
  ylim(range(event_study_diff_ny$conflow,
             event_study_diff_ny$confhigh,
             event_study_diff_ny2$conflow,
             event_study_diff_ny2$confhigh,
             twfe_es$point.conf.low,
             twfe_es$point.conf.high,
             event_study_diff_never$conf.low,
             event_study_diff_never$conf.high
  ))+
  #scale_y_continuous(breaks = seq(-600, 200, 200), limits = c(-700,200))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title = element_text(color="black",  size = 12))+
  theme(plot.title=ggtext::element_markdown(size = 12,
                                            #face = "bold",
                                            color="black",
                                            hjust=0,
                                            lineheight=1.2)
  )+
  annotate(geom="text", x=3, y=-0.01, label="Average of post-treatment ES coef's:
          0.0758 (0.0124)",
           color="black")

p_es_ny2

