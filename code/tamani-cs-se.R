#  program:  tamani-bacon.R
#  task:     run bacon decomposition for TAMANI analyses
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
#  SECTION 1: 
#  Notes: */
#-----------------------------------------------------------------------


# CS approach
# Use not-yet-treated as comparison group, aggregate data, panel structure, no clustering, no bootstrap
atts_1 <- did::att_gt(yname = "psba", tname = "time", 
  idname = "dist_id", gname = "group", data = d_sba, 
  xformla = NULL, weightsname = "tpop", est_method = "reg",
  control_group = "notyettreated", bstrap = FALSE, 
  biters = 1000, print_details = FALSE, panel = TRUE)

# Use not-yet-treated as comparison group, aggregate data, panel structure, no clustering, bootstrap
atts_2 <- did::att_gt(yname = "psba", tname = "time", 
  idname = "dist_id", gname = "group", data = d_sba, 
  xformla = NULL, weightsname = "tpop", est_method = "reg",
  control_group = "notyettreated", bstrap = TRUE, 
  biters = 1000, print_details = FALSE, panel = TRUE)

# Use not-yet-treated as comparison group, aggregate data, panel structure, clustering, bootstrap
atts_3 <- did::att_gt(yname = "psba", tname = "time", 
  idname = "dist_id", gname = "group", data = d_sba, 
  xformla = NULL, weightsname = "tpop", est_method = "reg",
  control_group = "notyettreated", bstrap = TRUE, 
  biters = 1000, print_details = FALSE, panel = TRUE, 
  clustervars = "dist_id")

# Use not-yet-treated as comparison group, aggregate data, repeated cross section structure, no clustering, no bootstrap
atts_4 <- did::att_gt(yname = "psba", tname = "time", 
  idname = "dist_id", gname = "group", data = d_sba, 
  xformla = NULL, weightsname = "tpop", est_method = "reg",
  control_group = "notyettreated", bstrap = FALSE, 
  biters = 1000, print_details = FALSE, panel = FALSE)

# Use not-yet-treated as comparison group, aggregate data, repeated cross section structure, no clustering, bootstrap
atts_5 <- did::att_gt(yname = "psba", tname = "time", 
  idname = "dist_id", gname = "group", data = d_sba, 
  xformla = NULL, weightsname = "tpop", est_method = "reg",
  control_group = "notyettreated", bstrap = TRUE, 
  biters = 1000, print_details = FALSE, panel = FALSE)

# Use not-yet-treated as comparison group, aggregate data, repeated cross section structure, clustering, bootstrap
atts_6 <- did::att_gt(yname = "psba", tname = "time", 
  idname = "dist_id", gname = "group", data = d_sba, 
  xformla = NULL, weightsname = "tpop", est_method = "reg",
  control_group = "notyettreated", bstrap = TRUE, 
  biters = 1000, print_details = FALSE, panel = FALSE, 
  clustervars = "dist_id")

# df <- vector("list", length(atts_1))
# x <- 1:6
# for (i in 1:6) {
#   z <- i
#   mean(atts_z$att)
# }
# 
# for (i in seq_along(x)) {
#   n <- sample(100, 1)
#   out[[i]] <- rnorm(n, means[[i]])
# }
# , 
#     type = ifelse(atts_[i]$DIDparams$panel,"panel","rc"), 
#     method = ifelse(atts_[i]$DIDparams$clustervars=="dist_id",
#                     ifelse(atts_[i]$DIDparams$bstrap,"bootc",
#                            "bootc"), "asm"), 
#     att = atts_[i]$att, se = atts_[i]$se)
# }
# 

df1 <- data.frame(group = atts_1$group, 
  time = atts_1$t, type = "panel", method = "asm", 
  att = atts_1$att, se = atts_1$se)

df2 <- data.frame(group = atts_2$group, 
  time = atts_2$t, type = "panel", method = "boot", 
  att = atts_2$att, se = atts_2$se)

df3 <- data.frame(group = atts_3$group, 
  time = atts_3$t, type = "panel", method = "bootc", 
  att = atts_2$att, se = atts_2$se)

df4 <- data.frame(group = atts_4$group, 
  time = atts_4$t, type = "rc", method = "asm", 
  att = atts_4$att, se = atts_4$se)

df5 <- data.frame(group = atts_5$group, 
  time = atts_5$t, type = "rc", method = "boot", 
  att = atts_5$att, se = atts_5$se)

df6 <- data.frame(group = atts_6$group, 
  time = atts_6$t, type = "rc", method = "bootc", 
  att = atts_6$att, se = atts_6$se)

atts <- bind_rows(df1,df2,df3,df4,df5,df6)

atts %>% 
  select(-se) %>% 
  pivot_wider(names_from=c(type, method), values_from = att) %>%
  mutate_if(is.numeric, round, 3)


## Individual-level data
# Use not-yet-treated as comparison group, individual data, repeated cross section structure, no clustering, no bootstrap
atts_7 <- did::att_gt(yname = "sba_birth", tname = "time", 
  idname = "p_id", gname = "group", data = d_ind, 
  xformla = NULL, weightsname = NULL, est_method = "reg",
  control_group = "notyettreated", bstrap = FALSE, 
  biters = 1000, print_details = FALSE, panel = FALSE)

# Use not-yet-treated as comparison group, individual data, repeated cross section structure, no clustering, bootstrap
atts_8 <- did::att_gt(yname = "sba_birth", tname = "time", 
  idname = "p_id", gname = "group", data = d_ind, 
  xformla = NULL, weightsname = NULL, est_method = "reg",
  control_group = "notyettreated", bstrap = TRUE, 
  biters = 1000, print_details = FALSE, panel = FALSE,
  cband=FALSE)

# Use not-yet-treated as comparison group, individual data, repeated cross section structure, clustering, bootstrap
atts_9 <- did::att_gt(yname = "sba_birth", tname = "time", 
  idname = "p_id", gname = "group", data = d_ind, 
  xformla = NULL, weightsname = NULL, est_method = "reg",
  control_group = "notyettreated", bstrap = TRUE, 
  biters = 1000, print_details = FALSE, panel = FALSE, 
  clustervars = "dist_id")

df7 <- data.frame(group = atts_7$group, 
  time = atts_7$t, type = "rc", method = "asm", 
  att = atts_7$att, se = atts_7$se)

df8 <- data.frame(group = atts_8$group, 
  time = atts_8$t, type = "rc", method = "boot", 
  att = atts_8$att, se = atts_8$se)

df9 <- data.frame(group = atts_9$group, 
  time = atts_9$t, type = "rc", method = "bootc", 
  att = atts_9$att, se = atts_9$se)

atts_i <- bind_rows(df7,df8,df9)

atts_i %>% 
  select(-se) %>% 
  pivot_wider(names_from=c(type, method), values_from = att) %>%
  mutate_if(is.numeric, round, 3)

atts_i %>% 
  select(-att) %>% 
  pivot_wider(names_from=c(type, method), values_from = se) %>%
  mutate_if(is.numeric, round, 3)



