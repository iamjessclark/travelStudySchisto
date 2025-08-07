#Making the DAG
#

rm(list = ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/HPGH/_PhD/Jess_data/data_manipulation")

library(dagitty)
library(ggdag)
library(brms)
library(dplyr)
library(ggplot2)

#install data
#
data <-read.csv("spreadsheets/Travel_prevEDIT.csv")
nrow(data) #586


#What about travel is a categorical variable

data$travel_frequency <- as.factor(data$travel_frequency)


# Define DAG

dag7 <- dagitty('
dag {
bb="0,0,1,1"
Act [pos="0.722,0.563"]
Age [pos="0.165,0.174"]
Dur [pos="0.587,0.795"]
Inf [outcome,pos="0.354,0.845"]
Location [pos="0.492,0.321"]
MDA [pos="0.305,0.550"]
Sex [pos="0.745,0.165"]
Travel [exposure,pos="0.466,0.578"]
Act -> Dur
Act -> Inf
Age -> Act
Age -> Dur
Age -> Inf
Age -> Location
Age -> MDA
Age -> Travel
Dur -> Inf
Location -> Act
Location -> Dur
Location -> Inf
Location -> MDA
Location -> Travel
MDA -> Inf
Sex -> Act
Sex -> Dur
Sex -> Inf
Sex -> Location
Sex -> Travel
Travel -> Act
Travel -> Inf
Travel -> MDA
}

')

plot(dag7)

# Define the DAG (assuming you have already defined `dag` as in the previous response)
implied_independencies <- impliedConditionalIndependencies(dag7)

print(implied_independencies)

# Ac_N _||_ MDA | FCON, a__N, lc_N
# DrMn _||_ FCON | Ac_N, Sx_x, a__N, lc_N
# DrMn _||_ MDA | FCON, a__N, lc_N
# DrMn _||_ MDA | Ac_N, Sx_x, a__N, lc_N
# MDA _||_ Sx_x | FCON, a__N, lc_N
# Sx_x _||_ a__N
# 

#use bayesian and glms to test independencies

data$ActNEW <- as.factor(data$ActNEW)
data$ActNEW <- relevel(data$ActNEW, ref = "None") 
data$Location <- as.factor(data$Location)
data$Location <- relevel(data$Location, ref = "NAAWA") 
data$infection <- as.factor(data$infection)
data$infection <- relevel(data$infection, ref = "0") 


# Act _||_ MDA | Age, Freq, Lctn

model1 <- brm((MDA~ActNEW+age_class+travel_frequency+Location),
              data=data,
              family = bernoulli(link = "logit"),
              iter = 10000, warmup = 3000, chains = 4, cores = 8,
              control = list(max_treedepth = 25))

summary(model1) #Yes independent

mod1 <- glm(MDA~ActNEW+age_class+travel_frequency+Location,
            data=data, family="binomial")

summary(mod1) #yes independent


# Age _||_ Sex
model3 <- brm(Sex_bin~age_class,
              data=data,
              family = bernoulli(),
              iter = 10000, warmup = 3000, chains = 4, cores = 8,
              control = list(max_treedepth = 15))

mod3 <- glm(Sex_bin ~ age_class, data=data, family="binomial")

summary(mod3) #

summary(model3) #No but that is because we had more PSAC boys than girls and this is an exogenous variable so ok


# DurL _||_ Freq | Act, Age, Sex, location

model4<-brm(travel_frequency ~ DurMin+ActNEW + age_class + Sex_bin+Location, 
            data = data, 
            family = categorical(),
            iter = 20000, warmup = 6000, chains = 4, cores = 8,
            control = list(max_treedepth = 25))

summary(model4) #yes independent

library(nnet)
mod4<-multinom(travel_frequency ~ DurMin+ActNEW + age_class + Sex_bin+Location, 
               data = data)

mod4_sum<-summary(mod4) 

# Extract coefficients and standard errors
coefs <- mod4_sum$coefficients
ses <- mod4_sum$standard.errors

# Calculate z-values and p-values manually
z_values <- coefs / ses
p_values <- 2 * (1 - pnorm(abs(z_values)))  # Two-tailed

# Combine all into a labelled data frame
results_list <- lapply(1:nrow(coefs), function(i) {
  data.frame(
    Outcome_Level = rownames(coefs)[i],
    Variable = colnames(coefs),
    Estimate = coefs[i, ],
    Std.Error = ses[i, ],
    z.value = z_values[i, ],
    p.value = p_values[i, ],
    row.names = NULL
  )
})

results_df <- do.call(rbind, results_list)

print(results_df) #yes independent

# DurL _||_ MDA | Age, Freq, Lctn

#have a look at the raw data
table(data$DurL, data$MDA)

model5 <-brm(MDA~DurMin+age_class+travel_frequncy+Location,
             data=data,
             family = bernoulli(),
             iter = 10000, warmup = 3000, chains = 4, cores = 8,
             control = list(max_treedepth = 15))

summary(model5) #yes independent

mod5 <-glm(MDA~DurMin+age_class+travel_frequency+Location,
           data=data,
           family = binomial)

summary(mod5) #yes independent



# DrMn _||_ MDA | Ac_N, Sx_x, a__N, lc_N

# 
model6 <-brm(MDA~DurMin+ActNEW+age_class+Location+Sex_bin,
             data=data,
             family = bernoulli(),
             iter = 10000, warmup = 3000, chains = 4, cores = 8,
             control = list(max_treedepth = 15))

summary(model6) #Yes 

mod6 <-glm(MDA~DurMin+ActNEW+age_class+Location+Sex_bin,
           data=data,
           family = binomial)

summary(mod6) #Yes 


# MDA _||_ Sex | Age, Freq, Lctn
model7 <-brm(Sex_bin~MDA+age_class+travel_frequency+Location,
             data=data,
             family = bernoulli(),
             iter = 10000, warmup = 3000, chains = 4, cores = 8,
             control = list(max_treedepth = 15))

summary(model7) #Yes 

# MDA _||_ Sex | Age, Freq, Lctn
mod7 <-glm(Sex_bin~MDA+age_class+travel_frequency+Location,
           data=data,
           family = binomial)

summary(mod7) #Yes 


#Create SCM #########
#SCM

# Age = U_Age
# Sex = U_Sex
# 
# Location = f0(Sex) + +f0(Age)+U_Location  
# Act = f1(Age, Sex, Location) + U_Act  
# DurL = f2(Act, Age, Sex, Location) + U_DurL  
# Freq = f3(Act, Age, Sex, Location) + U_Freq  
# MDA = f4(Age, Freq, Location) + U_MDA  
# Inf = f5(Act, DurL, Freq, MDA, Age, Sex, Location) + U_Inf  
# 
# Assuming linearity
# 
# Age = U_Age
# Sex = U_Sex
# 
# Location = α1 * Sex + α2 * Age + U_Location
# 
# Act = β1 * Age + β2 * Sex + β3 * Location + U_Act
# 
# DurL = γ1 * Act + γ2 * Age + γ3 * Sex + γ4 * Location + U_DurL
# 
# Freq = δ1 * Act + δ2 * Age + δ3 * Sex + δ4 * Location + U_Freq
# 
# MDA = θ1 * Age + θ2 * Freq + θ3 * Location + U_MDA
# 
# Inf = φ1 * Act + φ2 * DurL + φ3 * Freq + φ4 * MDA + φ5 * Age + φ6 * Sex + φ7 * Location + U_Inf


###############################################################################
####### Adjustment Sets #######################################################
####### ######################################################################
#===========================================================================#
#                     Edges                                                 #
#===========================================================================#

#===========================================================================#
#               Act -> DurL #####                                           #
#===========================================================================# 

## Total causal effect of Activity on Duration
adjustmentSets(dag7, "Act", "Dur", effect="total") #{age and sex and location}

## Direct causal effect of Activity on Duration
adjustmentSets(dag7, "Act", "Dur", effect="direct") #{age and sex and location}
#=============================================================#
#                 Act -> Freq #####                           #
#=============================================================#
### Total causal effect of Activity on frequency
adjustmentSets(dag7, "Act", "Travel", effect="total") # {age and sex and location}

## Direct causal effect of Activity on frequency
adjustmentSets(dag7, "Act", "Travel", effect="direct") # {age and sex and location}


#==================================================================#
#                     Act -> Inf                                #####
#==================================================================#
# ## Total causal effect of Activity on Infection
adjustmentSets(dag7, "Act", "Inf", effect="total") # {age and sex and location}

## Direct causal effect of Activity on Infection
adjustmentSets(dag7, "Act", "Inf",  effect="direct") # {age, duration, frequency, sex and location}


#==========================================================================
#                       Age -> Act                                    #####
#===========================================================================
# ## Total causal effect of Age on Activity
adjustmentSets(dag7, "Age", "Act", effect="total") #{} none

## Direct causal effect of Age and Activity
adjustmentSets(dag7, "Age", "Act", effect="direct") #{Sex and location}


#=====================================================================#
#                 Age -> DurL                                     #####
# ===================================================================#
# # ## Total causal effect of Age on Duration
adjustmentSets(dag7, "Age", "Dur", effect="total") #{} none

## Direct causal effect of Age on Duration
adjustmentSets(dag7, "Age", "Dur", effect="direct") #{Act, sex, location} 

#=====================================================================#
#           Age -> Freq                                           ######
# ====================================================================#
# # # ## Total causal effect of Age on frequency
adjustmentSets(dag7, "Age", "Travel", effect="total") #{} none

## Direct causal effect of Age on frequency
adjustmentSets(dag7, "Age", "Travel", effect="direct") #{Act, sex, location} 

#==========================================================================
# Age -> Inf #####
#========================================================================= 
# # # ## Total causal effect of Age on infection
adjustmentSets(dag7, "Age", "Inf", effect="total") #{} none

## Direct causal effect of Age on Infection
adjustmentSets(dag7, "Age", "Inf", effect="direct") #{Act, DurL, Freq, MDA, Sex and Location} 

#=====================================================================
# Age -> MDA ####
# ===================================================================
# # # ## Total causal effect of Age on MDA
adjustmentSets(dag7, "Age", "MDA", effect="total") #{} none

## Direct causal effect of age on MDA
adjustmentSets(dag7, "Age", "MDA", effect="direct") #{Freq and location} 

#=====================================================================
# Age -> Location ####
# ===================================================================
# # # ## Total causal effect of Age on MDA
adjustmentSets(dag7, "Age", "Location", effect="total") #{} none

## Direct causal effect of age on MDA
adjustmentSets(dag7, "Age", "Location", effect="direct") #{} 

#=====================#
# DurL -> Inf ####
# #=====================#
# # # ## Total causal effect of Duration on infection
adjustmentSets(dag7, "Dur", "Inf", effect="total") #{Act, Age, Sex, location} 

## Direct causal effect of Duration on infection
adjustmentSets(dag7, "Dur", "Inf", effect="direct") #{Act, Age, sex, location} 

#=====================#
# Freq -> Inf####
# #=====================#
# # # ## Total causal effect of frequency on Infection
adjustmentSets(dag7, "Travel", "Inf", effect="total") #{Age, Sex, location} 

## Direct causal effect of frequency on Infection
adjustmentSets(dag7, "Travel", "Inf",  effect="direct") #{Act, Activity, Age, sex, location, MDA} 

#=====================#
# Freq -> MDA####
# #=====================#
# # # ## Total causal effect of frequency on Infection
adjustmentSets(dag7, "Travel", "MDA", effect="total") #{Age, location} 

## Direct causal effect of frequency on Infection
adjustmentSets(dag7, "Travel", "MDA", effect="direct") #{Age,location} 

#=====================#
# MDA -> Inf####
# #=====================#
# # ## Total causal effect of of MDA on Infection
adjustmentSets(dag7, "MDA", "Inf", effect="total") #{Freq, Age, location} 

## Direct causal effect of MDA on Infection
adjustmentSets(dag7, "MDA", "Inf", effect="direct") #{Freq, Age, location} 

#=====================#
# Sex -> Act####
# #=====================#
# # ## Total causal effect of of Sex on activity
adjustmentSets(dag7, "Sex", "Act", effect="total") #{} none

## Direct causal effect of Sex on Activity
adjustmentSets(dag7, "Sex", "Act", effect="direct") #{Age, location} none

#=====================#
# Sex -> DurL ####
# #=====================#
# 
# # # ## Total causal effect of of Sex on Duration
adjustmentSets(dag7, "Sex", "Dur", effect="total") #{} none

## Direct causal effect of Sex on Duration
adjustmentSets(dag7, "Sex", "Dur", effect="direct") #{Act, Age, Location} 

#=====================#
# Sex -> Freq ####
# #=====================#
# # # ## Total causal effect of of Sex on Frequency 
adjustmentSets(dag7, "Sex", "Travel", effect="total") #{} none

## Direct causal effect of Sex on Frequency 
adjustmentSets(dag7, "Sex", "Travel", effect="direct") #{Act, Age, Location} 

#=====================#
# Sex -> Inf ####
# #=====================#

# # ## Total causal effect of of Sex on Infection
adjustmentSets(dag7, "Sex", "Inf", effect="total") #{} none

## Direct causal effect of Sex on Infection
adjustmentSets(dag7, "Sex", "Inf", effect="direct") #{Act, Age, DurL, Freq, Location} 

#=====================#
# Sex -> Location ####
# #=====================#

# # ## Total causal effect of of Sex on Infection
adjustmentSets(dag7, "Sex", "Location", effect="total") #{} none

## Direct causal effect of Sex on Infection
adjustmentSets(dag7, "Sex", "Location", effect="direct") #{} 

#=================================================================
# loc_Num ->Act_Num
# ===============================================================

# # ## Total causal effect of of Sex on Infection
adjustmentSets(dag7, "Location", "Act", effect="total") #{Sex and age} 

## Direct causal effect of Sex on Infection
adjustmentSets(dag7, "Location", "Act", effect="direct") #{Sex and Age} 

#========================================================================
# loc_Num ->DurMin
# =======================================================================

# # ## Total causal effect of of Sex on Infection
adjustmentSets(dag7, "Location", "Dur", effect="total") #{Sex and age} 

## Direct causal effect of Sex on Infection
adjustmentSets(dag7, "Location", "Dur", effect="direct") #{Act, Sex and Age}

#===========================================================================
# loc_Num ->MDA
# =========================================================================
# 
# # ## Total causal effect of of Sex on Infection
adjustmentSets(dag7, "Location", "MDA", effect="total") #{Sex and age} 

## Direct causal effect of Sex on Infection
adjustmentSets(dag7, "Location", "MDA", effect="direct") #{Freq and Age}

#===========================================================================
# loc_Num ->infection
# ==========================================================================
# 
# # ## Total causal effect of of Sex on Infection
adjustmentSets(dag7, "Location", "Inf", effect="total") #{Sex and age} 

## Direct causal effect of Sex on Infection
adjustmentSets(dag7, "Location", "Inf", effect="direct") #{Act, Dur, Freq, MDA, Sex and Age}

#==========================================================================
# loc_Num ->FreqCONT
# =========================================================================
# 
# # ## Total causal effect of of Sex on Infection
adjustmentSets(dag7, "Location", "Travel", effect="total") #{Sex and age} 

## Direct causal effect of Sex on Infection
adjustmentSets(dag7, "Location", "Travel", effect="direct") #{Act, Sex and Age}



#============================================================================
#   Models of direct and total effects to parameterise the SCM   
# ===========================================================================


#===========================================================================#
#               Act -> DurL #####                                           #
#===========================================================================#

priors <- get_custom_prior("negbinomial")

#total and direct the same

Act_dur <- brm(DurMin ~ ActNEW + age_class+Sex_bin+Location,
               data=data,
               family="negbinomial",
               prior=priors,
               iter = 10000, warmup = 3000, chains = 4, cores = 8,
               control = list(max_treedepth = 15))

#model fit
pp_check(Act_dur)  # shows dens_overlay plot by default
pp_check(Act_dur, type = "error_hist", ndraws = 11)
pp_check(Act_dur, type = "scatter_avg", ndraws = 100)
pp_check(Act_dur, type = "stat_2d") #slight underestimation but probs best I can get
pp_check(Act_dur, type = "rootogram")
pp_check(Act_dur, type = "loo_pit_overlay")


summary(Act_dur) #all significant and all positive

#make a table of results
# Get summary object
sum_fit_dur <- summary(Act_dur)

# Extract and convert the main coefficients
coef_table_dur <- as.data.frame(sum_fit_dur$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table_dur <- as.data.frame(sum_fit_dur$spec_pars)   # for 'shape', etc.
dist_table_dur <- as.data.frame(sum_fit_dur$distributional)  # for 'hu' predictors

# Combine all
full_table_dur <- rbind(coef_table_dur, hu_table_dur, dist_table_dur)
# 
# ###########################################################################
#=============================================================#
#                 Act -> Freq #####                           #
#=============================================================#


#Total and direct the same
#Age, sex and loaction

Act_freq <- brm(travel_frequency ~ ActNEW + age_class+Sex_bin+Location,
                data=data,
                family="categorical",
                iter = 10000, warmup = 3000, chains = 4, cores = 8,
                control = list(max_treedepth = 15))

#model fit
pp_check(Act_freq )  # shows dens_overlay plot by default
pp_check(Act_freq , type = "error_hist", ndraws = 11)
pp_check(Act_freq , type = "scatter_avg", ndraws = 100)
pp_check(Act_freq , type = "stat_2d") #
pp_check(Act_freq , type = "rootogram")
pp_check(Act_freq , type = "loo_pit_overlay")


summary(Act_freq ) #



#==================================================================#
#                     Act -> Inf                                #####
#==================================================================#

#This has its own separate code for model as going to make graphs of it for

#==========================================================================
#                       Age -> Act                                    #####
#===========================================================================

#DIRECT
Act_age <- brm(ActNEW ~age_class+Sex_bin+Location,
               data=data,
               family="categorical",
               iter = 20000, warmup = 4000, chains = 4, cores = 8,
               control = list(max_treedepth = 25))

#model fit
pp_check(Act_age)  # shows dens_overlay plot by default
pp_check(Act_age, type = "error_hist", ndraws = 11)
pp_check(Act_age, type = "stat_2d") 

summary(Act_age) #

#make a table of results
# Get summary object
sum_fit_act <- summary(Act_age)

# Extract and convert the main coefficients
coef_table_act <- as.data.frame(sum_fit_act$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table_act <- as.data.frame(sum_fit_act$spec_pars)   # for 'shape', etc.
dist_table_act <- as.data.frame(sum_fit_act$distributional)  # for 'hu' predictors

# Combine all
full_table_act <- rbind(coef_table_act, hu_table_act, dist_table_act)

full_table_act


#TOTAL
Act_age_tot <- brm(ActNEW ~age_class,
                   data=data,
                   family="categorical",
                   iter = 20000, warmup = 4000, chains = 4, cores = 8,
                   control = list(max_treedepth = 25))

#model fit
pp_check(Act_age_tot)  # shows dens_overlay plot by default
pp_check(Act_age_tot, type = "error_hist", ndraws = 11)
pp_check(Act_age_tot, type = "stat_2d") 

summary(Act_age_tot) #


# 
# #=====================================================================#
#                 Age -> DurL                                     #####
# ===================================================================#

priors <- get_custom_prior("negbinomial")

#DIRECT
Age_dur <- brm(DurMin ~age_class+Sex_bin+Location+ActNEW,
               data=data,
               prior=priors,
               family="negbinomial",
               iter = 10000, warmup = 3000, chains = 4, cores = 8,
               control = list(max_treedepth = 15))

#model fit
pp_check(Age_dur)  # shows dens_overlay plot by default
pp_check(Age_dur, type = "error_hist", ndraws = 11)
pp_check(Age_dur, type = "stat_2d") #underestimation but probs ok

summary(Age_dur) #not significant effects

#make a table of results
# Get summary object
sum_fit_act_ag <- summary(Act_dur)

# Extract and convert the main coefficients
coef_table_ag <- as.data.frame(sum_fit_act_ag$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table_ag<- as.data.frame(sum_fit_act_ag$spec_pars)   # for 'shape', etc.
dist_table_ag <- as.data.frame(sum_fit_act_ag$distributional)  # for 'hu' predictors

# Combine all
full_table_ag <- rbind(coef_table_ag, hu_table_ag, dist_table_ag)

#TOTAL
Age_dur_tot <- brm(DurMin ~age_class,
                   data=data,
                   prior=priors,
                   family="negbinomial",
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 15))

#model fit
pp_check(Age_dur_tot)  # shows dens_overlay plot by default
pp_check(Age_dur_tot, type = "error_hist", ndraws = 11)
pp_check(Age_dur_tot, type = "stat_2d") #underestimation but probs ok

summary(Age_dur_tot) #not significant effects

#make a table of results
# Get summary object
sum_fit_act_t <- summary(Age_dur_tot)

# Extract and convert the main coefficients
coef_table_t <- as.data.frame(sum_fit_act_t$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table_t<- as.data.frame(sum_fit_act_t$spec_pars)   # for 'shape', etc.
dist_table_t<- as.data.frame(sum_fit_act_t$distributional)  # for 'hu' predictors

# Combine all
full_table_ag <- rbind(coef_table_t, hu_table_t, dist_table_t)


# #=====================================================================#
#           Age -> Freq                                           ######
# ====================================================================#


Age_freq <- brm(travel_frequency~age_classNEW+Sex_bin+Location,
                data=data,
                family="categorical",
                iter = 10000, warmup = 3000, chains = 4, cores = 8,
                control = list(max_treedepth = 15))

#model fit
pp_check(Age_freq)  # shows dens_overlay plot by default
pp_check(Age_freq, type = "error_hist", ndraws = 11)
pp_check(Age_freq, type = "stat_2d") 

summary(Age_freq) #not significant effects


# Travel Freq	      PSAC	                SAC	              YoungAdult
# Once	          -1.05 (↓ less likely)	+0.62 (slightly ↑)	+1.08 (↑ more likely)
# 3×	            -0.60 (↓ slightly less)	-0.20	            +0.93
# 6×	            -2.56 (↓↓ much less)	-15.71 (!!)	        +0.13
# 13×	            +0.47	                  +1.17 (↑)	        +0.51
# 26×	            -2.31 (↓↓)	          +0.44	              +0.60
# 92×	            -1.36 (↓↓)	          +0.52	              +0.41
# 
# #==========================================================================
# Age -> Inf #####
#========================================================================= 
priors <- get_custom_prior("bernoulli")

#DIRECT####
#Bernoulli first

Age_inf_dir_bern <- brm(infected_bin ~age_class+ActNEW+DurMin+travel_frequency+MDA+Sex_bin+Location,
                        data=data,
                        family=bernoulli(),
                        prior=priors,
                        iter = 10000, warmup = 3000, chains = 4, cores = 8,
                        control = list(max_treedepth = 15))


#model fit
pp_check(Age_inf_dir_bern )  # shows dens_overlay plot by default
pp_check(Age_inf_dir_bern , type = "error_hist", ndraws = 11)
pp_check(Age_inf_dir_bern , type = "stat_2d") 

summary(Age_inf_dir_bern) #


#gamma of egg postive only

filt <- data %>%
  mean_sm >0

Age_inf_dir <- brm(mean_sm ~age_class+ActNEW+DurMin+travel_frequency+MDA+Sex_bin+Location,
                   data=filt,
                   family=gamma(),
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 15))


#model fit
pp_check(Age_inf_dir)  # shows dens_overlay plot by default
pp_check(Age_inf_dir, type = "error_hist", ndraws = 11)
pp_check(Age_inf_dir, type = "stat_2d") 

summary(Age_inf_dir) #not significant effects

#make a table of results
# Get summary object
sum_fit_age <- summary(Age_inf_dir)

# Extract and convert the main coefficients
coef_table_age <- as.data.frame(sum_fit_age$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table_age<- as.data.frame(sum_fit_age$spec_pars)   # for 'shape', etc.
dist_table_age<- as.data.frame(sum_fit_age$distributional)  # for 'hu' predictors

# Combine all
full_table_age <- rbind(coef_table_age, hu_table_age, dist_table_age)



#TOTAL####
#bernoulli
#
Age_inf_tot_bern <- brm(infected_bin ~age_class,
                        data=data,
                        family=bernoulli(),
                        prior=priors,
                        iter = 10000, warmup = 3000, chains = 4, cores = 8,
                        control = list(max_treedepth = 15))

Age_inf_tot <- brm(mean_sm ~age_class,
                   data=filt,
                   family=gamma(),
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 15))


#model fit
pp_check(Age_inf_tot)  # shows dens_overlay plot by default
pp_check(Age_inf_tot, type = "error_hist", ndraws = 11)
pp_check(Age_inf_tot, type = "stat_2d") 

summary(Age_inf_tot)


#make a table of results
# Get summary object
sum_fit <- summary(Age_inf_tot)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table <- rbind(coef_table, hu_table, dist_table)

#=====================================================================
# Age -> MDA ####
# ===================================================================
priors <- get_custom_prior("bernoulli")
#DIRECT 
Age_MDA <- brm(MDA ~age_class+travel_frequency+Location,
               data=data,
               family="bernoulli",
               prior=priors,
               iter = 10000, warmup = 3000, chains = 4, cores = 8,
               control = list(max_treedepth = 25))

#model fit
pp_check(Age_MDA)  # shows dens_overlay plot by default
pp_check(Age_MDA, type = "error_hist", ndraws = 11)
pp_check(Age_MDA, type = "stat_2d") 

summary(Age_MDA) #significant effects

#make a table of results
# Get summary object
sum_fit <- summary(Age_MDA)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table <- rbind(coef_table, hu_table, dist_table)


#TOTAL
Age_MDA_tot <- brm(MDA ~age_class,
                   data=data,
                   prior=priors,
                   family="bernoulli",
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 25))

#model fit
pp_check(Age_MDA_tot)  # shows dens_overlay plot by default
pp_check(Age_MDA_tot, type = "error_hist", ndraws = 11)
pp_check(Age_MDA_tot, type = "stat_2d") 

summary(Age_MDA_tot) #significant effects

#make a table of results
# Get summary object
sum_fit <- summary(Age_MDA_tot)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table <- rbind(coef_table, hu_table, dist_table)
#=====================================================================
#       Age -> Location ####
# ===================================================================
# TOTAL and DIRECT
Age_loc <- brm(Location ~age_class,
               data=data,
               family="categorical",
               iter = 10000, warmup = 3000, chains = 4, cores = 8,
               control = list(max_treedepth = 15))

#model fit
pp_check(Age_loc)  # shows dens_overlay plot by default
pp_check(Age_loc, type = "error_hist", ndraws = 11)
pp_check(Age_loc, type = "stat_2d") 

summary(Age_loc) #significant effects 

#make a table of results
# Get summary object
sum_fit <- summary(Age_loc)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table <- rbind(coef_table, hu_table, dist_table)
# 
# #===================================================================================#
#               DurL -> Inf ####
# #===================================================================================#



#This has its own code as I made graphs of it too
# 



#======================================================================================#
#             Freq -> Inf                                                           ####
# #====================================================================================#

#This has its own code as I made graphs of it too
# 
# #===============================================================================================#
#               Freq -> MDA                                                                     ####
# #===============================================================================================#


priors <- get_custom_prior("bernoulli")

#DIRECT and TOTAL
Freq_MDA <- brm(MDA ~ travel_frequency+age_class+Location,  
                data=data,
                family="bernoulli",
                prior=priors,
                iter = 10000, warmup = 3000, chains = 4, cores = 8,
                control = list(max_treedepth = 25))

#model fit
pp_check(Freq_MDA)  # shows dens_overlay plot by default
pp_check(Freq_MDA, type = "error_hist", ndraws = 11)
pp_check(Freq_MDA, type = "stat_2d") 

summary(Freq_MDA) #no significant

#make a table of results
# Get summary object
sum_fit <- summary(Freq_MDA)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table <- rbind(coef_table, hu_table, dist_table)
#======================================================================================#
#               MDA -> Inf                                                          ####
# #====================================================================================#
#This has its own code as I made graphs of it too


##=====================#
# Sex -> Act####
# #=====================#

#DIRECT
Sex_Act <- brm(ActNEW ~ Sex_bin+age_class+Location,
               data=data,
               family="categorical",
               iter = 10000, warmup = 3000, chains = 4, cores = 8,
               control = list(max_treedepth = 15))

#model fit
pp_check(Sex_Act)  # shows dens_overlay plot by default
pp_check(Sex_Act, type = "error_hist", ndraws = 11)
pp_check(Sex_Act, type = "stat_2d") 

summary(Sex_Act)


#make a table of results
# Get summary object
sum_fit_dir <- summary(Sex_Act)

# Extract and convert the main coefficients
coef_table_dir <- as.data.frame(sum_fit_dir$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table_dir<- as.data.frame(sum_fit_dir$spec_pars)   # for 'shape', etc.
dist_table_dir<- as.data.frame(sum_fit_dir$distributional)  # for 'hu' predictors

# Combine all
full_table_dir <- rbind(coef_table_dir, hu_table_dir, dist_table_dir)
view(full_table_dir)

#TOTAL


Sex_Act_tot <- brm(ActNEW ~ Sex_bin,
                   data=data,
                   family="categorical",
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 15))

#model fit
pp_check(Sex_Act_tot)  # shows dens_overlay plot by default
pp_check(Sex_Act_tot, type = "error_hist", ndraws = 11)
pp_check(Sex_Act_tot, type = "stat_2d") 

summary(Sex_Act_tot) 

#make a table of results
# Get summary object
sum_fit_tot <- summary(Sex_Act_tot)

# Extract and convert the main coefficients
coef_table_tot <- as.data.frame(sum_fit_tot$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table_tot<- as.data.frame(sum_fit_tot$spec_pars)   # for 'shape', etc.
dist_table_tot<- as.data.frame(sum_fit_tot$distributional)  # for 'hu' predictors

# Combine all
full_table_tot <- rbind(coef_table_tot, hu_table_tot, dist_table_tot)
#=====================#
# Sex -> DurL ####
# #=====================#
#
priors <- get_custom_prior("negbinomial")  

#Direct
Sex_Dur <- brm(DurMin ~ Sex_bin+ActNEW+age_class + Location,
               data=data,
               family="negbinomial",
               prior=priors,
               iter = 10000, warmup = 3000, chains = 4, cores = 8,
               control = list(max_treedepth = 15))

#model fit
pp_check(Sex_Dur)  # shows dens_overlay plot by default
pp_check(Sex_Dur, type = "error_hist", ndraws = 11)
pp_check(Sex_Dur, type = "stat_2d")  #slight over estimation but ok

summary(Sex_Dur) 

#make a table of results
# Get summary object
sum_fit<- summary(Sex_Dur)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table <- rbind(coef_table, hu_table, dist_table)

#TOTAL
Sex_Dur_tot <- brm(DurMin ~ Sex_bin,
                   data=data,
                   family="negbinomial",
                   prior=priors,
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 15))

#model fit
pp_check(Sex_Dur_tot)  # shows dens_overlay plot by default
pp_check(Sex_Dur_tot, type = "error_hist", ndraws = 11)
pp_check(Sex_Dur_tot, type = "stat_2d")  #slight over estimation but ok

summary(Sex_Dur_tot) #not signifiant

#make a table of results
# Get summary object
sum_fit_tot <- summary(Sex_Dur_tot)

# Extract and convert the main coefficients
coef_table_tot <- as.data.frame(sum_fit_tot$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table_tot<- as.data.frame(sum_fit_tot$spec_pars)   # for 'shape', etc.
dist_table_tot<- as.data.frame(sum_fit_tot$distributional)  # for 'hu' predictors

# Combine all
full_table_tot <- rbind(coef_table_tot, hu_table_tot, dist_table_tot)
#==========================================================================#
#                       Sex -> Freq ####
# #========================================================================#


Sex_Freq <- brm(travel_frequency ~ Sex_bin+ActNEW+age_classNEW + Location,
                data=data,
                family="categorical",
                iter = 20000, warmup = 3000, chains = 4, cores = 8,
                control = list(max_treedepth = 25))

#model fit
pp_check(Sex_Freq)  # shows dens_overlay plot by default
pp_check(Sex_Freq, type = "error_hist", ndraws = 11)
pp_check(Sex_Freq, type = "stat_2d")  #slight over estimation but ok

summary(Sex_Freq) #

#===========================================================================#
#           Sex -> Inf                                                    ####
# #=========================================================================#

priors <_get_custom_priors("bernoulli")
#DIRECT 
#bernoulli first
#
Sex_Inf_dir_bern <- brm(infected_bin ~ Sex_bin+ActNEW + age_class + DurMin + travel_frequency + Location,
                        data=data,
                        prior=priors,
                        family="bernoulli",
                        iter = 10000, warmup = 3000, chains = 4, cores = 8,
                        control = list(max_treedepth = 15))

#model fit
pp_check(Sex_Inf_dir_bern )  # shows dens_overlay plot by default
pp_check(Sex_Inf_dir_bern , type = "error_hist", ndraws = 11)
pp_check(Sex_Inf_dir_bern , type = "stat_2d")  #slight over estimation but ok

summary(Sex_Inf_dir_bern) #

#now postive only
Sex_Inf_dir <- brm(mean_sm ~ Sex_bin+ActNEW + age_class + DurMin + travel_frequency + Location,
                   data=filt,
                   family="gamma",
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 15))

#model fit
pp_check(Sex_Inf_dir )  # shows dens_overlay plot by default
pp_check(Sex_Inf_dir , type = "error_hist", ndraws = 11)
pp_check(Sex_Inf_dir , type = "stat_2d")  #slight over estimation but ok

summary(Sex_Inf_dir) #

#make a table of results
# Get summary object
sum_fit<- summary(Sex_Inf_dir)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table<- rbind(coef_table, hu_table, dist_table)

#TOTAL
#bernoulli first

Sex_Inf_tot_bern<- brm(infected_bin ~ Sex_bin,
                       data=data,
                       prior=priors,
                       family="bernoulli",
                       iter = 10000, warmup = 3000, chains = 4, cores = 8,
                       control = list(max_treedepth = 15))

#model fit
pp_check(Sex_Inf_tot_bern )  # shows dens_overlay plot by default
pp_check(Sex_Inf_tot_bern , type = "error_hist", ndraws = 11)
pp_check(Sex_Inf_tot_bern , type = "stat_2d")  #slight over estimation but ok


summary(Sex_Inf_tot) #

#positive only
Sex_Inf_tot<- brm(mean_sm ~ Sex_bin,
                  data=filt,
                  family="gamma",
                  iter = 10000, warmup = 3000, chains = 4, cores = 8,
                  control = list(max_treedepth = 15))

#model fit
pp_check(Sex_Inf_tot )  # shows dens_overlay plot by default
pp_check(Sex_Inf_tot , type = "error_hist", ndraws = 11)
pp_check(Sex_Inf_tot , type = "stat_2d")  #slight over estimation but ok

summary(Sex_Inf_tot) #

#make a table of results
# Get summary object
sum_fit<- summary(Sex_Inf_tot)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table<- rbind(coef_table, hu_table, dist_table)







#===========================================================================#
#                 Sex -> Location ####
# #========================================================================## 

#DIRECT and TOTAL

Sex_loc <- brm(Location ~ Sex_bin,, 
               data=data,
               family="categorical",
               iter = 10000, warmup = 3000, chains = 4, cores = 8,
               control = list(max_treedepth = 15))

#model fit
pp_check(Sex_loc)  # shows dens_overlay plot by default
pp_check(Sex_loc, type = "error_hist", ndraws = 11)
pp_check(Sex_loc, type = "stat_2d")  #perfect

summary(Sex_loc) #

#make a table of results
# Get summary object
sum_fit<- summary(Sex_loc)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table<- rbind(coef_table, hu_table, dist_table)

#=======================================================================================================#
#             Location ->Act_Num
# ======================================================================================================# 

#TOTAL AND DIRECT

Act_loc <- brm(ActNEW ~ Location +Sex_bin + age_class,
               data=data,
               family="categorical",
               iter = 20000, warmup = 3000, chains = 4, cores = 8,
               control = list(max_treedepth = 25))

#model fit
pp_check(Act_loc)  # shows dens_overlay plot by default
pp_check(Act_loc, type = "error_hist", ndraws = 11)
pp_check(Act_loc, type = "stat_2d")  #good
summary(Act_loc) #

#make a table of results
# Get summary object
sum_fit<- summary(Act_loc)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table<- rbind(coef_table, hu_table, dist_table)

#========================================================================
#           Location ->DurMin
# =======================================================================# 

priors <- get_custom_prior("negbinomial") 

#DIRECT
Dur_loc <- brm(DurMin ~ Location +ActNEW+ Sex_bin + age_class, 
               data=data,
               family="negbinomial",
               prior=priors,
               iter = 10000, warmup = 3000, chains = 4, cores = 8,
               control = list(max_treedepth = 15))

#model fit
pp_check(Dur_loc)  # shows dens_overlay plot by default
pp_check(Dur_loc, type = "error_hist", ndraws = 11)
pp_check(Dur_loc, type = "stat_2d")  #slight over estimation but ok

summary(Dur_loc) 

#make a table of results
# Get summary object
sum_fit<- summary(Dur_loc)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table<- rbind(coef_table, hu_table, dist_table)

#TOTAL
Dur_loc_tot <- brm(DurMin ~ Location +Sex_bin + age_class, 
                   data=data,
                   family="negbinomial",
                   prior=priors,
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 15))

#model fit
pp_check(Dur_loc_tot)  # shows dens_overlay plot by default
pp_check(Dur_loc_tot, type = "error_hist", ndraws = 11)
pp_check(Dur_loc_tot, type = "stat_2d")  #slight over estimation but ok

summary(Dur_loc_tot) 

#make a table of results
# Get summary object
sum_fit<- summary(Dur_loc_tot)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table<- rbind(coef_table, hu_table, dist_table)
#===========================================================================
#         Location ->MDA
# =========================================================================# 
# 
priors <- get_custom_prior("bernoulli")


# DIRECT
Loc_MDA <- brm(MDA ~ Location +travel_frequency+age_class,
               data=data,
               family="bernoulli",
               prior=priors,
               iter = 10000, warmup = 3000, chains = 4, cores = 8,
               control = list(max_treedepth = 25))

#model fit
pp_check(Loc_MDA)  # shows dens_overlay plot by default
pp_check(Loc_MDA, type = "error_hist", ndraws = 11)
pp_check(Loc_MDA, type = "stat_2d")  #slight over estimation but ok

summary(Loc_MDA) 

#make a table of results
# Get summary object
sum_fit<- summary(Loc_MDA)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table<- rbind(coef_table, hu_table, dist_table)

# TOTAL
Loc_MDA_tot <- brm(MDA ~ Location +Sex_bin+age_class,
                   data=data,
                   family="bernoulli",
                   prior=priors,
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 25))

#model fit
pp_check(Loc_MDA_tot)  # shows dens_overlay plot by default
pp_check(Loc_MDA_tot, type = "error_hist", ndraws = 11)
pp_check(Loc_MDA_tot, type = "stat_2d")  #slight over estimation but ok

summary(Loc_MDA_tot) 

#make a table of results
# Get summary object
sum_fit<- summary(Loc_MDA_tot)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table<- rbind(coef_table, hu_table, dist_table)



#===========================================================================
#           Location ->infection
# ==========================================================================# 
# 

#DIRECT bernoulli first 
#
loc_inf_bern <- brm(infected_bin ~ Location +ActNEW + DurMin + travel_frequency + MDA +Sex_bin + age_class,
                    data=data,
                    family="bernoulli",
                    prior=priors,
                    iter = 10000, warmup = 3000, chains = 4, cores = 8,
                    control = list(max_treedepth = 15))

#model fit
pp_check(loc_inf_bern)  # shows dens_overlay plot by default
pp_check(loc_inf_bern, type = "error_hist", ndraws = 11)
pp_check(loc_inf_bern, type = "stat_2d")  

summary(loc_inf_bern) 

#now positive only
loc_inf <- brm(mean_sm ~ Location +ActNEW + DurMin + travel_frequency + MDA +Sex_bin + age_class,
               data=filt,
               family="gamma",
               iter = 10000, warmup = 3000, chains = 4, cores = 8,
               control = list(max_treedepth = 15))

#model fit
pp_check(loc_inf)  # shows dens_overlay plot by default
pp_check(loc_inf, type = "error_hist", ndraws = 11)
pp_check(loc_inf, type = "stat_2d")  

summary(loc_inf) #

#make a table of results
# Get summary object
sum_fit<- summary(loc_inf)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table<- rbind(coef_table, hu_table, dist_table)


#TOTAL
#bern first
#
loc_inf_tot_bern <- brm(infected_bin ~ Location +Sex_bin + age_class,
                        data=data,
                        family="bernoulli",
                        prior=priors,
                        iter = 10000, warmup = 3000, chains = 4, cores = 8,
                        control = list(max_treedepth = 15))

#model fit
pp_check(loc_inf_tot_bern)  # shows dens_overlay plot by default
pp_check(loc_inf_tot_bern, type = "error_hist", ndraws = 11)
pp_check(loc_inf_tot_bern, type = "stat_2d")  

summary(loc_inf_tot_bern) #

#positive only
loc_inf_tot <- brm(mean_sm ~ Location +Sex_bin + age_class,
                   data=filt,
                   family="gamma",
                   prior=priors,
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 15))

#model fit
pp_check(loc_inf_tot)  # shows dens_overlay plot by default
pp_check(loc_inf_tot, type = "error_hist", ndraws = 11)
pp_check(loc_inf_tot, type = "stat_2d")  

summary(loc_inf_tot) #

#make a table of results
# Get summary object
sum_fit<- summary(loc_inf_tot)

# Extract and convert the main coefficients
coef_table <- as.data.frame(sum_fit$fixed)

# Optional: If you have distributional parameters (like `hu`, `shape`), add them too
hu_table<- as.data.frame(sum_fit$spec_pars)   # for 'shape', etc.
dist_table<- as.data.frame(sum_fit$distributional)  # for 'hu' predictors

# Combine all
full_table<- rbind(coef_table, hu_table, dist_table)




#==========================================================================
#         Location ->FreqCONT
# =========================================================================# 
loc_freq <- brm(travel_frequency~ Location +ActNEW + Sex_bin + age_classNEW,
                data=data,
                family="categorical",
                iter = 10000, warmup = 3000, chains = 4, cores = 8,
                control = list(max_treedepth = 15))

#model fit
pp_check(loc_freq)  # shows dens_overlay plot by default
pp_check(loc_freq, type = "error_hist", ndraws = 11)
pp_check(loc_freq, type = "stat_2d")  #slight over estimation but ok

summary(loc_freq) #