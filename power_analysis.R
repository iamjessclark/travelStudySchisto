# pilot data #

require(tidyverse)

ra_buikwe <- read.csv("RA_buikwe.csv")

standard_q <- read.csv("Standard_Travel_Survey.csv")
standard_q$child_id <- as.factor(substr(standard_q$Unique.identifier.for.participant,1,5))
standard_q <- standard_q %>%
  dplyr::select(child_id, age, sex, location, travel_next_3m)

nyenga <- read.csv("Travel_Survey_NYENGA.csv")
nyenga$child_id <- as.factor(substr(nyenga$Unique.identifier.for.participant,1,5))
nyenga <- nyenga %>%
  dplyr::select(child_id, age, sex, location, travel_next_3m)

senyi <- read.csv("Travel_Survey_SENYI.csv")
senyi$child_id <- substr(senyi$Unique.identifier.for.participant,1,5)
senyi[1,ncol(senyi)] <- "B1529"
senyi$child_id <- as.factor(senyi$child_id)
senyi <- senyi %>%
  dplyr::select(child_id, age, sex, location, travel_next_3m)

ra_buikwe <- ra_buikwe %>%
  mutate(infected = ifelse(mean>0, 1, 0), 
         child_id = as.factor(child_id), 
         location = as.factor(location),
         infected = as.factor(infected))
  
ra_buikwe %>%
  filter(!is.na(mean))%>%
  group_by(location, infected) %>%
  summarise(n = n())%>%
  mutate(freq = (n / sum(n))*100)%>%
  complete(infected, fill=list(infected=1))%>%
  filter(infected!=0)%>%
  ggplot()+
  geom_col(aes(x=infected, y=freq))+
  facet_grid(.~location)
ggsave("prevalence_RA_buikwe.pdf")

pa_df <- bind_rows(standard_q, nyenga, senyi)

pa_df %>%
  left_join(ra_buikwe, by=c("child_id", "age", "sex", "location"))
  filter(!is.na(mean))%>%
  group_by(location, infected) %>%
  summarise(n = n())%>%
  mutate(freq = (n / sum(n))*100)%>%
  complete(infected, fill=list(infected=1))%>%
  filter(infected!=0)%>%
  ggplot()+
  geom_col(aes(x=infected, y=freq))+
  facet_grid(.~location)

pa_df %>% 
  ggplot()+
  geom_point(aes(x=location, y=travel_next_3m))
  
pa_df %>%
  group_by(location) %>%
  summarise(var=var(travel_next_3m))

# no difference between locations in infection intensity 
TukeyHSD(aov(mean ~ location, data=pa_df))

# travel distance between location
summary(glm.nb(travel_next_3m ~ location, data=pa_df))

pa_df$mean <- round(pa_df$mean)
glmer.nb(mean~1+(1|location), data=pa_df) # this gives me town level variance in epg

# difference in infection prevalence 
prev <- glm(infected ~ location, data=pa_df, family=binomial(link='logit'))
anova(prev, test="Chisq")
summary(prev)

slope <- 1
pop_mean_eggs <- 0.89 # mean egg count
pop_sd_eggs <- 6.28 # sd egg count

pop_mean_journeys <- 9.74 # mean journeys
pop_sd_journeys <- 35.10 # sd of number of journeys

# use a log normal distribution to generate number of journey counts. 
# this can be rounded to make integers
location <- log((pop_mean_journeys^2)/sqrt(pop_sd_journeys^2+pop_mean_journeys^2))
shape <- sqrt(log(1+(pop_sd_journeys^2/pop_mean_journeys^2)))
 
n_locations <- 5 
reps <- 200
p_list <- vector()
town.sd <- 3.636
runs <- 1000

for(experimentloop in 1:runs){ 
  RE <- (towns = rep(LETTERS[1:n_locations], each = reps))
  children <- 1:reps
  Rand_Eff_Vals <- (town_eff = rnorm(n_locations, 0, 3.636)) # this is the size of the effect not the epg
  Rand_Eff_Vals <- (town_eff = rep(town_eff, each = reps)) # the effect of location is the same across all children in that location
  Travel_Eff <- round(rlnorm(reps*n_locations, location, shape))
  residuals <- rnorm(reps*n_locations,0,pop_sd_eggs)
  dat <- data.frame(children, Travel_Eff, as.factor(RE), Rand_Eff_Vals, residuals)
  dat$sim_epg <- with(dat, log(pop_mean_eggs+(slope*Travel_Eff)+Rand_Eff_Vals)+residuals)
  dat$sim_epg[which(is.nan(dat$sim_epg)==T)] <- 0
  dat$sim_epg[which(dat$sim_epg<0)] <- 0
  dat$sim_epg <- round(dat$sim_epg)                    
  dat$scale_travel <- scale(dat$Travel_Eff)                    
  #journeys <- vector()
  #eggs <- vector()
  # this is based on "true egg count"
  #for(r in 1:reps){
    #thisx <- round(rlnorm(1, location, shape))
    #journeys <- append(journeys, thisx)
    #predicted <- (intercept + (slope*thisx))
    #residual <- rnorm(1,0,stdd)
    #eggs <- append(eggs, (log(predicted)+residual)) # logged so it is suitable for the NB analysis
  #}
  #eggs[which(eggs<0)] <- 0 # this isn't perfect but not sure how else to do this
  #eggs <- round(eggs,digits=0) # rounding to make the values integers
  #dataset <- data.frame(journeys, eggs)
  #analysis <- anova(glm.nb(eggs~journeys, data = dataset))
  analysis <- summary(glmmTMB(sim_epg ~ (1|RE) + scale_travel, data=dat,family = nbinom2))
  p <- analysis$coefficients$cond[2,4]
  p_list <- append(p_list, p)
}  
  
significant <- length(which(p_list<0.05))
print(paste("power=",(significant/runs)))

# 125 people x 8 locations = 91% power to detect an effect size of a single egg per gram increase with each unit of travel (based on number of upcoming journeys)
# 20% drop out = 83.7% power to detect the same effect size 

# Banga II
# Bulutwe B
# Kirukwe-Gaba
# Kiwugi
# Mawoloba
# Nambeta
# Kavule
# BBogo
# Namukuma 

#### Trying with fewer locations #### 

town_choice <- ra_buikwe %>%
  filter(location=="BANGA II" | location=="NAAWA CELL1" | location=="SSESE" |
           location=="NANSO A"| location=="KAVULE" | location=="KIDOKOLO")%>%
  mutate(mean=round(mean))

town_choice%>%
  summarise(mean_epg=mean(mean, na.rm=T), 
            std_dev=sd(mean, na.rm=T))


# in just the above location
pop_mean_epg <- 1.59
std_dev <- 7.1

# taken from the wider survey - not all the above locations had surveys done
pop_mean_journeys <- 9.74 # mean journeys
pop_sd_journeys <- 35.10 # sd of number of journeys

# use a log normal distribution to generate number of journey counts. 
# this can be rounded to make integers
location <- log((pop_mean_journeys^2)/sqrt(pop_sd_journeys^2+pop_mean_journeys^2))
shape <- sqrt(log(1+(pop_sd_journeys^2/pop_mean_journeys^2)))
slope <- 1
n_locations <- 6 
reps <- 200
p_list <- vector()
runs <- 1000

summary(glmmTMB(mean ~ (1|location), data=town_choice,family = nbinom2))

town.sd <- 1.857

for(experimentloop in 1:runs){ 
  RE <- (towns = rep(LETTERS[1:n_locations], each = reps))
  children <- 1:reps
  Rand_Eff_Vals <- (town_eff = rnorm(n_locations, 0, town.sd)) # this is the size of the effect not the epg, 
  Rand_Eff_Vals <- (town_eff = rep(town_eff, each = reps)) # the effect of location is the same across all children in that location
  Travel_Eff <- round(rlnorm(reps*n_locations, location, shape))
  residuals <- rnorm(reps*n_locations,0,std_dev)
  dat <- data.frame(children, Travel_Eff, as.factor(RE), Rand_Eff_Vals, residuals)
  dat$sim_epg <- with(dat, log(pop_mean_epg+(slope*Travel_Eff)+Rand_Eff_Vals)+residuals)
  dat$sim_epg[which(is.nan(dat$sim_epg)==T)] <- 0
  dat$sim_epg[which(dat$sim_epg<0)] <- 0
  dat$sim_epg <- round(dat$sim_epg)                    
  dat$scale_travel <- scale(dat$Travel_Eff)                    
  #journeys <- vector()
  #eggs <- vector()
  # this is based on "true egg count"
  #for(r in 1:reps){
  #thisx <- round(rlnorm(1, location, shape))
  #journeys <- append(journeys, thisx)
  #predicted <- (intercept + (slope*thisx))
  #residual <- rnorm(1,0,stdd)
  #eggs <- append(eggs, (log(predicted)+residual)) # logged so it is suitable for the NB analysis
  #}
  #eggs[which(eggs<0)] <- 0 # this isn't perfect but not sure how else to do this
  #eggs <- round(eggs,digits=0) # rounding to make the values integers
  #dataset <- data.frame(journeys, eggs)
  #analysis <- anova(glm.nb(eggs~journeys, data = dataset))
  analysis <- summary(glmmTMB(sim_epg ~ (1|RE) + scale_travel, data=dat,family = nbinom2))
  p <- analysis$coefficients$cond[2,4]
  p_list <- append(p_list, p)
}  

significant <- length(which(p_list<0.05))
print(paste("power=",(significant/runs)))