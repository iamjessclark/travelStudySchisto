#Process data
#Prevalence and intensity calculations

require(tidyverse)
require(leaflet)
require(parzer)
require(sp)

require(sf)
require(tmap)
require(tmaptools)
require(tidycensus)
require(rgdal)
require(broom)
require(rgeos)
require(maptools)
require(viridis)
require(gap)
require(ggmap)
require(spData)
require(spDataLarge)
require(raster)

require(lme4)
require(lmerTest)




# load data ----

all_data <- read.csv("anon_Travel_survey_extra.csv")

#correct the sex mistake
all_data$Sex[which(all_data$Sex == "M (REPEATED)")] <- "M"

# demographics ----
all_data %>% 
  filter(Sex!="")%>%
  mutate(age_group = fct_relevel(age_group, 
                                 "0-5", "6-10", "11-15", "16-20",
                                 "21-25", "26-30", "31-35", "36-40",
                                 "41-45", "46-50", "51-55", "56-60",
                                 "61-65", "66-70", "71-75", "76-80", 
                                 "81-85"))%>%
  filter(!is.na(age_group))%>%
  ggplot()+
  geom_bar(aes(x=age_group, fill=Sex), position="dodge")+
  facet_grid(Testing_loc~.)+
  theme(axis.text.x = element_text(angle = 90))
ggsave("demographics.pdf")  
 
# infection data  ---- #this will introduce NAs where there is writing in the cvs. This is ok
selection <- quo(c(schAd1,schBd1, schAd2,schBd2,schAd3,schBd3))
prevalence <-  all_data %>%
  dplyr::select("ID","CCA_day1", "CCA_day2", "CCA_day3", 
         "schAd1","schBd1", "schAd2","schBd2","schAd3","schBd3", "Testing_loc", "age_group","age_class","Age", "Sex")%>%
  mutate(schAd1=as.numeric(schAd1), schBd1=as.numeric(schBd1),
         schAd2=as.numeric(schAd2), schBd2=as.numeric(schBd2), 
         schAd3=as.numeric(schAd3), schBd3=as.numeric(schBd3))%>%
  filter(Sex!="")%>%
  rowwise()%>%
  mutate(mean_sm= mean(!!selection, na.rm=T))%>%
  ungroup()%>%
  filter(!is.na(mean_sm))%>%
  mutate(infection=case_when(
    mean_sm>0 ~ 1,
    mean_sm==0 & CCA_day1 >= 3 ~ 1,
    mean_sm==0 & CCA_day1 < 3 ~ 0))%>%
  mutate(infectionKK=case_when(
    mean_sm>0 ~ 1,
    mean_sm==0 ~0))%>%
  mutate(infectionCCA=case_when(
    CCA_day1 >= 3 ~ 1,
    CCA_day1 < 3 ~ 0))

write.csv(prevalence, "prevalence.csv")
 
# prevalence by location - infection with KK and CCA

prevalence %>%  
  filter(!is.na(infection)) %>%
  group_by(Testing_loc) %>%  # Group by Testing_loc first to compute total individuals
  mutate(total_n = n()) %>%  # Compute total sample size per Testing_loc
  group_by(infection, Testing_loc, total_n) %>%  # Now group by infection status too
  summarise(n = n(), .groups = "drop") %>%  # Count infected/not infected
  mutate(freq = (n / total_n) * 100) %>%  # Compute true prevalence
  filter(infection == 1) %>%  # Keep only infected individuals
  ggplot() +
  geom_col(aes(y = freq, x = Testing_loc, fill = Testing_loc), position = "dodge", colour = "black") +  # Map fill to Testing_loc
  ylab("Prevalence (%)")+
  theme_bw() +
  scale_fill_manual(values = c("B" = "#663171FF", "Nw" = "#CF3A36FF", 
                               "Nam" = "#EA7428FF", "S" = "#E2998AFF", 
                               "Nan" = "#0C7156FF")) + 
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(size = 20, face = "bold"), 
    legend.position = "none") 

ggsave("PrevLocaleRivka.pdf")

#prev and age group

prevalence %>%  
  filter(!is.na(infection)) %>%
  
  # Compute total sample size per age group & Testing_loc
  group_by(age_group, Testing_loc) %>%
  summarise(total_n = n(), .groups = "drop") %>%  
  
  # Join total sample size back to main dataset
  left_join(
    prevalence %>%
      filter(!is.na(infection)) %>%
      group_by(Testing_loc, age_group, infection) %>%  
      summarise(n = n(), .groups = "drop"),
    by = c("age_group", "Testing_loc")
  ) %>%
  
  # Compute true prevalence
  mutate(freq = (n / total_n) * 100) %>%
  
  # Keep only infected individuals
  filter(infection == 1) %>%
  
  # Plot
  ggplot() +
  geom_col(aes(y = freq, x = age_group, fill = Testing_loc), position = "dodge") +  
  facet_grid(Testing_loc ~ .) +
  theme_bw() +
  ylab("Prevalence (%)") +  # Fixed missing parenthesis
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 16, face = "bold"), 
        legend.position = "none") +
  scale_fill_manual(values = c("B" = "#663171FF", "N" = "#CF3A36FF", 
                               "Nam" = "#EA7428FF", "S" = "#E2998AFF", 
                               "Nan" = "#0C7156FF")) 

ggsave("PrevByAge.pdf")

#prev and age class ####

prevalence %>%  
  filter(!is.na(infection)) %>%
  
  # Compute total sample size per age group & Testing_loc
  group_by(age_class, Testing_loc) %>%
  summarise(total_n = n(), .groups = "drop") %>%  
  
  # Join total sample size back to main dataset
  left_join(
    prevalence %>%
      filter(!is.na(infection)) %>%
      group_by(Testing_loc, age_class, infection) %>%  
      summarise(n = n(), .groups = "drop"),
    by = c("age_class", "Testing_loc")
  ) %>%
  
  # Compute true prevalence
  mutate(freq = (n / total_n) * 100) %>%
  
  # Reorder age_class levels
  mutate(age_class = factor(age_class, levels = c("PSAC", "SAC", "Adult"))) %>%
  
  # Keep only infected individuals
  filter(infection == 1) %>%
  
  # Plot
  ggplot() +
  geom_col(aes(y = freq, x = age_class, fill = Testing_loc), position = "dodge") +  
  facet_grid(Testing_loc ~ .) +
  theme_bw() +
  ylab("Prevalence (%)") +
  theme(
        text = element_text(size = 16, face = "bold"), 
        legend.position = "none") +
  scale_fill_manual(values = c("B" = "#663171FF", "N" = "#CF3A36FF", 
                               "Nam" = "#EA7428FF", "S" = "#E2998AFF", 
                               "Nan" = "#0C7156FF")) 


ggsave("PrevByAgeClass.pdf")





#dataframes of prevalence by community, SAC and each each group

#Community
prevalnceDF <-prevalence %>%  
  filter(!is.na(infection)) %>%
  group_by(Testing_loc) %>%  # Group by Testing_loc first to compute total individuals
  mutate(total_n = n()) %>%  # Compute total sample size per Testing_loc
  group_by(infection, Testing_loc, total_n) %>%  # Now group by infection status too
  summarise(n = n(), .groups = "drop") %>%  # Count infected/not infected
  mutate(freq = (n / total_n) * 100) %>%  # Compute true prevalence
  filter(infection == 1)  # Keep only infected individuals
prevalnceDF


#SAC
prevalnceDFSAC2 <- prevalence %>%  
  filter(!is.na(infection)) %>%
  filter(age_group %in% c("6-10", "11-15")) %>%  # Keep only 6-15 year olds
  
  # Compute total sample size for ALL 6-15 year olds per Testing_loc
  group_by(Testing_loc) %>%
  mutate(total_n = n()) %>%
  
  # Count number of infected individuals per Testing_loc
  group_by(infection, Testing_loc, total_n) %>%  
  summarise(n = n(), .groups = "drop") %>%
  
  # Compute prevalence correctly
  mutate(freq = (n / total_n) * 100) %>%  
  
  # Keep only the "infected" category
  filter(infection == 1)

prevalnceDFSAC2

#Age groups
prevalnceDF_AGE<-prevalence %>%
  group_by(age_group, Testing_loc) %>%
  summarise(total_n = n(), .groups = "drop") %>%  
  
  # Join total sample size back to main dataset
  left_join(
    prevalence %>%
      filter(!is.na(infection)) %>%
      group_by(Testing_loc, age_group, infection) %>%  
      summarise(n = n(), .groups = "drop"),
    by = c("age_group", "Testing_loc")
  ) %>%
  
  # Compute true prevalence
  mutate(freq = (n / total_n) * 100) %>%
  
  # Keep only infected individuals
  filter(infection == 1)

prevalnceDF_AGE

####KK only

#Community
prevalnceDF_KK <-prevalence %>%  
  filter(!is.na(infectionKK)) %>%
  group_by(Testing_loc) %>%  # Group by Testing_loc first to compute total individuals
  mutate(total_n = n()) %>%  # Compute total sample size per Testing_loc
  group_by(infectionKK, Testing_loc, total_n) %>%  # Now group by infection status too
  summarise(n = n(), .groups = "drop") %>%  # Count infected/not infected
  mutate(freq = (n / total_n) * 100) %>%  # Compute true prevalence
  filter(infectionKK == 1)  # Keep only infected individuals
prevalnceDF_KK


#SAC
prevalnceDFSAC_KK <- prevalence %>%  
  filter(!is.na(infectionKK)) %>%
  filter(age_group %in% c("6-10", "11-15")) %>%  # Keep only 6-15 year olds
  
  # Compute total sample size for ALL 6-15 year olds per Testing_loc
  group_by(Testing_loc) %>%
  mutate(total_n = n()) %>%
  
  # Count number of infected individuals per Testing_loc
  group_by(infectionKK, Testing_loc, total_n) %>%  
  summarise(n = n(), .groups = "drop") %>%
  
  # Compute prevalence correctly
  mutate(freq = (n / total_n) * 100) %>%  
  
  # Keep only the "infected" category
  filter(infectionKK == 1)

prevalnceDFSAC_KK


#Sex
prevalnceSex <- prevalence %>%  
  filter(!is.na(infection)) %>%
  
  # Compute total sample size 
  group_by(Sex, Testing_loc) %>%
  mutate(total_n = n()) %>%
  
  # Count number of infected individuals per Sex
  group_by(infection, Sex, total_n, Testing_loc) %>%  
  summarise(n = n(), .groups = "drop") %>%
  
  # Compute prevalence correctly
  mutate(freq = (n / total_n) * 100) %>%  
  
  # Keep only the "infected" category
  filter(infection == 1)

prevalnceSex

# Plot
ggplot(prevalnceSex) +
  geom_col(
    aes(x = Testing_loc, y = freq, fill = Testing_loc),
    position = position_dodge(width = 0.8), colour = "black"
  ) +
  facet_wrap(~ Sex) +
  theme_bw() +
  ylab("Prevalence (%)") +
  xlab("Testing location") +
  theme(
    text = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c(
    "B" = "#663171FF",
    "N" = "#CF3A36FF",
    "Nam" = "#EA7428FF",
    "S" = "#E2998AFF",
    "Nan" = "#0C7156FF"
  ))


#Age groups
prevalnceDF_AGE_KK<-prevalence %>%
  group_by(age_group, Testing_loc) %>%
  summarise(total_n = n(), .groups = "drop") %>%  
  
  # Join total sample size back to main dataset
  left_join(
    prevalence %>%
      filter(!is.na(infectionKK)) %>%
      group_by(Testing_loc, age_group, infectionKK) %>%  
      summarise(n = n(), .groups = "drop"),
    by = c("age_group", "Testing_loc")
  ) %>%
  
  # Compute true prevalence
  mutate(freq = (n / total_n) * 100) %>%
  
  # Keep only infected individuals
  filter(infectionKK == 1)

prevalnceDF_AGE_KK

####CCA only

#Community
prevalnceDF_CCA <-prevalence %>%  
  filter(!is.na(infectionCCA)) %>%
  group_by(Testing_loc) %>%  # Group by Testing_loc first to compute total individuals
  mutate(total_n = n()) %>%  # Compute total sample size per Testing_loc
  group_by(infectionCCA, Testing_loc, total_n) %>%  # Now group by infection status too
  summarise(n = n(), .groups = "drop") %>%  # Count infected/not infected
  mutate(freq = (n / total_n) * 100) %>%  # Compute true prevalence
  filter(infectionCCA == 1)  # Keep only infected individuals

prevalnceDF_CCA

#SAC
prevalnceDFSAC2_CCA <- prevalence %>%  
  filter(!is.na(infectionCCA)) %>%
  filter(age_group %in% c("6-10", "11-15")) %>%  # Keep only 6-15 year olds
  
  # Compute total sample size for ALL 6-15 year olds per Testing_loc
  group_by(Testing_loc) %>%
  mutate(total_n = n()) %>%
  
  # Count number of infected individuals per Testing_loc
  group_by(infectionCCA, Testing_loc, total_n) %>%  
  summarise(n = n(), .groups = "drop") %>%
  
  # Compute prevalence correctly
  mutate(freq = (n / total_n) * 100) %>%  
  
  # Keep only the "infected" category
  filter(infectionCCA == 1)

prevalnceDFSAC2_CCA

#Age groups CCA
prevalnceDF_AGE_CCA<-prevalence %>%
  group_by(age_group, Testing_loc) %>%
  summarise(total_n = n(), .groups = "drop") %>%  
  
  # Join total sample size back to main dataset
  left_join(
    prevalence %>%
      filter(!is.na(infectionCCA)) %>%
      group_by(Testing_loc, age_group, infectionCCA) %>%  
      summarise(n = n(), .groups = "drop"),
    by = c("age_group", "Testing_loc")
  ) %>%
  
  # Compute true prevalence
  mutate(freq = (n / total_n) * 100) %>%
  
  # Keep only infected individuals
  filter(infectionCCA == 1)

prevalnceDF_AGE_CCA

#Rapid assesment results 

rapid <- data.frame(Testing_loc = c(c("B", "N",
                                      "Nam", "S",
                                      "Nan")),
                    Prevalence = c(12, 36, 0, 20,0 ))


ggplot(rapid) +
  geom_col(aes(y = Prevalence, x = Testing_loc, fill = Testing_loc), position = "dodge", colour="black") +
  theme_bw() +
  ylab("Prevalence (%)") + 
  theme(text = element_text(size = 16, face = "bold"), 
        legend.position = "none") +
  scale_fill_manual(values = c("B" = "#663171FF", "N" = "#CF3A36FF", 
                               "Nam" = "#EA7428FF", "S" = "#E2998AFF", 
                               "Nan" = "#0C7156FF")) 


#Comparative KK only resulst in SAC from baseline


KKbl <- data.frame(Testing_loc = c(c("B", "N",
                                      "Nam", "S",
                                      "Nan")),
                    Prevalence = c(10.8, 8.5, 1.88, 19.75,6.66667 ))


ggplot(KKbl) +
  geom_col(aes(y = Prevalence, x = Testing_loc, fill = Testing_loc), position = "dodge", color="black") +
  theme_bw() +
  ylab("Prevalence (%)") + 
  theme(text = element_text(size = 16, face = "bold"), 
        legend.position = "none") +
  scale_fill_manual(values = c("B" = "#663171FF", "N" = "#CF3A36FF", 
                               "Nam" = "#EA7428FF", "S" = "#E2998AFF", 
                               "Nan" = "#0C7156FF")) +
  ylim(0,30)


  
