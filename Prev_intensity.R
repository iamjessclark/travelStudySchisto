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



#set the function to sort the files out - remove unneeded headers and sort out columns

ParasiteData <- function(object, Testing_loc) {
  # Remove the top rows and unneeded columns
  object <- object[c(-1, -2, -3), c(1:4, 7:9, 11:38)]
  
  # Rename columns
  colnames(object) <- c("ID", "Name", "Age", "Sex", "CCA_day1", "CCA_day2", "CCA_day3",
                        "schAd1", "schBd1", "AscLumAd1", "AscLumBd1", "TTAd1", "TTBd1", 
                        "HWAd1", "HWBd1", "otherd1", "schAd2", "schBd2", "AscLumAd2", "AscLumBd2", 
                        "TTAd2", "TTBd2", "HWAd2", "HWBd2", "otherd2", "schAd3", "schBd3", 
                        "AscLumAd3", "AscLumBd3", "TTAd3", "TTBd3", "HWAd3", "HWBd3", "otherd3", "Testing_loc")
  
  # Ensure ID is a character so it can join
  object$ID <- as.character(object$ID)
  
  # Filter by Testing_loc
  object <- dplyr::filter(object, Testing_loc == !!Testing_loc)
  
  return(object)
}


# load data ----

TravelSurvey <- read.csv("Standard_travel_survey_20221606.csv")




# Process bangaII
bangaII <- read.csv("banga11.csv", na.strings = "NS")
bangaII$Testing_loc <- "BANGA_II"
#bangaII <- bangaII %>% relocate(Testing_loc, .after = 4)


# Process naawa
naawa <- read.csv("naawa_baseline.csv", na.strings = "NS")
naawa$Testing_loc <- "NAAWA"
#naawa <- naawa %>% relocate(Testing_loc, .after = 4)

# Process namukuma
namukuma <- read.csv("namukuma.csv", na.strings = "NS")
namukuma$Testing_loc <- "NAMUKUMA"
#namukuma <- namukuma %>% relocate(Testing_loc, .after = 4)

# Process nanso
nanso <- read.csv("nanso_baseline.csv", na.strings = "NS")
nanso$Testing_loc <- "NANSO_A"
nanso <- nanso %>% dplyr::select(-c(38, 39)) #remove two empty rows
#nanso <- nanso %>% relocate(Testing_loc, .after = 4)

# Process ssese
ssese <- read.csv("ssese_baseline.csv", na.strings = "NS")
ssese$Testing_loc <- "SSESE"
#ssese <- ssese %>% relocate(Testing_loc, .after = 4)


#use this function on our baseline docs

bangaII <- ParasiteData(bangaII, "BANGA_II")
naawa <- ParasiteData(naawa, "NAAWA")
namukuma <- ParasiteData(namukuma, "NAMUKUMA")
nanso <- ParasiteData(nanso, "NANSO_A")
ssese <- ParasiteData(ssese, "SSESE")

#put them all together
base_data <- bind_rows(bangaII, naawa, namukuma, nanso, ssese)


#add age groups and age categories
all_data <- bind_rows(bangaII, naawa, namukuma, nanso, ssese)%>%
  mutate(Age=as.numeric(Age), 
         age_group = case_when(
    Age > 85 ~ "over 85",
    Age >=80 & Age <=85 ~ "81-85",
    Age >=76 & Age <=80 ~ "76-80",
    Age >=71 & Age <=75 ~ "71-75",
    Age >=66 & Age <=70 ~ "66-70",
    Age >=61 & Age <=65 ~ "61-65",
    Age >=56 & Age <=60 ~ "56-60",
    Age >=51 & Age <=55 ~ "51-55",
    Age >=46 & Age <=50 ~ "46-50",
    Age >=41 & Age <=45 ~ "41-45",
    Age >=36 & Age <=40 ~ "36-40",
    Age >=31 & Age <=35 ~ "31-35",
    Age >=26 & Age <=30 ~ "26-30",
    Age >=21 & Age <=25 ~ "21-25",
    Age >=16 & Age <=20 ~ "16-20", 
    Age >=11 & Age <=15 ~ "11-15",
    Age >5  & Age <=10 ~ "06-10",
    Age >0  & Age <=5 ~ "0-5", 
  )
) %>%
  mutate(age_class = case_when(
    Age <7 ~ "PSAC",
    Age >6 & Age <16 ~ "SAC",
    Age >15 ~ "Adult"))


filter <- all_data %>%
    dplyr::select(Age, age_group, age_class)

all_data$Sex[which(all_data$Sex=="M (REPEATED)")] <- "M"

#save.image(file='myEnvironment.RData')

# demographics ----
all_data %>% 
  filter(Sex!="")%>%
  mutate(age_group = fct_relevel(age_group, 
                                 "0-5", "06-10", "11-15", "16-20",
                                 "21-25", "26-30", "31-35", "36-40",
                                 "41-45", "46-50", "51-55", "56-60",
                                 "61-65", "66-70", "71-75", "76-80", 
                                 "81-85", "over 85"))%>%
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
  filter(infection == "infected") %>%  # Keep only infected individuals
  ggplot() +
  geom_col(aes(y = freq, x = Testing_loc, fill = Testing_loc), position = "dodge", colour = "black") +  # Map fill to Testing_loc
  ylab("Prevalence (%)")+
  theme_bw() +
  scale_fill_manual(values = c("BANGA_II" = "#663171FF", "NAAWA" = "#CF3A36FF", 
                               "NAMUKUMA" = "#EA7428FF", "SSESE" = "#E2998AFF", 
                               "NANSO_A" = "#0C7156FF")) + 
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
  scale_fill_manual(values = c("BANGA_II" = "#663171FF", "NAAWA" = "#CF3A36FF", 
                               "NAMUKUMA" = "#EA7428FF", "SSESE" = "#E2998AFF", 
                               "NANSO_A" = "#0C7156FF")) 

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
  scale_fill_manual(values = c("BANGA_II" = "#663171FF", "NAAWA" = "#CF3A36FF", 
                               "NAMUKUMA" = "#EA7428FF", "SSESE" = "#E2998AFF", 
                               "NANSO_A" = "#0C7156FF")) 


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
  filter(infection == "infected")  # Keep only infected individuals

#SAC
prevalnceDFSAC2 <- prevalence %>%  
  filter(!is.na(infection)) %>%
  filter(age_group %in% c("06-10", "11-15")) %>%  # Keep only 6-15 year olds
  
  # Compute total sample size for ALL 6-15 year olds per Testing_loc
  group_by(Testing_loc) %>%
  mutate(total_n = n()) %>%
  
  # Count number of infected individuals per Testing_loc
  group_by(infection, Testing_loc, total_n) %>%  
  summarise(n = n(), .groups = "drop") %>%
  
  # Compute prevalence correctly
  mutate(freq = (n / total_n) * 100) %>%  
  
  # Keep only the "infected" category
  filter(infection == "infected")



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
  filter(infection == "infected")

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

#SAC
prevalnceDFSAC_KK <- prevalence %>%  
  filter(!is.na(infectionKK)) %>%
  filter(age_group %in% c("06-10", "11-15")) %>%  # Keep only 6-15 year olds
  
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

#SAC
prevalnceDFSAC2_CCA <- prevalence %>%  
  filter(!is.na(infectionCCA)) %>%
  filter(age_group %in% c("06-10", "11-15")) %>%  # Keep only 6-15 year olds
  
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


#Rapid assesment results 

rapid <- data.frame(Testing_loc = c(c("BANGA_II", "NAAWA",
                                      "NAMUKUMA", "SSESE",
                                      "NANSO_A")),
                    Prevalence = c(12, 36, 0, 20,0 ))


ggplot(rapid) +
  geom_col(aes(y = Prevalence, x = Testing_loc, fill = Testing_loc), position = "dodge", colour="black") +
  theme_bw() +
  ylab("Prevalence (%)") + 
  theme(text = element_text(size = 16, face = "bold"), 
        legend.position = "none") +
  scale_fill_manual(values = c("BANGA_II" = "#663171FF", "NAAWA" = "#CF3A36FF", 
                               "NAMUKUMA" = "#EA7428FF", "SSESE" = "#E2998AFF", 
                               "NANSO_A" = "#0C7156FF")) 


#Comparative KK only resulst in SAC from baseline


KKbl <- data.frame(Testing_loc = c(c("BANGA_II", "NAAWA",
                                      "NAMUKUMA", "SSESE",
                                      "NANSO_A")),
                    Prevalence = c(10.8, 8.5, 1.88, 19.75,6.66667 ))


ggplot(KKbl) +
  geom_col(aes(y = Prevalence, x = Testing_loc, fill = Testing_loc), position = "dodge", color="black") +
  theme_bw() +
  ylab("Prevalence (%)") + 
  theme(text = element_text(size = 16, face = "bold"), 
        legend.position = "none") +
  scale_fill_manual(values = c("BANGA_II" = "#663171FF", "NAAWA" = "#CF3A36FF", 
                               "NAMUKUMA" = "#EA7428FF", "SSESE" = "#E2998AFF", 
                               "NANSO_A" = "#0C7156FF")) +
  ylim(0,30)


  
