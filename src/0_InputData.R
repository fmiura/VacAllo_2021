#######################################################################################################
# Input data used for (1)age-SIR simulation and (2)Dutch data simulation.
# 0. package 
# 1. population 
# 2. seroprevalence
# 3. incidence of notified cases
# 4. vaccine efficacy 
# 5. maximum vaccine uptake
# 6. infection hospitalization rate
# 7. infection mortality rate
# 8. contact matrix (Backer 2021 Eurosurveillance)
#######################################################################################################

##0. package -----
library(deSolve)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(hrbrthemes)
library(readxl)
library(readr)
library(styler)

##1. population -----
#nine age groups (for age-SIR simulation)
cbs_excel_9gro <- read_xlsx("CBS_DutchPop2020.xlsx") %>% #https://www.cbs.nl/en-gb/visualisations/dashboard-population/population-pyramid
  mutate(All = Men + Women) %>%
  mutate(Group = case_when(Age <= 10 ~ 1,
                           10 < Age & Age<=20 ~ 2,
                           20 < Age & Age<=30 ~ 3,
                           30 < Age & Age<=40 ~ 4,
                           40 < Age & Age<=50 ~ 5,
                           50 < Age & Age<=60 ~ 6,
                           60 < Age & Age<=70 ~ 7,
                           70 < Age & Age<=80 ~ 8,
                           80 < Age           ~ 9,
                           TRUE ~ 0))
csb_pop_9gro <- cbs_excel_9gro %>% group_by(Group) %>% summarise(sum=sum(All))
n_i <- csb_pop_9gro$sum
pop <- n_i/sum(n_i)

#six age groups (for Dutch data simulation)
cbs_excel_6gro <- read_xlsx("CBS_DutchPop2020.xlsx") %>% #https://www.cbs.nl/en-gb/visualisations/dashboard-population/population-pyramid
  mutate(All = Men + Women) %>%
  mutate(Group = case_when(Age <= 20          ~ 1,
                           20 < Age & Age<=30 ~ 2,
                           30 < Age & Age<=40 ~ 3,
                           40 < Age & Age<=50 ~ 4,
                           50 < Age & Age<=60 ~ 5,
                           60 < Age           ~ 6,
                           TRUE ~ 0))
csb_pop_6gro <- cbs_excel_6gro %>% group_by(Group) %>% summarise(sum=sum(All))

#2.seroprevalence -----
PICO3sero_data <- read_csv("PICO3_seroprev_age_05012020.csv") %>% #https://www.rivm.nl/pienter-corona-studie/resultaten
  mutate(Group = case_when(Age <= 20          ~ 1,
                           20 < Age & Age<=30 ~ 2,
                           30 < Age & Age<=40 ~ 3,
                           40 < Age & Age<=50 ~ 4,
                           50 < Age & Age<=60 ~ 5,
                           60 < Age           ~ 6,
                           TRUE ~ 0))

PICO3_sero <- PICO3sero_data %>% group_by(Group) %>% summarise(avg=mean(seroprevalence))

#3.incidence of notified cases -----
inci_OSRIS <- read.csv("NL_testpos.csv", header = F) #No. of infectious individuals in group i 
inci_OSRIS <- round(inci_OSRIS$V2,0)
inci_OSRIS <- c(sum(inci_OSRIS[1:4]), #20<
         sum(inci_OSRIS[5:6]), #21-30
         sum(inci_OSRIS[7:8]), #31-40
         sum(inci_OSRIS[9:10]),#41-50
         sum(inci_OSRIS[11:12]), #51-60
         sum(inci_OSRIS[13:17]) #60+
)

#4.vaccine efficacy -----
gro <- 9 #no. of age groups
q_i <- c(
  rep(0.948,gro),#Pfizer
  rep(0.941,gro),#Moderena
  rep(0.621,gro),#Astrazeneca
  rep(0.663,gro) #Jansen 
)

#5.maximum vaccine uptake (assumption based on questionnaires in the Netherlands) -----
Uptake_rate <- rep(0.8, gro) #80% for all age groups, "gro" is no. of age groups 

#6.infection hospitalization rate (Walker 2020 Science) -----
h_i <- c(0.10,0.10,0.10,0.20,0.50,1.0,1.6,2.3,2.9,3.9,5.8,7.2,10.2,11.7,14.6,17.7,17.7)*10^(-2)
h_i <- c(mean(h_i[1:2]), #-10y
         mean(h_i[3:4]), #11-20
         mean(h_i[5:6]), #21-30
         mean(h_i[7:8]), #31-40
         mean(h_i[9:10]), #41-50
         mean(h_i[11:12]), #51-60
         mean(h_i[13:14]), #61-70
         mean(h_i[15:16]), #71-80
         mean(h_i[17]) #80+
)

#7.infection mortality rate (O'Driscoll 2020 Nature) -----
d_i <- c(0.0030,0.00069,0.0011,0.0026,0.0079,0.017,0.033,0.055,0.11,0.17,0.30,0.46,0.60,1.5,2.4,4.3,4.3)*10^(-2)
d_i <- c(mean(d_i[1:2]), #-10y
         mean(d_i[3:4]), #11-20
         mean(d_i[5:6]), #21-30
         mean(d_i[7:8]), #31-40
         mean(d_i[9:10]), #41-50
         mean(d_i[11:12]), #51-60
         mean(d_i[13:14]), #61-70
         mean(d_i[15:16]), #71-80
         mean(d_i[17]) #80+
)

#8.contact matrix (Backer 2021 Eurosurveillance) -----
contact10y <- read.table("S2_contact_matrices_withPico3_10y.tsv", sep = '\t', header = TRUE)
a <- contact10y %>% filter(survey=="June 2020", contact_type == "all") 
contactmatrix <-  as.matrix(cbind(a$m_est[1:9], 
                                  a$m_est[10:18],
                                  a$m_est[19:27],
                                  a$m_est[28:36],
                                  a$m_est[37:45],
                                  a$m_est[46:54],
                                  a$m_est[55:63],
                                  a$m_est[64:72],
                                  a$m_est[73:81]))