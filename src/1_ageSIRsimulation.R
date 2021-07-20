#######################################################################################################
# Source codes for age-structured SIR simulations
###0. Package 
###1. Input data 
###2. Age-structured SIR model
###3. Immediate vaccination with age-structured SIR model
###4. Run models
###Reference: https://kinglab.eeb.lsa.umich.edu/480/nls/de.html
#######################################################################################################

###0. Package -----
library(deSolve)
library(dplyr)
library(ggplot2)
library(patchwork)
library(hrbrthemes)
library(readxl)
library(readr)

###1. Input data -----
#Run "0_InputData.R", and below variables are defined there
#n_i 
#pop
#contactmatrix
#h_i
#d_i
Vac40dist_minR <- read_rds("Vac40dist_minR") #Vaccine distribution optimized for min(R), output from "2_Run_VacAllo_9gro1vac.R" 
Vac40dist_minH <- read_rds("Vac40dist_minH") #Vaccine distribution optimized for min(H), output from "2_Run_VacAllo_9gro1vac.R"
Vac40dist_minD <- read_rds("Vac40dist_minD") #Vaccine distribution optimized for min(D), output from "2_Run_VacAllo_9gro1vac.R"

###2. Age-structured SIR model -----
closed.sir.age.model <- function (t, y, parms) {
  ## compartment
  S <- y[grepl("sus", names(y))] #Note: this "S" is the proportion of susceptibles, which is in R code "n_i/sum(n_i)"
  I <- y[grepl("inf", names(y))]
  R <- y[grepl("rec", names(y))]
  H <- y[grepl("hos", names(y))]
  D <- y[grepl("dea", names(y))]
  ## parameter
  beta <- parms$beta
  gamma <- parms$gamma
  pop <- parms$pop
  contactmatrix <- parms$C
  ## ODE
  dSdt <- -beta*rowSums(S*t(I*contactmatrix))
  dIdt <-  beta*rowSums(S*t(I*contactmatrix)) - gamma*I
  dRdt <-                                       gamma*I
  dHdt <- (beta*rowSums(S*t(I*contactmatrix)) - gamma*I)*h_i
  dDdt <- (beta*rowSums(S*t(I*contactmatrix)) - gamma*I)*d_i
  ## combine results into a single vector
  dxdt <- c(dSdt,dIdt,dRdt,dHdt,dDdt)
  ## return result as a list!
  list(dxdt)
}

###3. Immediate vaccination with age-structured SIR model -----
immed.vac.model <- function(timeVac=50, VE=0, Coverage=0){ 
  #timeVac = the timing of depletion of susceptibles and infecteds due to vaccination
  #VE = vaccine efficacy
  #coverage = vaccination coverage 
  
  ##initial value
  xstart = c(sus = pop - (1/n_i), #fully susceptible, except for 1 infected individual in each age group
             inf = (1/n_i), # 1 infected individual in each age group
             rec = rep(0,9),
             hos = (1/n_i)*h_i, 
             dea = (1/n_i)*d_i) 
  parms <- list(beta=1.2*0.2, #R0 * gamma = 1.2 * 0.2
                gamma=0.2, #inverse of infectious period of 5 days
                C = contactmatrix, #No of contacts per stratum per day
                pop=pop) #population (density, Susceptible/Total)
  times <- seq(from=0,to=365,by=1) #time steps, from day 0 to 365
  
  ##trajectory before vaccination
  res <-
    ode(
      func=closed.sir.age.model,
      y=xstart,
      times=times,
      parms=parms
    ) %>%
    as.data.frame()
  
  ##trajectory after vaccination
  res_after <- ode(
    func=closed.sir.age.model,
    y=c(sus=as.numeric(res[timeVac,2:10]*(1-VE)*Coverage + res[timeVac,2:10]*(1-Coverage)),
        inf=as.numeric(res[timeVac,11:19]*(1-VE)*Coverage + res[timeVac,11:19]*(1-Coverage)),
        rec=as.numeric(res[timeVac,20:28]),
        hos=as.numeric(res[timeVac,29:37]*(1-VE)*Coverage + res[timeVac,29:37]*(1-Coverage)),
        dea=as.numeric(res[timeVac,38:46]*(1-VE)*Coverage + res[timeVac,38:46]*(1-Coverage))
        ),
    times=times,
    parms=parms
  ) %>%
    as.data.frame()
  
  ##Store the simulated results
  #combine the results of before/after vaccination
  res_int <- rbind(res[1:(timeVac),],res_after[2:(length(times) - timeVac +1),])
  #store the result as ggplot format
  res_int_data <- as.data.frame(list(
    time=rep(res$time,10),
    sus=c(res_int$sus1,
          res_int$sus2,
          res_int$sus3,
          res_int$sus4,
          res_int$sus5,
          res_int$sus6,
          res_int$sus7,
          res_int$sus8,
          res_int$sus9,
          (res_int$sus1 + res_int$sus2 + res_int$sus3 + res_int$sus4 + res_int$sus5 + res_int$sus6 + res_int$sus7 + res_int$sus8 + res_int$sus9) ),
    inf=c(res_int$inf1,
          res_int$inf2,
          res_int$inf3,
          res_int$inf4,
          res_int$inf5,
          res_int$inf6,
          res_int$inf7,
          res_int$inf8,
          res_int$inf9,
          (res_int$inf1 + res_int$inf2 + res_int$inf3 + res_int$inf4 + res_int$inf5 + res_int$inf6 + res_int$inf7 + res_int$inf8 + res_int$inf9) ),
    rec=c(res_int$rec1,
          res_int$rec2,
          res_int$rec3,
          res_int$rec4,
          res_int$rec5,
          res_int$rec6,
          res_int$rec7,
          res_int$rec8,
          res_int$rec9,
          (res_int$rec1 + res_int$rec2 + res_int$rec3 + res_int$rec4 + res_int$rec5 + res_int$rec6 + res_int$rec7 + res_int$rec8 + res_int$rec9) ),
    hosp=c(res_int$hos1,
           res_int$hos2,
           res_int$hos3,
           res_int$hos4,
           res_int$hos5,
           res_int$hos6,
           res_int$hos7,
           res_int$hos8,
           res_int$hos9,
           (res_int$hos1 + res_int$hos2 + res_int$hos3 + res_int$hos4 + res_int$hos5 + res_int$hos6 + res_int$hos7 + res_int$hos8 + res_int$hos9) ),
    death=c(res_int$dea1,
            res_int$dea2,
            res_int$dea3,
            res_int$dea4,
            res_int$dea5,
            res_int$dea6,
            res_int$dea7,
            res_int$dea8,
            res_int$dea9,
            (res_int$dea1 + res_int$dea2 + res_int$dea3 + res_int$dea4 + res_int$dea5 + res_int$dea6 + res_int$dea7 + res_int$dea8 + res_int$dea9) ),
    age=c(rep("1",length(res_int$inf1)),
          rep("2",length(res_int$inf2)),
          rep("3",length(res_int$inf3)),
          rep("4",length(res_int$inf4)),
          rep("5",length(res_int$inf5)),
          rep("6",length(res_int$inf6)),
          rep("7",length(res_int$inf7)),
          rep("8",length(res_int$inf8)),
          rep("9",length(res_int$inf9)),
          rep("all",length(res_int$inf9)))
  ))
  res_int_data$age <- as.factor(res_int_data$age)
  
  ##Return the simulated results
  return(res_int_data)
}

###4. Run models -----
## Non-optimal vaccination strategies (no vaccination, at random, from young to old, from old to young)
#Without vaccination
nonvac_result <- immed.vac.model(timeVac=50, VE=rep(0,9), Coverage=0)
#Random vaccination
random_result <- immed.vac.model(timeVac=50, VE=rep(0.95,9), Coverage=3.4*pop)
#Young to old 
YoungToOld_result <- immed.vac.model(timeVac=50, VE=rep(0.95,9), Coverage=c(0.8,0.8,0.8,0.8,0.12105726872,0,0,0,0)) #0.12105726872 = (sum(n_i*0.4)-sum(0.8*n_i[1:4]))/n_i[5] 
#Old to young 
OldToYoung_result <- immed.vac.model(timeVac=50, VE=rep(0.95,9), Coverage=c(0,0,0,0,0.67894273128,0.8,0.8,0.8,0.8)) #0.67894273128 = (sum(n_i*0.4)-sum(0.8*n_i[6:9]))/n_i[5] 

##Optimal vaccination strategy
#Store the observed numbers of infected and susceptible individuals from time 15 to 45
nonvac_time45 <- nonvac_result %>% filter(age!="all") %>% filter(time == 45)
nonvac_time15 <- nonvac_result %>% filter(age!="all") %>% filter(time == 15)
incidence_time15to45 <- round((nonvac_time45$rec - nonvac_time15$rec)*sum(n_i),0)
susceptible_time15to45 <- -round((nonvac_time45$sus - nonvac_time15$sus)*sum(n_i),0)
write_rds(incidence_time15to45, "incidence_time15to45")
write_rds(susceptible_time15to45, "susceptible_time15to45")

#Simulate strategies determined by the algorithm
optim_inf_result  <- immed.vac.model(timeVac=50, VE=rep(0.95,9), Coverage=Vac40dist_minR) #for minimizing infection: Vac40dist_minR is the output from "2_Run_VacAllo_9gro1vac.R"
optim_hosp_result  <- immed.vac.model(timeVac=50, VE=rep(0.95,9), Coverage=Vac40dist_minH) #for minimizing hospitalization: Vac40dist_minH is the output from "2_Run_VacAllo_9gro1vac.R"
optim_death_result  <- immed.vac.model(timeVac=50, VE=rep(0.95,9), Coverage=Vac40dist_minD) #for minimizing death: Vac40dist_minD is the output from "2_Run_VacAllo_9gro1vac.R"

#Store results of total incidences ("all") for each strategy
all_nonvac_inf <- nonvac_result %>% filter(age=="all")
all_random_inf <- random_result %>% filter(age=="all")
all_optimR_inf  <- optim_inf_result %>% filter(age=="all")
all_optimH_inf  <- optim_hosp_result %>% filter(age=="all")
all_optimD_inf  <- optim_death_result %>% filter(age=="all")
all_YtoO_inf  <- YoungToOld_result %>% filter(age=="all")
all_OtoY_inf  <- OldToYoung_result %>% filter(age=="all")

summary_inf_sim <- as.data.frame(list(
  time = rep(all_nonvac_inf$time, 5), #5 scenarios
  inf = c(all_nonvac_inf$inf,
          all_random_inf$inf,
          all_YtoO_inf$inf,
          all_OtoY_inf$inf,
          all_optimR_inf$inf),
  hosp = c(all_nonvac_inf$hosp,
          all_random_inf$hosp,
          all_YtoO_inf$hosp,
          all_OtoY_inf$hosp,
          all_optimH_inf$hosp),
  death = c(all_nonvac_inf$death,
          all_random_inf$death,
          all_YtoO_inf$death,
          all_OtoY_inf$death,
          all_optimD_inf$death),
  Strategy =c(rep("No vaccination",length(all_nonvac_inf$time)),
          rep("Random allocation",length(all_nonvac_inf$time)),
          rep("Young-to-old allocation",length(all_nonvac_inf$time)),
          rep("Old-to-young allocation",length(all_nonvac_inf$time)),
          rep("Optimized allocation",length(all_nonvac_inf$time)))
))
summary_inf_sim$Strategy <- as.factor(summary_inf_sim$Strategy)
write_rds(summary_inf_sim, "summary_inf_sim")

#Store results of non-vaccination strategy
test_nonvac <- nonvac_result %>% mutate(FOI = inf/sus)
levels(test_nonvac$age) <- c(levels(factor(contact10y$part_age)),"all")
write_rds(test_nonvac, "test_nonvac")
