#######################################################################################################
# Source codes for the allocation algorithm (6 age groups, 4 vaccine type)
###0. Run "0_InputData.R" 
###1. Define function to optimize vaccine allocations for different objectives
###2. Run simulations for minimizing R/H/D
###3. Store the results for visualization
#######################################################################################################

#####0. Run "0_InputData.R" ------

#####1. Define function to optimize vaccine allocations for different objectives ------
VacAllo_func_4vac  <- function(z_plot_length=100, target="R"){
#####1. input data#####
gro <- 6   #No. of age groups
VacType <- 4 #No. of vaccine types (Pfizer, Moderena, AstraZeneca, Jansen)

#Total population
pop <- csb_pop_6gro$sum #No. of total individuals in group i 
n_i <- rep(pop, VacType)

#Incidence
x_i <- inci_OSRIS #Incidence of notified cases in group i  
x_i <- rep(x_i, VacType)

#Seroprevalence
dat_sero <- PICO3_sero$avg #seroprevalence from PICO3
s_i <- round(pop*(1-dat_sero/100),0) #No. of susceptible individuals in group i 
s_i <- rep(s_i, VacType)

#per contact probability of acquiring infection for group i
a_i <- rep(1,gro) 
a_i <- rep(a_i, VacType)

#per contact infectiousness of group i 
c_i <- rep(1,gro)
c_i <- rep(c_i, VacType)

###Vaccine stock 
VacA_pro <- 0.46 #Pfizer 46%
VacB_pro <- 0.08 #Moderna 8%
VacC_pro <- 0.22 #AstraZeneca 22%
VacD_pro <- 0.24 #Jansen 24%

###Maximum uptake (willingness to vaccination)
Uptake_rate_value <- 0.8 #80% is the max uptake
Uptake_rate <- rep(0.8, gro) #80% for all age groups, "gro" age-groups 

###vaccine efficacy 
q_i <- c(
  rep(0.948,gro),#Pfizer
  rep(0.941,gro),#Moderena
  rep(0.621,gro),#Astrazeneca
  rep(0.663,gro) #Jansen 
) 

###Infection hospitalization rate (Meehan 2020, Table S3)
h_i <- c(0.10,0.10,0.10,0.20,0.50,1.0,1.6,2.3,2.9,3.9,5.8,7.2,10.2,11.7,14.6,17.7,17.7)*10^(-2)
h_i <- c(mean(h_i[1:4]), #20<
         mean(h_i[5:6]), #21-30
         mean(h_i[7:8]), #31-40
         mean(h_i[9:10]),#41-50
         mean(h_i[11:12]), #51-60
         mean(h_i[13:17]) #60+
)
h_i <- rep(h_i,VacType)

###Infection mortality rate (Meehan 2020, Table S3)
d_i <- c(0.0030,0.00069,0.0011,0.0026,0.0079,0.017,0.033,0.055,0.11,0.17,0.30,0.46,0.60,1.5,2.4,4.3,4.3)*10^(-2)
d_i <- c(mean(d_i[1:4]), #20<
              mean(d_i[5:6]), #21-30
              mean(d_i[7:8]), #31-40
              mean(d_i[9:10]),#41-50
              mean(d_i[11:12]), #51-60
              mean(d_i[13:17]) #60+
)
d_i <- rep(d_i,VacType)

#####2. Importance weight#####
### Define the approximated next generation matrix 
#top right and left eigenvectors
w_1 <- (1 / sum(x_i[1:gro])) * x_i[1:gro]
v_1 <- (sum(x_i[1:gro]) / sum(x_i[1:gro]^2/s_i[1:gro])) * (x_i[1:gro] / s_i[1:gro])
#approximated K
R_ini <- 1.2
approxK <- R_ini * w_1 %x% t(v_1)

### Target indexes
#initial values
R_ini <- 1.2
H_ini <- Re(eigen(diag(h_i[1:gro]) %*% approxK)$values[1])
D_ini <- Re(eigen(diag(d_i[1:gro]) %*% approxK)$values[1])

R <- R_ini
H <- H_ini
D <- D_ini
#normalization factor h 
h <- 1/sum((c_i[1:gro]/a_i[1:gro])*(x_i[1:gro]^2/s_i[1:gro])) 

#reproduction number 
dR_du_i <- R * (-h) * q_i * (c_i / a_i) * (x_i / s_i) * (x_i / n_i)
#hospitalization number
dH_du_i <- dR_du_i * h_i
#death number
dD_du_i <- dR_du_i * d_i

#####3. Run loop#####
#Minimize R0
u_i <- rep(0,gro*VacType)

z <- sum(pop*Uptake_rate)#Uptake_rate <- rep(0.8, gro) #80% for all age groups, "gro" age-groups 
z <- round(z, 0)
z_vacA <- round(z*VacA_pro, 0)  #Vaccine stock of A
z_vacB <- round(z*VacB_pro, 0)  #Vaccine stock of B
z_vacC <- round(z*VacC_pro, 0)  #Vaccine stock of C
z_vacD <- round(z*VacD_pro, 0)  #Vaccine stock of D

y_R_ini <-  R * (h) * q_i * (c_i / a_i) * (x_i / s_i) * (x_i / n_i)
y_H_ini <-  R * (h) * q_i * (c_i / a_i) * (x_i / s_i) * (x_i / n_i) * h_i
y_D_ini <-  R * (h) * q_i * (c_i / a_i) * (x_i / s_i) * (x_i / n_i) * d_i

dy_R_du <- rep(0, gro * VacType)
dy_H_du <- rep(0, gro * VacType)
dy_D_du <- rep(0, gro * VacType)

dy_du <- rep(0, gro * VacType) ### target index ###

y_loop <- matrix(0, z, gro * VacType) ### target index ###
u_loop <- matrix(0, z, gro * VacType)
R_loop <- rep(0,z)
H_loop <- rep(0,z)
D_loop <- rep(0,z)

if(target=="R"){
  y <- y_R_ini
}else if(target=="H"){
  y <- y_H_ini
}else if(target=="D"){
  y <- y_D_ini
}else {
  print("no target index is specified")
}

if(target=="R"){
  for (l in 1:z) {
    y_loop[l,] <- y #vector, gro*VacType
    u_loop[l,] <- u_i #vector, gro*VacType
    R_loop[l] <- R #scalar, 1 variable
    H_loop[l] <- H #scalar, 1 variable
    D_loop[l] <- D #scalar, 1 variable
    
    #check remaining vaccine stocks
    if(sum(u_i[1:gro])>=z_vacA){y[1:gro] <- 0}
    if(sum(u_i[(gro+1):(gro*2)])>=z_vacB){y[(gro+1):(gro*2)] <- 0}
    if(sum(u_i[(gro*2+1):(gro*3)])>=z_vacC){y[(gro*2+1):(gro*3)] <- 0}
    if(sum(u_i[(gro*3+1):(gro*4)])>=z_vacD){y[(gro*3+1):(gro*4)] <- 0}
    
    for (s in 1:(gro*VacType)) {
      ### choose the set of same age-group with different vaccines  (17 age groups)
      if(s == gro | s == gro*2 | s == gro*3 | s == gro*4){
        s_set <- c(gro,gro*2,gro*3,gro*4)
      }else {s_set <- c( (s%%gro), (s%%gro)+gro, (s%%gro)+gro*2, (s%%gro)+gro*3) }
      ### check the maximum vaccine uptake per age group. "sum(u_[s_set])" is equal to vacA + vacB + vacC + vacD in age-group "s".
      if(sum(u_i[s_set]) >= (n_i[s]*Uptake_rate_value)){y[s_set] <- 0} #to check if u_i[s] reaches the maximum vaccine uptake. if yes, y[s] <- 0
      if(max(y)==y[s]){
        ###
        u_i[s] <- u_i[s] + 1
        ###
        dy_R_du[s_set] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) ,VacType)
        dy_H_du[s_set] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) * h_i[s], VacType)
        dy_D_du[s_set] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) * d_i[s], VacType)
        ###
        dy_du[s_set] <- dy_R_du[s_set] #change the importance weights of group s for all vaccine types (and therefore "s_set" is used here)
        ###
        dR_du <- dR_du_i[s] #only s (not "s_set" here). This is a scalar (not vector)
        dH_du <- dH_du_i[s] #only s (not "s_set" here). This is a scalar (not vector)
        dD_du <- dD_du_i[s] #only s (not "s_set" here). This is a scalar (not vector)
      }
    }
    ###
    y <- y + dy_du
    ###
    R <- R + dR_du
    H <- H + dH_du
    D <- D + dD_du
    
    dy_R_du <- rep(0,gro*VacType)
    dy_H_du <- rep(0,gro*VacType)
    dy_D_du <- rep(0,gro*VacType)
    ###
    dy_du <- rep(0,gro*VacType)
    ###
    dR_du <- 0
    dH_du <- 0
    dD_du <- 0
  }
}
else if(target=="H"){
  for (l in 1:z) {
    y_loop[l,] <- y
    u_loop[l,] <- u_i
    R_loop[l] <- R
    H_loop[l] <- H
    D_loop[l] <- D
    
    #check remaining vaccine stocks
    if(sum(u_i[1:gro])>=z_vacA){y[1:gro] <- 0}
    if(sum(u_i[(gro+1):(gro*2)])>=z_vacB){y[(gro+1):(gro*2)] <- 0}
    if(sum(u_i[(gro*2+1):(gro*3)])>=z_vacC){y[(gro*2+1):(gro*3)] <- 0}
    if(sum(u_i[(gro*3+1):(gro*4)])>=z_vacD){y[(gro*3+1):(gro*4)] <- 0}
    
    for (s in 1:(gro*VacType)) {
      ### choose the set of same age-group with different vaccines  (17 age groups)
      if(s == gro | s == gro*2 | s == gro*3 | s == gro*4){
        s_set <- c(gro,gro*2,gro*3,gro*4)
      }else {s_set <- c( (s%%gro), (s%%gro)+gro, (s%%gro)+gro*2, (s%%gro)+gro*3) }
      ###
      if(sum(u_i[s_set]) >= (n_i[s]*Uptake_rate_value)){y[s_set] <- 0} #to check if u_i[s] reaches the maximum vaccine uptake. if yes, y[s] <- 0
      if(max(y)==y[s]){
        
        u_i[s] <- u_i[s] + 1
        ###
        dy_R_du[s_set] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) ,VacType)
        dy_H_du[s_set] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) * h_i[s], VacType)
        dy_D_du[s_set] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) * d_i[s], VacType)
        ###
        dy_du[s_set] <- dy_H_du[s_set] #change the importance weights of group s for all vaccine types (and therefore "s_set" is used here)
        ###
        dR_du <- dR_du_i[s]
        dH_du <- dH_du_i[s]
        dD_du <- dD_du_i[s]
      }
    }
    y <- y + dy_du
    R <- R + dR_du
    H <- H + dH_du
    D <- D + dD_du
    
    dy_R_du <- rep(0,gro*VacType)
    dy_H_du <- rep(0,gro*VacType)
    dy_D_du <- rep(0,gro*VacType)
    ###
    dy_du <- rep(0,gro*VacType)
    ###
    dR_du <- 0
    dH_du <- 0
    dD_du <- 0
  }
}
else if(target=="D"){
  for (l in 1:z) {
    y_loop[l,] <- y
    u_loop[l,] <- u_i
    R_loop[l] <- R
    H_loop[l] <- H
    D_loop[l] <- D
    
    #check remaining vaccine stocks
    if(sum(u_i[1:gro])>=z_vacA){y[1:gro] <- 0}
    if(sum(u_i[(gro+1):(gro*2)])>=z_vacB){y[(gro+1):(gro*2)] <- 0}
    if(sum(u_i[(gro*2+1):(gro*3)])>=z_vacC){y[(gro*2+1):(gro*3)] <- 0}
    if(sum(u_i[(gro*3+1):(gro*4)])>=z_vacD){y[(gro*3+1):(gro*4)] <- 0}
    
    for (s in 1:(gro*VacType)) {
      ### choose the set of same age-group with different vaccines  (17 age groups)
      if(s == gro | s == gro*2 | s == gro*3 | s == gro*4){
        s_set <- c(gro,gro*2,gro*3,gro*4)
      }else {s_set <- c( (s%%gro), (s%%gro)+gro, (s%%gro)+gro*2, (s%%gro)+gro*3) }
      ###
      if(sum(u_i[s_set]) >= (n_i[s]*Uptake_rate_value)){y[s_set] <- 0} #to check if u_i[s] reaches the maximum vaccine uptake. if yes, y[s] <- 0
      if(max(y)==y[s]){
        ###
        u_i[s] <- u_i[s] + 1
        ###
        dy_R_du[s_set] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) ,VacType)
        dy_H_du[s_set] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) * h_i[s], VacType)
        dy_D_du[s_set] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) * d_i[s], VacType)
        ###
        dy_du[s_set] <- dy_D_du[s_set] #change the importance weights of group s for all vaccine types (and therefore "s_set" is used here)
        ###
        dR_du <- dR_du_i[s] #only s (not "s_set" here). This is a scalar (not vector)
        dH_du <- dH_du_i[s] #only s (not "s_set" here). This is a scalar (not vector)
        dD_du <- dD_du_i[s] #only s (not "s_set" here). This is a scalar (not vector)
      }
    }
    ###
    y <- y + dy_du
    ###
    R <- R + dR_du
    H <- H + dH_du
    D <- D + dD_du
    
    dy_R_du <- rep(0,gro*VacType)
    dy_H_du <- rep(0,gro*VacType)
    dy_D_du <- rep(0,gro*VacType)
    ###
    dy_du <- rep(0,gro*VacType)
    ###
    dR_du <- 0
    dH_du <- 0
    dD_du <- 0
  }
}
else {
  print("no target index is specified, and thus no loop")
}

#####4. Store the result#####
z_plot <- round(seq(1, z, length=z_plot_length),0) #Points to be plotted. z_plot means 1:100 % of vaccine stock, when z_plot_length = 100.

result_plot <- as.data.frame(list(NoVac = z_plot, 
                                  PerVac = z_plot/sum(n_i),
                                  Vac_PerAge_PerType = u_loop[z_plot,],
                                  ProVac_PerAge_PerType = u_loop[z_plot,]/n_i,
                                  Imp_PerAge_PerType = y_loop[z_plot,],
                                  R_del=R_loop[z_plot]/R_ini,
                                  H_del=H_loop[z_plot]/H_ini,
                                  D_del=D_loop[z_plot]/D_ini
                                  )
)

#####5. Return the result#####
return(list(res = result_plot, VacCov = u_loop[z_plot,]))

}

#####2. Run simulations for minimizing R/H/D ------
run_R_6gro4vac <- VacAllo_func_4vac(z_plot_length=100, target="R")
run_H_6gro4vac <- VacAllo_func_4vac(z_plot_length=100, target="H")
run_D_6gro4vac <- VacAllo_func_4vac(z_plot_length=100, target="D")

write_rds(run_R_6gro4vac ,"run_R_6gro4vac")
write_rds(run_H_6gro4vac ,"run_H_6gro4vac")
write_rds(run_D_6gro4vac ,"run_D_6gro4vac")

#####3. Store the results for visualization ------
result_objective_6gro4vac <- as.data.frame(list(
  AlloVac = rep((run_R_6gro4vac$res$NoVac/(sum(n_i)/3)),3), 
  R_del = c(run_R_6gro4vac$res$R_del,
            run_H_6gro4vac$res$R_del,
            run_D_6gro4vac$res$R_del
  ),
  H_del = c(run_R_6gro4vac$res$H_del,
            run_H_6gro4vac$res$H_del,
            run_D_6gro4vac$res$H_del
  ),
  D_del = c(run_R_6gro4vac$res$D_del,
            run_H_6gro4vac$res$D_del,
            run_D_6gro4vac$res$D_del
  ),
  Objective = c(rep("Min(infection)",100),
            rep("Min(hospitalization)",100),
            rep("Min(death)",100))
))
result_objective_6gro4vac$Objective <- factor(result_objective_6gro4vac$Objective, levels = c("Min(infection)", "Min(hospitalization)", "Min(death)"))
result_objective_6gro4vac <- dplyr::mutate(result_objective_6gro4vac, AlloVac100 = rep(1:100,3), R_del100=(1-R_del)*100, H_del100=(1-H_del)*100, D_del100=(1-D_del)*100)

write_rds(result_objective_6gro4vac, "result_objective_6gro4vac")

#Vac_PerAge_PerType_R_4vac ----
VacType <- 4
n_i_vec <- rep(csb_pop_6gro$sum,VacType)
result_Vac_R_6gro4vac <- rbind((run_R_6gro4vac$res$Vac_PerAge_PerType.1/n_i_vec[1]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.2/n_i_vec[2]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.3/n_i_vec[3]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.4/n_i_vec[4]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.5/n_i_vec[5]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.6/n_i_vec[6]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.7/n_i_vec[7]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.8/n_i_vec[8]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.9/n_i_vec[9]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.10/n_i_vec[10]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.11/n_i_vec[11]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.12/n_i_vec[12]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.13/n_i_vec[13]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.14/n_i_vec[14]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.15/n_i_vec[15]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.16/n_i_vec[16]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.17/n_i_vec[17]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.18/n_i_vec[18]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.19/n_i_vec[19]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.20/n_i_vec[20]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.21/n_i_vec[21]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.22/n_i_vec[22]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.23/n_i_vec[23]),
                                   (run_R_6gro4vac$res$Vac_PerAge_PerType.24/n_i_vec[24]))
#Vac_PerAge_PerType_H_4vac ----
result_Vac_H_6gro4vac <- rbind((run_H_6gro4vac$res$Vac_PerAge_PerType.1/n_i_vec[1]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.2/n_i_vec[2]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.3/n_i_vec[3]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.4/n_i_vec[4]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.5/n_i_vec[5]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.6/n_i_vec[6]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.7/n_i_vec[7]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.8/n_i_vec[8]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.9/n_i_vec[9]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.10/n_i_vec[10]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.11/n_i_vec[11]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.12/n_i_vec[12]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.13/n_i_vec[13]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.14/n_i_vec[14]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.15/n_i_vec[15]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.16/n_i_vec[16]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.17/n_i_vec[17]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.18/n_i_vec[18]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.19/n_i_vec[19]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.20/n_i_vec[20]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.21/n_i_vec[21]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.22/n_i_vec[22]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.23/n_i_vec[23]),
                                   (run_H_6gro4vac$res$Vac_PerAge_PerType.24/n_i_vec[24]))
#Vac_PerAge_PerType_D_4vac ----
result_Vac_D_6gro4vac <- rbind((run_D_6gro4vac$res$Vac_PerAge_PerType.1/n_i_vec[1]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.2/n_i_vec[2]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.3/n_i_vec[3]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.4/n_i_vec[4]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.5/n_i_vec[5]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.6/n_i_vec[6]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.7/n_i_vec[7]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.8/n_i_vec[8]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.9/n_i_vec[9]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.10/n_i_vec[10]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.11/n_i_vec[11]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.12/n_i_vec[12]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.13/n_i_vec[13]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.14/n_i_vec[14]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.15/n_i_vec[15]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.16/n_i_vec[16]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.17/n_i_vec[17]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.18/n_i_vec[18]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.19/n_i_vec[19]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.20/n_i_vec[20]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.21/n_i_vec[21]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.22/n_i_vec[22]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.23/n_i_vec[23]),
                                   (run_D_6gro4vac$res$Vac_PerAge_PerType.24/n_i_vec[24]))
#Imp_PerAge_PerType_R_4vac ----
result_Imp_R_6gro4vac <- rbind((run_R_6gro4vac$res$Imp_PerAge_PerType.1/n_i_vec[1]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.2/n_i_vec[2]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.3/n_i_vec[3]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.4/n_i_vec[4]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.5/n_i_vec[5]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.6/n_i_vec[6]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.7/n_i_vec[7]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.8/n_i_vec[8]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.9/n_i_vec[9]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.10/n_i_vec[10]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.11/n_i_vec[11]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.12/n_i_vec[12]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.13/n_i_vec[13]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.14/n_i_vec[14]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.15/n_i_vec[15]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.16/n_i_vec[16]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.17/n_i_vec[17]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.18/n_i_vec[18]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.19/n_i_vec[19]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.20/n_i_vec[20]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.21/n_i_vec[21]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.22/n_i_vec[22]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.23/n_i_vec[23]),
                                   (run_R_6gro4vac$res$Imp_PerAge_PerType.24/n_i_vec[24]))
#Imp_PerAge_PerType_H_4vac ----
result_Imp_H_6gro4vac <- rbind((run_H_6gro4vac$res$Imp_PerAge_PerType.1/n_i_vec[1]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.2/n_i_vec[2]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.3/n_i_vec[3]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.4/n_i_vec[4]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.5/n_i_vec[5]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.6/n_i_vec[6]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.7/n_i_vec[7]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.8/n_i_vec[8]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.9/n_i_vec[9]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.10/n_i_vec[10]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.11/n_i_vec[11]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.12/n_i_vec[12]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.13/n_i_vec[13]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.14/n_i_vec[14]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.15/n_i_vec[15]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.16/n_i_vec[16]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.17/n_i_vec[17]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.18/n_i_vec[18]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.19/n_i_vec[19]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.20/n_i_vec[20]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.21/n_i_vec[21]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.22/n_i_vec[22]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.23/n_i_vec[23]),
                                   (run_H_6gro4vac$res$Imp_PerAge_PerType.24/n_i_vec[24]))
#Imp_PerAge_PerType_D_4vac ----
result_Imp_D_6gro4vac <- rbind((run_D_6gro4vac$res$Imp_PerAge_PerType.1/n_i_vec[1]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.2/n_i_vec[2]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.3/n_i_vec[3]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.4/n_i_vec[4]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.5/n_i_vec[5]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.6/n_i_vec[6]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.7/n_i_vec[7]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.8/n_i_vec[8]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.9/n_i_vec[9]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.10/n_i_vec[10]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.11/n_i_vec[11]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.12/n_i_vec[12]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.13/n_i_vec[13]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.14/n_i_vec[14]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.15/n_i_vec[15]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.16/n_i_vec[16]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.17/n_i_vec[17]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.18/n_i_vec[18]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.19/n_i_vec[19]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.20/n_i_vec[20]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.21/n_i_vec[21]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.22/n_i_vec[22]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.23/n_i_vec[23]),
                                   (run_D_6gro4vac$res$Imp_PerAge_PerType.24/n_i_vec[24]))
#write rds----
write_rds(result_Vac_R_6gro4vac, "result_Vac_R_6gro4vac")
write_rds(result_Vac_H_6gro4vac, "result_Vac_H_6gro4vac")
write_rds(result_Vac_D_6gro4vac, "result_Vac_D_6gro4vac")

write_rds(result_Imp_R_6gro4vac, "result_Imp_R_6gro4vac")
write_rds(result_Imp_H_6gro4vac, "result_Imp_H_6gro4vac")
write_rds(result_Imp_D_6gro4vac, "result_Imp_D_6gro4vac")