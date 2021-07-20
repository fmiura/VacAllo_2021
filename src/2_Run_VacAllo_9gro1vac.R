#######################################################################################################
# Source codes for the allocation algorithm (9 age groups, 1 vaccine type)
###0. Run "0_InputData.R" and read output data 
###1. Define function to optimize vaccine allocations for different objectives
###2. Run simulations for minimizing R/H/D
###3. Store the results for age-SIR simulations (40% coverage)
#######################################################################################################

#####0. Run "0_InputData.R" and read output data ------
incidence_time15to45 <- read_rds("incidence_time15to45")
susceptible_time15to45 <- read_rds("susceptible_time15to45")

#####1. Define function to optimize vaccine allocations for different objectives ------
VacAllo_func_gro9vac1 <- function(z_plot_length=100, target="R"){
#####1. input data#####
gro <- 9   #No. of stratified group
VacType <- 1 #No. of vaccine types (Pfizer, Moderena, AstraZeneca, Jansen)

#Total population
n_i <- n_i
#Incidence
x_i <- incidence_time15to45
#Susceptibles
s_i <- susceptible_time15to45
#per contact probability of acquiring infection for group i
a_i <- rep(1,gro) 
#per contact infectiousness of group i 
c_i <- rep(1,gro)
#Maximum uptake (willingness to vaccination)
Uptake_rate <- rep(0.8, gro) #40% for all age groups, "gro" age-groups 
#vaccine efficacy 
q_i <- rep(0.95,gro)
#Infection hospitalization rate
h_i <- h_i
#Infection mortality rate
d_i <- d_i

#####2. Importance weight#####
### Define the approximated next generation matrix 
#top right and left eigenvectors
w_1 <- (1 / sum(x_i)) * x_i
v_1 <- (sum(x_i) / sum(x_i^2/s_i)) * (x_i / s_i)
#approximated K
R_ini <- 1.2
approxK <- R_ini * w_1 %x% t(v_1)

### Target indexes
#initial values
R_ini <- 1.2
H_ini <- Re(eigen(diag(h_i) %*% approxK)$values[1])
D_ini <- Re(eigen(diag(d_i) %*% approxK)$values[1])

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
#allocated vaccine for group i
u_i <- rep(0,gro*VacType)
#total vaccine stock
z <- sum(n_i*Uptake_rate)
z <- round(z, 0)

#initial values 
y_R_ini <-  R * (h) * q_i * (c_i / a_i) * (x_i / s_i) * (x_i / n_i) #((q_i * (c_i / a_i) * (n_i / s_i) )^(1/2)) * (x_i / n_i) 
y_H_ini <-  R * (h) * q_i * (c_i / a_i) * (x_i / s_i) * (x_i / n_i) * h_i#(h_i^(1/2))
y_D_ini <-  R * (h) * q_i * (c_i / a_i) * (x_i / s_i) * (x_i / n_i) * d_i#(d_i^(1/2))

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
    
    for (s in 1:(gro*VacType)) {
      if(u_i[s] >= (n_i[s]*Uptake_rate[s])){y[s] <- 0} #to check if u_i[s] reaches the maximum vaccine uptake. if yes, y[s] <- 0
      if(max(y)==y[s]){
        ###
        u_i[s] <- u_i[s] + 1
        ###
        dy_R_du[s] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) ,VacType)
        dy_H_du[s] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) * h_i[s], VacType)
        dy_D_du[s] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) * d_i[s], VacType)
        ###
        dy_du[s] <- dy_R_du[s] 
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
    
    for (s in 1:(gro*VacType)) {
      if(u_i[s] >= (n_i[s]*Uptake_rate[s])){y[s] <- 0} #to check if u_i[s] reaches the maximum vaccine uptake. if yes, y[s] <- 0
      if(max(y)==y[s]){
        u_i[s] <- u_i[s] + 1
        ###
        dy_R_du[s] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) ,VacType)
        dy_H_du[s] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) * h_i[s], VacType)
        dy_D_du[s] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) * d_i[s], VacType)
        ###
        dy_du[s] <- dy_H_du[s] 
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

    
    for (s in 1:(gro*VacType)) {
      if(u_i[s] >= (n_i[s]*Uptake_rate[s])){y[s] <- 0} #to check if u_i[s] reaches the maximum vaccine uptake. if yes, y[s] <- 0
      if(max(y)==y[s]){
        ###
        u_i[s] <- u_i[s] + 1
        ###
        dy_R_du[s] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) ,VacType)
        dy_H_du[s] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) * h_i[s], VacType)
        dy_D_du[s] <- rep(- (q_i[s] / n_i[s]) * R * (h) * q_i[s] * (c_i[s] / a_i[s]) * (x_i[s] / s_i[s]) * (x_i[s] / n_i[s]) * d_i[s], VacType)
        ###
        dy_du[s] <- dy_D_du[s] 
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
z_plot <- round(seq(1, z, length=z_plot_length),0) #Points to be plotted

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
run_R_gro9vac1 <- VacAllo_func_gro9vac1(z_plot_length=100, target="R")
run_H_gro9vac1 <- VacAllo_func_gro9vac1(z_plot_length=100, target="H")
run_D_gro9vac1 <- VacAllo_func_gro9vac1(z_plot_length=100, target="D")

#####3. Store the results for age-SIR simulations (40% coverage) ------
Vac40dist_minR <- run_R_gro9vac1$VacCov[51,]/n_i
Vac40dist_minH <- run_H_gro9vac1$VacCov[51,]/n_i
Vac40dist_minD <- test_run_D_gro9vac1$VacCov[51,]/n_i

write_rds(Vac40dist_minR, "Vac40dist_minR")
write_rds(Vac40dist_minH, "Vac40dist_minH")
write_rds(Vac40dist_minD, "Vac40dist_minD")