#######################################################################################################
# Source codes for visualization
### 0. package & color pallet
### 1. read simulation outputs
### 2. plot all results 
### 3. make figures
#######################################################################################################

## 0. package & color pallet -----
library(ggplot2)
library(patchwork)
library(readr)
custom.col.10 <- c("#001219","#005F73","#0A9396","#94D2BD","#E9D8A6","#EE9B00","#CA6702","#BB3E03","#AE2012","#9B2226") #https://coolors.co/001219-005f73-0a9396-94d2bd-e9d8a6-ee9b00-ca6702-bb3e03-ae2012-9b2226 #barplot(1:10, col = custom.col.10)

## 1. read simulation outputs -----
result_list_R_6gro4vac <- read_rds("result_list_R_6gro4vac")
result_list_H_6gro4vac <- read_rds("result_list_H_6gro4vac")
result_list_D_6gro4vac <- read_rds("result_list_D_6gro4vac")

result_objective_6gro4vac <- read_rds("result_objective_6gro4vac")
result_objective_6gro3vac <- read_rds("result_objective_6gro3vac")

VE_dat_list <- read_rds("VE_dat_list")
Input_dat_list <- read_rds("Input_dat_list")

summary_inf_sim <- read_rds("summary_inf_sim")
test_nonvac <- read_rds("test_nonvac")
## 2. plot all results -----
#plot all results in minimize(R) -----
plot_R_vac <- list()
j <- 1
for (i in c("<20", "21-30", "31-40", "41-50", "51-60", "60+")) {
  plot_R_vac[[j]] <- ggplot(
    data = result_list_R_6gro4vac %>% filter(Age == i),
    mapping = aes(x = Time, y = value, fill = VaccineType)
  ) +
    geom_bar(stat = "identity") +
    ggtitle(paste("Age", i)) +
    labs(x = "Allocated vaccine [%]", y = "Proportion vaccinated") +
    scale_fill_manual(values = rev(custom.col.10[1:4])) +
    theme_minimal() #+
  # theme(legend.position = "none")
  j <- j + 1
}

plot_R_vac_all <- ((plot_R_vac[[1]] | plot_R_vac[[2]] | plot_R_vac[[3]]) / (plot_R_vac[[4]] | plot_R_vac[[5]] | plot_R_vac[[6]]))
###
# "proportion vaccinated" vs "allocated vaccine[%]" when minimizing R
plot_R_6gro4vac_S3 <- ggplot(
  data = result_list_R_6gro4vac %>% group_by(Age, Time) %>% mutate(CumDose = sum(value)),
  mapping = aes(x = Time, y = CumDose, col = Age)
) +
  geom_point() +
  geom_line() +
  labs(x = "Allocated vaccine [%]", y = "Proportion vaccinated") +
  scale_color_manual(values = custom.col.10[5:10]) + 
  theme_minimal()
###
# barplot: "cumulative proportion vaccinated" vs "allocated vaccine[%]" when minimizing R, by vaccine type (A) and by age-group()
plot_R_6gro4vac_S2A <- ggplot(
  data = result_list_R_6gro4vac,
  mapping = aes(x = Time, y = RelativeVac, fill = VaccineType)
) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = rev(custom.col.10[1:4])) +
  labs(x = "Allocated vaccine [%]", y = "Proportion vaccinated") + 
  theme_minimal()

plot_R_6gro4vac_S2B <- ggplot(
  data = result_list_R_6gro4vac,
  mapping = aes(x = Time, y = RelativeVac, fill = Age)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom.col.10[5:10]) +
  labs(x = "Allocated vaccine [%]", y = "Proportion vaccinated") + 
  theme_minimal()

plot_R_6gro4vac_S2 <- plot_R_6gro4vac_S2A | plot_R_6gro4vac_S2B
#plot all results in minimize(H) -----
plot_H_vac <- list()
j <- 1
for (i in c("<20", "21-30", "31-40", "41-50", "51-60", "60+")) {
  plot_H_vac[[j]] <- ggplot(
    data = result_list_H_6gro4vac %>% filter(Age == i),
    mapping = aes(x = Time, y = value, fill = VaccineType)
  ) +
    geom_bar(stat = "identity") +
    ggtitle(paste("Age", i)) +
    labs(x = "Allocated vaccine [%]", y = "Proportion vaccinated") +
    scale_fill_manual(values = rev(custom.col.10[1:4])) +
    theme_minimal()
  j <- j + 1
}

plot_H_vac_all <- ((plot_H_vac[[1]] | plot_H_vac[[2]] | plot_H_vac[[3]]) / (plot_H_vac[[4]] | plot_H_vac[[5]] | plot_H_vac[[6]]))

###
# "proportion vaccinated" vs "allocated vaccine[%]" when minimizing H
plot_H_6gro4vac_S3 <- ggplot(
  data = result_list_H_6gro4vac %>% group_by(Age, Time) %>% mutate(CumDose = sum(value)),
  mapping = aes(x = Time, y = CumDose, col = Age)
) +
  geom_point() +
  geom_line() +
  labs(x = "Allocated vaccine [%]", y = "Proportion vaccinated") +
  scale_color_manual(values = custom.col.10[5:10]) + 
  theme_minimal()

###
# barplot: "cumulative proportion vaccinated" vs "allocated vaccine[%]" when minimizing H, by vaccine type (A) and by age-group()
plot_H_6gro4vac_S2A <- ggplot(
  data = result_list_H_6gro4vac,
  mapping = aes(x = Time, y = RelativeVac, fill = VaccineType)
) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = rev(custom.col.10[1:4])) +
  labs(x = "Allocated vaccine [%]", y = "Proportion vaccinated") + 
  theme_minimal()

plot_H_6gro4vac_S2B <- ggplot(
  data = result_list_H_6gro4vac,
  mapping = aes(x = Time, y = RelativeVac, fill = Age)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom.col.10[5:10]) + 
  labs(x = "Allocated vaccine [%]", y = "Proportion vaccinated") + 
  theme_minimal()

plot_H_6gro4vac_S2 <- plot_H_6gro4vac_S2A | plot_H_6gro4vac_S2B


#plot all results in minimize(D) -----
plot_D_vac <- list()
j <- 1
for (i in c("<20", "21-30", "31-40", "41-50", "51-60", "60+")) {
  plot_D_vac[[j]] <- ggplot(
    data = result_list_D_6gro4vac %>% filter(Age == i),
    mapping = aes(x = Time, y = value, fill = VaccineType)
  ) +
    geom_bar(stat = "identity") +
    ggtitle(paste("Age", i)) +
    labs(x = "Allocated vaccine [%]", y = "Proportion vaccinated") +
    scale_fill_manual(values = rev(custom.col.10[1:4])) +
    theme_minimal()
  j <- j + 1
}

plot_D_vac_all <- ((plot_D_vac[[1]] | plot_D_vac[[2]] | plot_D_vac[[3]]) / (plot_D_vac[[4]] | plot_D_vac[[5]] | plot_D_vac[[6]]))

###
# "proportion vaccinated" vs "allocated vaccine[%]" when minimizing D
plot_D_6gro4vac_S3 <- ggplot(
  data = result_list_D_6gro4vac %>% group_by(Age, Time) %>% mutate(CumDose = sum(value)),
  mapping = aes(x = Time, y = CumDose, col = Age)
) +
  geom_point() +
  geom_line() +
  labs(x = "Allocated vaccine [%]", y = "Proportion vaccinated") +
  scale_color_manual(values = custom.col.10[5:10]) + 
  theme_minimal()
###
# barplot: "cumulative proportion vaccinated" vs "allocated vaccine[%]" when minimizing D, by vaccine type (A) and by age-group()

plot_D_6gro4vac_S2A <- ggplot(
  data = result_list_D_6gro4vac,
  mapping = aes(x = Time, y = RelativeVac, fill = VaccineType)
) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = rev(custom.col.10[1:4])) +
  labs(x = "Allocated vaccine [%]", y = "Proportion vaccinated") + 
  theme_minimal()

plot_D_6gro4vac_S2B <- ggplot(
  data = result_list_D_6gro4vac,
  mapping = aes(x = Time, y = RelativeVac, fill = Age)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom.col.10[5:10]) + 
  labs(x = "Allocated vaccine [%]", y = "Proportion vaccinated") + 
  theme_minimal()

plot_D_6gro4vac_S2 <- plot_D_6gro4vac_S2A | plot_D_6gro4vac_S2B


## 3. make figures -----
#Fig-1 ----
#Fig-1(A)(B)(C)
Epicurve_panel <-  ggplot(data = summary_inf_sim %>% filter(Strategy=="No vaccination", time <201),
                          mapping = aes(x = time, y = inf, col = Strategy)) +
  geom_line() +
  labs(x="Time [day]", y = "Incidence of Infection") +
  theme_minimal() +
  geom_line(size=1.2)

incidence_panel <- ggplot(data = test_nonvac %>% filter(age!="all", time <201),#, time >41
                          mapping = aes(x = time, y =age)) +
  geom_raster(aes(fill=inf)) +
  scale_fill_gradient(low="grey90", high="red", name="Incidence of infection") +
  labs(x="Time [day]", y="Age") +
  theme_minimal() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                          axis.text.y=element_text(size=9),
                          plot.title=element_text(size=11))

FOI_panel <- ggplot(data = test_nonvac %>% filter(age!="all", time <201),
                    mapping = aes(x = time, y =age)) +
  geom_raster(aes(fill=FOI)) +
  scale_fill_gradient(low="grey90", high="red", name="Force of infection") +
  labs(x="Time [day]", y="Age") +
  theme_minimal() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                          axis.text.y=element_text(size=9),
                          plot.title=element_text(size=11))

#Fig-1(D)(E)(F)
inf_summary_age_sir_plot <- ggplot(data = summary_inf_sim %>% filter(time <201),
                                   mapping = aes(x = time, y = inf, col = Strategy)) +
  geom_line() +
  labs(x="Time [day]", y = "Incidence of Infection") +
  theme_minimal() +
  geom_line(size=1.2)

hosp_summary_age_sir_plot <- ggplot(data = summary_inf_sim %>% filter(time <201),
                                    mapping = aes(x = time, y = hosp, col = Strategy)) +
  geom_line() +
  labs(x="Time [day]", y = "Incidence of Hospitalization") +
  theme_minimal() +
  geom_line(size=1.2)

death_summary_age_sir_plot <- ggplot(data = summary_inf_sim %>% filter(time <201),
                                     mapping = aes(x = time, y = death, col = Strategy)) +
  geom_line() +
  labs(x="Time [day]", y = "Incidence of Death") +
  theme_minimal() +
  geom_line(size=1.2)

#Fig-1 Merged
SimEpi_marge <- (Epicurve_panel/incidence_panel/FOI_panel)
Fig_1 <- ((SimEpi_marge)|(inf_summary_age_sir_plot / hosp_summary_age_sir_plot / death_summary_age_sir_plot)) + 
  plot_annotation(tag_levels = 'A') #1_ageSIRsimulation.R: H8.13 x W17
#Fig-2 ----
Fig2_R <- (plot_R_vac[[1]] + theme(axis.title.y = element_blank()) |
  plot_R_vac[[2]] + theme(axis.text.y = element_blank(), axis.title.y = element_blank()) |
  plot_R_vac[[3]] + theme(axis.text.y = element_blank(), axis.title.y = element_blank()) |
  plot_R_vac[[4]] + theme(axis.text.y = element_blank(), axis.title.y = element_blank()) |
  plot_R_vac[[5]] + theme(axis.text.y = element_blank(), axis.title.y = element_blank()) |
  plot_R_vac[[6]] + theme(axis.text.y = element_blank(), axis.title.y = element_blank())) +
  plot_layout(guides = "collect") & theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank())
Fig2_H <- (plot_H_vac[[1]] |
  plot_H_vac[[2]] + theme(axis.text.y = element_blank(), legend.position = "none") |
  plot_H_vac[[3]] + theme(axis.text.y = element_blank(), legend.position = "none") |
  plot_H_vac[[4]] + theme(axis.text.y = element_blank(), legend.position = "none") |
  plot_H_vac[[5]] + theme(axis.text.y = element_blank(), legend.position = "none") |
  plot_H_vac[[6]] + theme(axis.text.y = element_blank(), legend.position = "none")) +
  plot_layout(guides = "collect") & theme(axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank())
Fig2_D <- (plot_D_vac[[1]] |
  plot_D_vac[[2]] + theme(axis.text.y = element_blank()) |
  plot_D_vac[[3]] + theme(axis.text.y = element_blank()) |
  plot_D_vac[[4]] + theme(axis.text.y = element_blank()) |
  plot_D_vac[[5]] + theme(axis.text.y = element_blank()) |
  plot_D_vac[[6]] + theme(axis.text.y = element_blank())) +
  plot_layout(guides = "collect") & theme(legend.position = "none", axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank())
Fig_2_merge <- (Fig2_R / Fig2_H / Fig2_D)
gt <- patchwork::patchworkGrob(Fig_2_merge)
Fig_2 <- gridExtra::grid.arrange(gt, left = "Proportion vaccinated", bottom = "Allocated Vaccines [%]")
#Fig-3 ----
plot_R_4vac <- ggplot(data = result_objective_6gro4vac, 
                                  mapping = aes(x = AlloVac100, y = R_del100, col = Objective)) +
  geom_line(size=1.25) + 
  labs(x="Allocated vaccine [%]", y = "Reduction in infection [%]") +
  scale_color_manual(values = c("#E63946", "#A8DADC", "#1D3557"))+
  ylim(0,100) +
  theme_minimal()

plot_H_4vac <- ggplot(data = result_objective_6gro4vac, 
                                  mapping = aes(x = AlloVac100, y = H_del100, col = Objective)) +
  geom_line(size=1.25) + 
  labs(x="Allocated vaccine [%]", y = "Reduction in hospitalization [%]") +
  scale_color_manual(values = c("#E63946", "#A8DADC", "#1D3557"))+
  ylim(0,100) +
  theme_minimal()

plot_D_4vac <- ggplot(data = result_objective_6gro4vac, 
                                  mapping = aes(x = AlloVac100, y = D_del100, col = Objective)) +
  geom_line(size=1.25) + 
  labs(x="Allocated vaccine [%]", y = "Reduction in death [%]") +
  scale_color_manual(values = c("#E63946", "#A8DADC", "#1D3557"))+
  ylim(0,100) +
  theme_minimal()

Fig_3 <- (plot_R_4vac | plot_H_4vac | plot_D_4vac) + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect") & theme(legend.position = "right")
#Fig-S1 -----
p1 <- ggplot(Input_dat_list, aes(x = age, y = pop)) +
  geom_bar(stat = "identity") +
  ggtitle("Population structure(2019)") +
  labs(x = "Age class", y = "Number of people (x 100,000)") +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1))

p2 <- ggplot(Input_dat_list, aes(x = age, y = sero)) +
  geom_bar(stat = "identity") +
  ggtitle("Seroprevalence (early June 2020)") +
  labs(x = "Age class", y = "Percentage seropositive") +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1))

p3 <- ggplot(Input_dat_list, aes(x = age, y = inci)) +
  geom_bar(stat = "identity") +
  ggtitle("Incidence of notified cases*") +
  labs(x = "Age class", y = "Notification per 100,000 inhabitats") +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1)) # plot.title = element_text(size=8),

p4 <- ggplot(VE_dat_list, aes(x = age, y = VE, fill = VaccineType)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Vaccine efficacy") +
  labs(x = "Age class", y = "Proportion protected") +
  scale_fill_manual(values = rev(custom.col.10[1:4])) +
  ylim(0, 1) +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1))

p5 <- ggplot(Input_dat_list, aes(x = age, y = will)) +
  geom_bar(stat = "identity") +
  ggtitle("Max vaccine uptake") +
  labs(x = "Age class", y = "Maximum proportion vaccinated") +
  ylim(0, 1) +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1))

p_hosp <- ggplot(Input_dat_list, aes(x = age, y = HospRate)) +
  geom_bar(stat = "identity") +
  ggtitle("Infection hospitalization rate") +
  labs(x = "Age class", y = "Infection hospitalization rate (%)") +
  # ylim(0,1) +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1))

p_death <- ggplot(Input_dat_list, aes(x = age, y = DeathRate)) +
  geom_bar(stat = "identity") +
  ggtitle("Infection fatality rate") +
  labs(x = "Age class", y = "Infection fatality rate (%)") +
  # ylim(0,1) +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1))

Fig_S1 <- (p1 | p2 | p3) / (p4 | p5 | p_hosp | p_death) + plot_annotation(tag_levels = "A")
#Fig-S2 ----
Fig_S2 <- plot_R_6gro4vac_S2 / plot_H_6gro4vac_S2 / plot_D_6gro4vac_S2 + plot_annotation(tag_levels = "A")
#Fig-S3 ----
Fig_S3 <- plot_R_6gro4vac_S3 / plot_H_6gro4vac_S3 / plot_D_6gro4vac_S3 + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect") & theme(legend.position = "right")
#Fig-S4 ----
plot_R_3vac <- ggplot(data = result_objective_6gro3vac, 
                      mapping = aes(x = AlloVac100, y = R_del100, col = Objective)) +
  geom_line(size=1.25) + 
  labs(x="Allocated vaccine [%]", y = "Reduction in infection [%]") +
  scale_color_manual(values = c("#E63946", "#A8DADC", "#1D3557"))+
  ylim(0,100) +
  theme_minimal()

plot_H_3vac <- ggplot(data = result_objective_6gro3vac, 
                      mapping = aes(x = AlloVac100, y = H_del100, col = Objective)) +
  geom_line(size=1.25) + 
  labs(x="Allocated vaccine [%]", y = "Reduction in hospitalization [%]") +
  scale_color_manual(values = c("#E63946", "#A8DADC", "#1D3557"))+
  ylim(0,100) +
  theme_minimal()

plot_D_3vac <- ggplot(data = result_objective_6gro3vac, 
                      mapping = aes(x = AlloVac100, y = D_del100, col = Objective)) +
  geom_line(size=1.25) + 
  labs(x="Allocated vaccine [%]", y = "Reduction in death [%]") +
  scale_color_manual(values = c("#E63946", "#A8DADC", "#1D3557"))+
  ylim(0,100) +
  theme_minimal()

Fig_S4 <- (plot_R_3vac | plot_H_3vac | plot_D_3vac) + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect") & theme(legend.position = "right")