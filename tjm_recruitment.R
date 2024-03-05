library(wham)
library(dplyr)
asap3 <- read_asap3_dat("haddock.dat")

# bias_correct_process     <- c(TRUE,FALSE)
# bias_correct_observation <- c(TRUE,FALSE)
# bias_correct_BRPs        <- c(TRUE,FALSE)
# df.mods <- expand.grid(bias_correct_process = bias_correct_process, 
#                        bias_correct_observation = bias_correct_observation, 
#                        bias_correct_BRPs = bias_correct_BRPs, 
#                        stringsAsFactors = FALSE)

# n.mods <- dim(df.mods)[1] #8 scenarios!
# df.mods$Model <- paste0("m_",1:n.mods)
# df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
nsim = 100
# sim_input = list()

##########################################
#recruitment only
inputFF <- prepare_wham_input(asap3, 
  recruit_model=2,
  selectivity=list(model=rep("age-specific",5), 
    re=rep("none",5), 
    initial_pars=list(c(0.01, 0.1, 0.3, 1, 1, 1, 1, 1, 1),
      c(0.01, 0.1, 0.3, 0.5, 0.8, 0.9, 1, 1, 1),
      c(0.01, 0.1, 0.3, 0.5, 0.8, 0.9, 1, 1, 1),
      c(0.01, 0.1, 0.3, 1, 1, 1, 1, 1, 1),
      c(0.01, 0.1, 0.3, 0.5, 0.8, 1, 1, 1, 1)), 
    fix_pars=list(4:9,7:9,7:9,4:9,6:9)),
  NAA_re = list(sigma="rec", cor="iid"),
  age_comp = "logistic-normal-miss0",
  model_name="FF",
  basic_info = list(
    bias_correct_process = FALSE,
    bias_correct_observation = FALSE,
    simulate_process_error = rep(TRUE, 5),
    #XSPR_R_avg_yrs = 50, #only need 1 year because expected recruitment is constant
    XSPR_R_opt = 5)
)
modFF <- fit_wham(inputFF, do.fit = TRUE, do.osa = F, do.retro = F) 
inputFF_true <- inputFF
inputFF_true$par <- modFF$parList
inputFF_true$par$log_NAA[] <- 10
modFF_true <- fit_wham(inputFF_true, do.fit = F, do.osa = F, do.retro = F) 
modFF_true$env$parList()$log_NAA #no inner optimization
modFF_true$fn()
modFF_true$env$parList()$log_NAA #inner optimization completed
range(modFF_true$env$parList()$log_NAA - modFF$parList$log_NAA) #essentially 0

modFF_true_est <- fit_wham(inputFF_true, do.fit = T, do.osa = F, do.retro = F) 

inputTT <- inputFF
inputTT$par <- modFF$parList
inputTT$data$bias_correct_pe <- 1 
inputTT$map <- lapply(inputTT$par, function(x) factor(rep(NA, length(x))))
inputTT$map$mean_rec_pars <- factor(1:length(inputTT$par$mean_rec_pars)) #just estimate mean par 
inputTT$map$log_NAA <- inputFF$map$log_NAA
modTT <- fit_wham(inputTT, do.fit = TRUE, do.osa = F, do.retro = F) 

#what if all fixed effects are estimated, starting at the MLEs
inputTT_est_all <- modTT$input
inputTT_est_all$par <- modTT$parList
inputTT_est_all$map <- inputFF$map
modTT_est_all <- fit_wham(inputTT_est_all, do.fit = TRUE, do.osa = F, do.retro = F) 
modTT_est_all$opt$obj  - modFF$opt$obj  #the same
modFF$opt$par - modTT_est_all$opt$par #the same except mean_rec_pars

exp(modTT$parList$mean_rec_pars)
exp(modFF$parList$mean_rec_pars) 

exp(modTT_est_all$parList$mean_rec_pars)
exp(modFF$parList$mean_rec_pars + 0.5*exp(modFF$parList$log_NAA_sigma)^2) 
exp(modTT_est_all$parList$mean_rec_pars - 0.5*exp(modTT_est_all$parList$log_NAA_sigma)^2)

#identical
modTT_est_all$parList$log_NAA[,1]
modFF$parList$log_NAA[,1]
range(modTT_est_all$parList$log_NAA[,1] - modFF$parList$log_NAA[,1])
#confirm RE are actually estimated
as.list(modTT_est_all$sdrep, "Std")$log_NAA[,1]
as.list(modFF$sdrep, "Std")$log_NAA[,1]

#use default initial values
inputTT_est_all_2 <- inputFF
inputTT_est_all_2$data$bias_correct_pe <- 1 
modTT_est_all_2 <- fit_wham(inputTT_est_all_2, do.fit = TRUE, do.osa = F, do.retro = F) 
modTT_est_all_2$opt$obj  - modFF$opt$obj  #the same
modFF$opt$par - modTT_est_all_2$opt$par #the same except mean_rec_pars

exp(modTT_est_all_2$parList$mean_rec_pars)
exp(modFF$parList$mean_rec_pars + 0.5*exp(modFF$parList$log_NAA_sigma)^2) 
exp(modTT_est_all_2$parList$mean_rec_pars - 0.5*exp(modTT_est_all_2$parList$log_NAA_sigma)^2)

#identical
modTT_est_all_2$parList$log_NAA[,1]
modFF$parList$log_NAA[,1]
range(modTT_est_all_2$parList$log_NAA[,1] - modFF$parList$log_NAA[,1])
#confirm RE are actually estimated
as.list(modTT_est_all_2$sdrep, "Std")$log_NAA[,1]
as.list(modFF$sdrep, "Std")$log_NAA[,1]

simsTT <- sapply(1:1000, function(x) modTT$simulate()$log_NAA[,1])
mean(exp(simsTT) - exp(modTT$parList$mean_rec_pars))
mean(exp(simsTT- modTT$parList$mean_rec_pars)) # = 1

simsFF <- sapply(1:100, function(x) modFF$simulate()$log_NAA[,1])
mean(exp(simsFF) - exp(modFF$parList$mean_rec_pars))
##########################################




nsims <- 100
#estimate given true Rs are simulated without bias-correction
estR_FT <- estR_FF <- trueR_F <- matrix(NA, inputFF$data$n_years_model-1, nsims)
#estimate given true Rs are simulated with bias-correction
estR_TT <- estR_TF <- trueR_T <- matrix(NA, inputFF$data$n_years_model-1, nsims)
set.seed(123)
for(OM_bc in c("T","F")) {
  if(OM_bc == "F") {
    OM_input <- modFF$input
    OM_input$par <- modFF$parList
  } else {
    OM_input <- modTT$input
    OM_input$par <- modTT$parList
  }
  OM <- fit_wham(OM_input, do.fit = F, do.osa = F, do.retro = F) 
  OM$fn()
  for(i in 1:nsims){
    sim <- OM$simulate(complete=T)
    if(OM_bc == "T"){
      trueR_T[,i] <- sim$log_NAA[,1]
    } else {
      trueR_F[,i] <- sim$log_NAA[,1]
    }
    input_i <- OM$input
    input_i$data <- sim
    for(EM_bc in c("T","F")){
      if(EM_bc == "T"){ 
        input_i$data$bias_correct_pe <- 1
      } else {
        input_i$data$bias_correct_pe <- 0
      }
      mod_i <- fit_wham(input_i, do.fit = F, do.osa = F, do.retro = F) 
      mod_i$fn()
      if(EM_bc == "T"){
        if(OM_bc == "T") estR_TT[,i] <- mod_i$env$parList()$log_NAA[,1]
        else estR_FT[,i] <- mod_i$env$parList()$log_NAA[,1]
      } else{
        if(OM_bc == "T") estR_TF[,i] <- mod_i$env$parList()$log_NAA[,1]
        else estR_FF[,i] <- mod_i$env$parList()$log_NAA[,1]
      }
    }
  }
}


#############################################
par(mfrow = c(2,2), oma = c(1,5,2,1), mar = c(5,3,1,3))
i <- 40
plot(exp(estR_FT[,i]), ylab = "Recruitment", xlab = "Year", type = "l")
lines(exp(trueR_F[,i]), col = 'red')
legend("topleft", col = c("black","red"), legend = c("Estimated", "Simulated"), lty = 1)

title("Simulated without bias-correction", outer=F, line = 1, xpd = NA)

plot(exp(estR_TT[,i]), ylab = "Recruitment", xlab = "Year", type = "l")
lines(exp(trueR_T[,i]), col = 'red')
mtext(side = 4, "Estimated with bias-correction", outer = F, line = 1)
title("Simulated with bias-correction", outer=F, line = 1, xpd = NA)

plot(exp(estR_FF[,i]), ylab = "Recruitment", xlab = "Year", type = "l")
lines(exp(trueR_F[,i]), col = 'red')

plot(exp(estR_TF[,i]), ylab = "Recruitment", xlab = "Year", type = "l")
lines(exp(trueR_T[,i]), col = 'red')
mtext(side = 4, "Estimated without bias-correction", outer = F, line = 1)

title("Examples of simulated and estimated recruitment", outer=T, line = 1)

#############################################
par(mfrow = c(1,2), oma = c(1,5,2,1), mar = c(5,3,1,1))
plot(exp(estR_FT[,i]-trueR_F[,i]), ylab = "Estimated/TRUE", xlab = "Year", type = "l")
lines(exp(estR_FF[,i]-trueR_F[,i]), col = 'red')
legend("topleft", col = c("black","red"), legend = c(" with b-c", "without b-c"), lty = 1)

title("Simulated without bias-correction", outer=F, line = 1, xpd = NA)

plot(exp(estR_TT[,i]-trueR_T[,i]), ylab = "Estimated/TRUE", xlab = "Year", type = "l")
lines(exp(estR_TF[,i]-trueR_T[,i]), col = 'red')
legend("topleft", col = c("black","red"), legend = c(" with b-c", "without b-c"), lty = 1)
title("Simulated with bias-correction", outer=F, line = 1, xpd = NA)

#title("Examples of simulated and estimated recruitment", outer=T, line = 1)


#############################################

apply(exp(trueR_F), 2, mean)
apply(exp(trueR_T), 2, mean)
#
mean(exp(estR_FT) - exp(trueR_F))
mean(exp(estR_FF) - exp(trueR_F))
mean(exp(estR_TT) - exp(trueR_T))
mean(exp(estR_TF) - exp(trueR_T))

apply(exp(estR_FT) - exp(trueR_F),1,mean)
apply(exp(estR_FF) - exp(trueR_F),1,mean)

apply(exp(estR_TT) - exp(trueR_T),1,mean)
apply(exp(estR_TF) - exp(trueR_T),1,mean)

par(mfrow = c(1,2))
plot(apply(exp(estR_FT) - exp(trueR_F), 1, mean), ylim =c(0.96,1.2), ylab = "mean ratio of Estimated (b-c) and TRUE (no b-c)", xlab = "Year")
grid(col = gray(0.7))
abline(h = 1)
lines(supsmu(1:NROW(estR_FF), apply(exp(estR_FT)/exp(trueR_F), 1, mean)))

plot(apply(exp(estR_FF)/exp(trueR_F), 1, mean), ylim =c(0.96,1.2), ylab = "mean ratio of Estimated (no b-c) and TRUE (no b-c)", xlab = "Year")
grid(col = gray(0.7))
abline(h = 1)
lines(supsmu(1:NROW(estR_FF), apply(exp(estR_FF)/exp(trueR_F), 1, mean)))

par(mfrow = c(1,2))
plot(apply(exp(estR_TT)/exp(trueR_T), 1, mean), ylim =c(0.96,1.2), ylab = "mean ratio of Estimated (b-c) and TRUE (b-c)", xlab = "Year")
grid(col = gray(0.7))
abline(h = 1)
lines(supsmu(1:NROW(estR_TT), apply(exp(estR_TT)/exp(trueR_T), 1, mean)))

plot(apply(exp(estR_TF)/exp(trueR_T), 1, mean), ylim =c(0.96,1.2), ylab = "mean ratio of Estimated (no b-c) and TRUE (b-c)", xlab = "Year")
grid(col = gray(0.7))
abline(h = 1)
lines(supsmu(1:NROW(estR_TT), apply(exp(estR_TF)/exp(trueR_T), 1, mean)))

apply(exp(estR_FT), 2, mean)
plot(apply(exp(trueR_F), 2, mean)/apply(exp(estR_FT), 2, mean))
plot(apply(exp(estR_FF), 2, mean)/apply(exp(estR_FT), 2, mean))
plot(apply(exp(trueR_T), 2, mean)/apply(exp(estR_TT), 2, mean))

