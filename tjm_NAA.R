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
#recruitment + survival
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
  NAA_re = list(sigma="rec+1", cor="iid"),
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

inputTT <- inputFF
inputTT$par <- modFF$parList
inputTT$data$bias_correct_pe <- 1 
inputTT$map <- lapply(inputTT$par, function(x) factor(rep(NA, length(x)))) #fix all the pars
inputTT$map$mean_rec_pars <- factor(1:length(inputTT$par$mean_rec_pars)) #just estimate mean par 
inputTT$map$log_NAA <- inputFF$map$log_NAA
modTT <- fit_wham(inputTT, do.fit = TRUE, do.osa = F, do.retro = F) 
modTT$opt$obj  - modFF$opt$obj  #NOT the same

#NOT identical
modTT$parList$log_NAA
modFF$parList$log_NAA
range(modTT$parList$log_NAA - modFF$parList$log_NAA)
#confirm RE are actually estimated
as.list(modTT$sdrep, "Std")$log_NAA[,1]
as.list(modFF$sdrep, "Std")$log_NAA[,1]

mean(exp(modTT$rep$NAA_devs))

#what if all fixed effects are estimated, starting at the MLEs
inputTT_est_all <- modTT$input
inputTT_est_all$par <- modTT$parList
inputTT_est_all$map <- inputFF$map
modTT_est_all <- fit_wham(inputTT_est_all, do.fit = TRUE, do.osa = F, do.retro = F) 
mean(exp(modTT_est_all$rep$NAA_devs))

mean(exp(modTT$simulate()$NAA_devs[,-1])) # = 1
mean(exp(modFF$simulate()$NAA_devs[,-1])) # != 1

exp(modTT$parList$mean_rec_pars)
exp(modFF$parList$mean_rec_pars) 

exp(modTT_est_all$parList$mean_rec_pars)
exp(modFF$parList$mean_rec_pars + 0.5*exp(modFF$parList$log_NAA_sigma[1])^2) 
exp(modTT_est_all$parList$mean_rec_pars - 0.5*exp(modTT_est_all$parList$log_NAA_sigma[1])^2)

#NOT identical
modTT_est_all$parList$log_NAA
modFF$parList$log_NAA
range(modTT$parList$log_NAA - modFF$parList$log_NAA)
range(modTT_est_all$parList$log_NAA - modFF$parList$log_NAA)
apply(modTT$parList$log_NAA - modFF$parList$log_NAA,2,range)
apply(modTT_est_all$parList$log_NAA - modFF$parList$log_NAA,2,range)

simTT <- modTT$simulate()
range(log(simTT$NAA/simTT$pred_NAA)[-1,] - simTT$NAA_devs) #same thing
mean(exp(simTT$NAA_devs[,-1])) #remove recruitment which has larger SD

simsTT <- sapply(1:100, function(x) modTT$simulate()$NAA_devs[,-1]) #remove recruitment which has larger SD
mean(exp(simsTT)) # = 1

simsFF <- sapply(1:100, function(x) modFF$simulate()$NAA_devs[,-1]) #remove recruitment which has larger SD
mean(exp(simsFF)) 
exp(0.5*exp(modFF$parList$log_NAA_sigma[2])^2) #equal to this

##########################################


