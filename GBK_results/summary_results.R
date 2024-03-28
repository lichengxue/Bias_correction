setwd("E:/GBK")
# generating OM with pe and oe on and off (should have 4 different mods)
bias_correct_process     <- c(1,0)
bias_correct_observation <- c(1,0)
sigma <- c("rec","rec+1")
cor <- c("iid","ar1_y")

df.mods1 <- expand.grid(bias_correct_process = bias_correct_process, 
                        bias_correct_observation = bias_correct_observation, 
                        sigma = sigma,
                        cor = cor,
                        stringsAsFactors = FALSE)

bias_correct_process     <- c(1,0)
bias_correct_observation <- c(1,0)
sigma <- c("rec+1")
cor <- c("ar1_a","2dar1")

df.mods2 <- expand.grid(bias_correct_process = bias_correct_process, 
                        bias_correct_observation = bias_correct_observation, 
                        sigma = sigma,
                        cor = cor,
                        stringsAsFactors = FALSE) 

df.mods <- rbind(df.mods1,df.mods2)

str <- df.mods[1:4,1:2]

df = readRDS("dat.RDS")
res = NULL
for (om in 1:24) {
  num = df[[om]]$ID 
  files_names <- paste0("GBK_full/OM",om,"/OM",om,"nsim",num,".RDS")
  mods_list <- sapply(files_names, try(readRDS), simplify = FALSE)
  em = as.numeric(rownames(str[om %% 4,])) # true model is j
  if (length(em) == 0) em = 4
  match = rep(FALSE,4)
  match[em] = TRUE
  # mod.tmp <- mods_list[[i]][[j]]
  for (i in 1:50) {
    mean_naa <- sapply(1:4,function(j) {
      mod.tmp <- mods_list[[i]][[j]]
      if(length(mod.tmp) != 0){
        sapply(1:6, function(a) {     
          mean.NAA <- mean(mod.tmp$rep$NAA[-1,a]/exp(mod.tmp$input$data$log_NAA[,a]))
        })
      } else {
        mean.NAA <- matrix(NA,6,1)
      }
    })
    
    # total
    mean_naa_total <- sapply(1:4,function(j) {
      mod.tmp <- mods_list[[i]][[j]]
      if(length(mod.tmp) != 0){
        mean_naa_total <- mean(mod.tmp$rep$NAA[-1,-1]/exp(mod.tmp$input$data$log_NAA[,-1]))
      } else {
        mean_naa_total <- NA
      }
    })
    
    #
    mean_catch <- sapply(1:4,function(j) {
      mod.tmp <- mods_list[[i]][[j]]
      if(length(mod.tmp) != 0){
        mean_catch <- mean(mod.tmp$rep$pred_catch/mod.tmp$input$data$agg_catch)
      } else {
        mean_catch <- NA
      }
    })
    
    Fmax <- sapply(1:4,function(j) {
      mod.tmp <- mods_list[[i]][[j]]
      if(length(mod.tmp) != 0){
        Fmax <- mean(mod.tmp$rep$F/mod.tmp$input$data$F)
      } else {
        Fmax <- NA
      }
    })
    
    mean_Index <- sapply(1:4,function(j) {
      mod.tmp <- mods_list[[i]][[j]]
      if(length(mod.tmp) != 0){
        sapply(1:3, function(a) {     
          mean.Index <- mean(mod.tmp$rep$pred_indices[,a]/mod.tmp$input$data$agg_indices[,a])
        })
      } else {
        mean.Index <- matrix(NA,3,1)
      }
    })
    
    mean_q <- sapply(1:4,function(j) {
      mod.tmp <- mods_list[[i]][[j]]
      if(length(mod.tmp) != 0){
        sapply(1:3, function(a) {     
          mean.q <- mean(mod.tmp$rep$q[,a]/mod.tmp$input$data$q[,a])
        })
      } else {
        mean.q <- matrix(NA,3,1)
      }
    })
    
    NAA_sigma <- sapply(1:4,function(j) {
      mod.tmp <- mods_list[[i]][[j]]
      if(length(mod.tmp) != 0){
        if (length(mod.tmp$parList$log_NAA_sigma) == 2) {
          NAA.sigma <- exp(mod.tmp$parList$log_NAA_sigma)
        } else {
          NAA.sigma <- exp(c(mod.tmp$parList$log_NAA_sigma,-Inf))
        }
      } else {
        NAA.sigma <- matrix(NA,2,1)
      }
    })
    
    NAA_rho <- sapply(1:4,function(j) {
      mod.tmp <- mods_list[[i]][[j]]
      if(length(mod.tmp) != 0){
        NAA.rho <- mod.tmp$parList$trans_NAA_rho
        NAA.rho <- 2/(1+exp(-2*NAA.rho)) - 1
      } else {
        NAA.rho <- matrix(NA,2,1)
      }
    })
    
    selfwt_par <- sapply(1:4,function(j) {
      mod.tmp <- mods_list[[i]][[j]]
      if(length(mod.tmp) != 0){
        tmp1 <- exp(as.list(mod.tmp$sdrep,"Est")$catch_paa_pars[,1])
        tmp2 <- exp(as.list(mod.tmp$sdrep,"Est")$index_paa_pars[,1])
        selfwt_par <- c(tmp1,tmp2)
      } else {
        selfwt_par <- matrix(NA,4,1)
      }
    })
    
    convergence <- sapply(1:4,function(j) {
      mod.tmp <- mods_list[[i]][[j]]
      if(length(mod.tmp) != 0){
        conv = as.logical(1-mod.tmp$opt$convergence)
        pdHess = as.logical(if(mod.tmp$na_sdrep==FALSE & !is.na(mod.tmp$na_sdrep)) 1 else 0)
        if (!conv | !pdHess) {
          convergence = 0
        } else {
          convergence = 1
        }
      } else {
        convergence <- NA
      }
    })
    
    AIC <- sapply(1:4,function(j) {
      mod.tmp <- mods_list[[i]][[j]]
      if(length(mod.tmp) != 0){
        aic <- 2*(mod.tmp$opt$obj + length(mod.tmp$opt$par)) 
      } else {
        aic <- NA
      }
    })
    
    dAIC = AIC-min(AIC)
    
    dat = data.frame(t(mean_naa),mean_naa_total, mean_catch, Fmax, t(mean_Index), t(mean_q), t(NAA_sigma),t(NAA_rho), 
                     t(selfwt_par), convergence, AIC, dAIC, OM = om, EM = 1:4, Match = match, nsim = i)
    res = rbind(res,dat)
  }
}

names(res)[1:26] <- c(paste0("Age_",1:6),"NAA","Catch", "Fmax", paste0("Index_",1:3),paste0("Catchability_",1:3),
                      "Rec_sigma","NAA_sigma","Rho_a","Rho_y",paste0("selfwt_par",1:4),"Convergence","AIC","dAIC")
write.csv(res,"res1.csv",row.names = F)

# res <- NULL
# for (om in 1:24) {
#   files_names <- paste0("GBK_full/OM",om,"nsim",1:100,".RDS")
#   mods_list <- sapply(files_names, try(readRDS), simplify = FALSE)
#   em = as.numeric(rownames(str[om %% 4,])) # true model is j
#   if (length(em) == 0) em = 4
#   match = rep(FALSE,4)
#   match[em] = TRUE
#   # mod.tmp <- mods_list[[i]][[j]]
#   for (i in 1:100) {
#     mean_naa <- sapply(1:4,function(j) {
#       mod.tmp <- mods_list[[i]][[j]]
#       if(length(mod.tmp) != 0){
#         sapply(1:6, function(a) {     
#           mean.NAA <- mean(mod.tmp$rep$NAA[-1,a]/exp(mod.tmp$input$data$log_NAA[,a]))
#         })
#       } else {
#         mean.NAA <- matrix(NA,6,1)
#       }
#     })
#     
#     # total
#     mean_naa_total <- sapply(1:4,function(j) {
#       mod.tmp <- mods_list[[i]][[j]]
#       if(length(mod.tmp) != 0){
#         mean_naa_total <- mean(mod.tmp$rep$NAA[-1,-1]/exp(mod.tmp$input$data$log_NAA[,-1]))
#       } else {
#         mean_naa_total <- NA
#       }
#     })
#     
#     #
#     mean_catch <- sapply(1:4,function(j) {
#       mod.tmp <- mods_list[[i]][[j]]
#       if(length(mod.tmp) != 0){
#         mean_catch <- mean(mod.tmp$rep$pred_catch/mod.tmp$input$data$agg_catch)
#       } else {
#         mean_catch <- NA
#       }
#     })
#     
#     Fmax <- sapply(1:4,function(j) {
#       mod.tmp <- mods_list[[i]][[j]]
#       if(length(mod.tmp) != 0){
#         Fmax <- mean(mod.tmp$rep$F/mod.tmp$input$data$F)
#       } else {
#         Fmax <- NA
#       }
#     })
#     
#     mean_Index <- sapply(1:4,function(j) {
#       mod.tmp <- mods_list[[i]][[j]]
#       if(length(mod.tmp) != 0){
#         sapply(1:3, function(a) {     
#           mean.Index <- mean(mod.tmp$rep$pred_indices[,a]/mod.tmp$input$data$agg_indices[,a])
#         })
#       } else {
#         mean.Index <- matrix(NA,3,1)
#       }
#     })
#     
#     mean_q <- sapply(1:4,function(j) {
#       mod.tmp <- mods_list[[i]][[j]]
#       if(length(mod.tmp) != 0){
#         sapply(1:3, function(a) {     
#           mean.q <- mean(mod.tmp$rep$q[,a]/mod.tmp$input$data$q[,a])
#         })
#       } else {
#         mean.q <- matrix(NA,3,1)
#       }
#     })
#     
#     NAA_sigma <- sapply(1:4,function(j) {
#       mod.tmp <- mods_list[[i]][[j]]
#       if(length(mod.tmp) != 0){
#         if (length(mod.tmp$parList$log_NAA_sigma) == 2) {
#           NAA.sigma <- exp(mod.tmp$parList$log_NAA_sigma)
#         } else {
#           NAA.sigma <- exp(c(mod.tmp$parList$log_NAA_sigma,-Inf))
#         }
#       } else {
#         NAA.sigma <- matrix(NA,2,1)
#       }
#     })
#     
#     NAA_rho <- sapply(1:4,function(j) {
#       mod.tmp <- mods_list[[i]][[j]]
#       if(length(mod.tmp) != 0){
#         NAA.rho <- mod.tmp$parList$trans_NAA_rho
#         NAA.rho <- 2/(1+exp(-2*NAA.rho)) - 1
#       } else {
#         NAA.rho <- matrix(NA,2,1)
#       }
#     })
#     
#     convergence <- sapply(1:4,function(j) {
#       mod.tmp <- mods_list[[i]][[j]]
#       if(length(mod.tmp) != 0){
#         conv = as.logical(1-mod.tmp$opt$convergence)
#         pdHess = as.logical(if(mod.tmp$na_sdrep==FALSE & !is.na(mod.tmp$na_sdrep)) 1 else 0)
#         if (!conv | !pdHess) {
#           convergence = 0
#         } else {
#           convergence = 1
#         }
#       } else {
#         convergence <- NA
#       }
#     })
#     
#     AIC <- sapply(1:4,function(j) {
#       mod.tmp <- mods_list[[i]][[j]]
#       if(length(mod.tmp) != 0){
#         aic <- 2*(mod.tmp$opt$obj + length(mod.tmp$opt$par)) 
#       } else {
#         aic <- NA
#       }
#     })
#     
#     dat = data.frame(t(mean_naa),mean_naa_total, mean_catch, Fmax, t(mean_Index), t(mean_q), t(NAA_sigma),t(NAA_rho), convergence, AIC, OM = om, EM = 1:4, Match = match, nsim = i)
#     res = rbind(res,dat)
#   }
# }
# 
# names(res)[1:21] <- c(paste0("Age_",1:6),"NAA","Catch", "Fmax", paste0("Index_",1:3),paste0("Catchability_",1:3),"Rec_sigma","NAA_sigma",
#                      "Rho_a","Rho_y","Convergence","AIC")
# write.csv(res,"res1.csv",row.names = F)

require(ggplot2)
require(dplyr)
source("bias_correct_plot.R")

names = names(res)[1:26]

for (i in 1:15){
  var = names[i]
  bias_correct_plot(res,var,base.line = TRUE, relative = FALSE)
}

for (i in 1:15){
  var = names[i]
  bias_correct_plot(res,var,base.line = TRUE, relative = TRUE)
}

for (i in 16:26){
  var = names[i]
  bias_correct_plot(res,var,base.line = FALSE)
}


