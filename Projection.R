lib.dir = "C:/Program Files/R/R-4.3.2/library/single_wham"
library(wham,lib.loc = lib.dir)
library(dplyr)
dir <- "C:/Users/spaswham/Desktop/Bias_correction"
setwd(dir)
Parallel = TRUE
if (Parallel) {
  library(doParallel)
  myCluster <- makeCluster(60, # number of cores to use
                           type = "PSOCK") # type of cluster
  print(paste0("Number of Cores: ",detectCores()))
  print(paste0("Cores Registered: ",myCluster))
  registerDoParallel(myCluster)
}
for(m in 1:4) {
  for (n in 1:8) {
    foreach (k = 1:100) %dopar% {
      lib.dir = "C:/Program Files/R/R-4.3.2/library/single_wham"
      library(wham,lib.loc = lib.dir)
      mod <- readRDS(paste0("Res_",'OM_',m,'_EM_',n,"_nsim_",k,".RDS"))
      proj <- project_wham(mod, proj.opts=list(n.yrs=50, use.FXSPR=TRUE, proj_R_opt = 1), do.sdrep=T)
      saveRDS(proj$rep,paste0("Proj_",'OM_',m,'_EM_',n,"_nsim_",k,".RDS"))
    }
  }
}
