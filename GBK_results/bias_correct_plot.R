bias_correct_plot <- function(res, var, base.line = TRUE, relative = FALSE) {
  require(ggplot2)
  require(dplyr)
  res <- na.omit(res)
  res <- res[res$Convergence==1,]
  res$OM = factor(res$OM)
  res$EM = factor(res$EM)
  
  if (relative) res[,1:15] <- res[,1:15] - 1
  
  res1 = res[res$Match,]
  
  p <- ggplot(res, aes(x = EM, y = !! rlang::sym(var))) + 
    geom_boxplot() +
    facet_wrap(~ OM, ncol = 4,
               # scales = "free_y", 
               labeller = labeller(OM = c("1" = "Rec (iid)\npe1_oe1",
                                          "2" = "Rec (iid)\npe0_oe1",
                                          "3" = "Rec (iid)\npe1_oe0",
                                          "4" = "Rec (iid)\npe0_oe0",
                                          "5" = "Rec+1 (iid)\npe1_oe1",
                                          "6" = "Rec+1 (iid)\npe0_oe1",
                                          "7" = "Rec+1 (iid)\npe1_oe0",
                                          "8" = "Rec+1 (iid)\npe0_oe0",
                                          "9" = "Rec (ar1_y)\npe1_oe1",
                                          "10" = "Rec (ar1_y)\npe0_oe1",
                                          "11" = "Rec (ar1_y)\npe1_oe0",
                                          "12" = "Rec (ar1_y)\npe0_oe0",
                                          "13" = "Rec+1 (ar1_y)\npe1_oe1",
                                          "14" = "Rec+1 (ar1_y)\npe0_oe1",
                                          "15" = "Rec+1 (ar1_y)\npe1_oe0",
                                          "16" = "Rec+1 (ar1_y)\npe0_oe0",
                                          "17" = "Rec+1 (ar1_a)\npe1_oe1",
                                          "18" = "Rec+1 (ar1_a)\npe0_oe1",
                                          "19" = "Rec+1 (ar1_a)\npe1_oe0",
                                          "20" = "Rec+1 (ar1_a)\npe0_oe0",
                                          "21" = "Rec+1 (2dar1)\npe1_oe1",
                                          "22" = "Rec+1 (2dar1)\npe0_oe1",
                                          "23" = "Rec+1 (2dar1)\npe1_oe0",
                                          "24" = "Rec+1 (2dar1)\npe0_oe0"))) + 
    labs(x = "EM", y = "Mean (Sim/True)") + 
    ggtitle(var) + 
    scale_x_discrete(labels = c("1" = "pe1_oe1","2" = "pe0_oe1","3" = "pe1_oe0","4" = "pe0_oe0")) + 
    geom_boxplot(data = res1, fill = 'red') + 
    theme_bw()
  
  p + theme(axis.text.x = element_text(face = 2, angle = 45, hjust = 1),
            axis.text.y = element_text(face = 2),
            strip.text = element_text(face = 2))
  
  if (base.line) {
    if (relative) {
      p + geom_hline(yintercept=0, linetype="dashed", 
                     color = "blue", size=1) +
        labs(y = "Mean (Sim/True - 1)") 
    } else {
      p + geom_hline(yintercept=1, linetype="dashed", 
                     color = "blue", size=1)
    }
  } else {
    p + labs(y = "Estimated Value") 
  }
  
  ggsave(paste0(var,".PNG"),width = 10,height = 10)
}