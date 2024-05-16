osmcf_plot_rho0 <- function(rho0.res){
  library(ggplot2)
  
  ####################### plot for rho0 #######################
  xlim = round(quantile(rho0.res$Y, c(0.05,0.95)),2)
  plot_rho0 = ggplot(rho0.res, aes(x = Y)) + 
    geom_line(aes(y = rho0.true, color = "True"), size = 0.8) +
    geom_line(aes(y = rho0.mean, color = "Fractional Polynomial"), size = 0.8) +
    geom_line(aes(y = rho0.cil, color = "Fractional Polynomial"), linetype = "dotted", size = 0.8) +
    geom_line(aes(y = rho0.ciu, color = "Fractional Polynomial"), linetype = "dotted", size = 0.8) +
    scale_x_continuous("t", limits = c(xlim[1],xlim[2])) +
    scale_y_continuous("", limits = c(-0.5,3)) +
    scale_color_manual("",
                       breaks = c("True", "Fractional Polynomial"),
                       values = c("True" = "black", "Fractional Polynomial" = "#F8766D"),
                       drop = FALSE) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.position = c(.4, 1), legend.justification = c("right", "top"), legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6))
  return(plot_rho0)
}