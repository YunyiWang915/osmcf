osmcf_res_sum <- function(est.res, vcov.res, Y = seq(0.01, 3, 0.01),
                          mu, lambda, nu, theta, beta, alpha, nsim = 1000){
  colnames(est.res) = c("gamma1.int", "gamma1.0",
                        "gamma1.1", "gamma1.2",
                        "gamma2.X1", "gamma2.X2",
                        "beta.X1", "beta.X2",
                        
                        "gamma1.int.se", "gamma1.0.se",
                        "gamma1.1.se", "gamma1.2.se",
                        "gamma2.X1.se", "gamma2.X2.se",
                        "beta.X1.se", "beta.X2.se")
  
  ############# rho0.hat, rho0.true, rho0.se #################
  Y.std = Y/0.6 # time points should be standardized
  Y.tilde = cbind(rep(1,length(Y)), log(Y.std), Y.std, Y.std^2)
  
  rho0.true = exp(log(nu) + nu*log(theta) - log(mu) - mu*log(lambda) + (nu-mu)*log(Y)) 
  
  rho0.hat.df = as.data.frame(matrix(nrow = nsim, ncol = length(Y)))
  rho0.se.df = as.data.frame(matrix(nrow = nsim, ncol = length(Y)))
  rho0.cil.df = as.data.frame(matrix(nrow = nsim, ncol = length(Y)))
  rho0.ciu.df = as.data.frame(matrix(nrow = nsim, ncol = length(Y)))
  for (i in 1:nsim) {
    rho0.hat.df[i,] = c(exp(Y.tilde%*%t(est.res[i,1:4])))
    rho0.se.df[i,] = sqrt(diag(Y.tilde %*% vcov.res[,,i] %*% t(Y.tilde)))
    rho0.cil.df[i,] = rho0.hat.df[i,] + qnorm(.025) * rho0.se.df[i,]
    rho0.ciu.df[i,] = rho0.hat.df[i,] + qnorm(.975) * rho0.se.df[i,]
  }
  
  rho0.mean = round(colMeans(rho0.hat.df),3)
  rho0.sd = round(apply(rho0.hat.df,2,sd),3)
  rho0.se = round(colMeans(rho0.se.df),3)
  rho0.cil = round(colMeans(rho0.cil.df),3)
  rho0.ciu = round(colMeans(rho0.ciu.df),3)
  rho0.res = as.data.frame(cbind(Y, rho0.true, rho0.mean, rho0.sd, rho0.se, rho0.cil, rho0.ciu))
  #mean((rho0.mean - rho0.true)^2) # MSE between true and estimated rho0
  
  ############# coefficients coverage probabilities #################
  gamma1.int_cil = est.res$gamma1.int + qnorm(.025) * est.res$gamma1.int.se
  gamma1.int_ciu = est.res$gamma1.int + qnorm(.975) * est.res$gamma1.int.se
  gamma1.0_cil = est.res$gamma1.0 + qnorm(.025) * est.res$gamma1.0.se
  gamma1.0_ciu = est.res$gamma1.0 + qnorm(.975) * est.res$gamma1.0.se
  gamma1.1_cil = est.res$gamma1.1 + qnorm(.025) * est.res$gamma1.1.se
  gamma1.1_ciu = est.res$gamma1.1 + qnorm(.975) * est.res$gamma1.1.se
  gamma1.2_cil = est.res$gamma1.2 + qnorm(.025) * est.res$gamma1.2.se
  gamma1.2_ciu = est.res$gamma1.2 + qnorm(.975) * est.res$gamma1.2.se
  
  gamma2.X1_cil = est.res$gamma2.X1 + qnorm(.025) * est.res$gamma2.X1.se
  gamma2.X1_ciu = est.res$gamma2.X1 + qnorm(.975) * est.res$gamma2.X1.se
  gamma2.X2_cil = est.res$gamma2.X2 + qnorm(.025) * est.res$gamma2.X2.se
  gamma2.X2_ciu = est.res$gamma2.X2 + qnorm(.975) * est.res$gamma2.X2.se
  
  beta.X1_cil = est.res$beta.X1 + qnorm(.025) * est.res$beta.X1.se
  beta.X1_ciu = est.res$beta.X1 + qnorm(.975) * est.res$beta.X1.se
  beta.X2_cil = est.res$beta.X2 + qnorm(.025) * est.res$beta.X2.se
  beta.X2_ciu = est.res$beta.X2 + qnorm(.975) * est.res$beta.X2.se
  
  
  gamma10 = c(log(nu) + nu*log(theta) - log(mu) - mu*log(lambda) + (nu-mu)*log(0.6), nu-mu, 0, 0) # true value
  gamma20 = alpha - beta # true value
  beta0 = beta # beta true value
  trues = round(c(gamma10, gamma20, beta0),3)
  
  
  gamma1.int.cp = (gamma1.int_cil < gamma10[1]) & (gamma10[1] < gamma1.int_ciu)
  gamma1.0.cp = (gamma1.0_cil < gamma10[2]) & (gamma10[2] < gamma1.0_ciu)
  gamma1.1.cp = (gamma1.1_cil < gamma10[3]) & (gamma10[3] < gamma1.1_ciu)
  gamma1.2.cp = (gamma1.2_cil < gamma10[4]) & (gamma10[4] < gamma1.2_ciu)
  gamma2.X1.cp = (gamma2.X1_cil < gamma20[1]) & (gamma20[1] < gamma2.X1_ciu)
  gamma2.X2.cp = (gamma2.X2_cil < gamma20[2]) & (gamma20[2] < gamma2.X2_ciu)
  beta.X1.cp = (beta.X1_cil < beta0[1]) & (beta0[1] < beta.X1_ciu)
  beta.X2.cp = (beta.X2_cil < beta0[2]) & (beta0[2] < beta.X2_ciu)
  est.res = cbind(est.res, gamma1.int.cp, gamma1.0.cp, gamma1.1.cp, gamma1.2.cp,
                  gamma2.X1.cp, gamma2.X2.cp, beta.X1.cp, beta.X2.cp)
  
  
  coeff.res = data.frame("Parameter" = c("gamma1.int", "gamma1.0", "gamma1.1", "gamma1.2",
                                          "gamma2.X1", "gamma2.X2", 
                                          "beta.X1", "beta.X2"),
                                  "True Value" = trues,
                                  "Estimates" = round(colMeans(est.res[,1:8]),3),
                                  "Bias" = round(colMeans(est.res[,1:8]),3) - trues,
                                  "SD" = round(apply(est.res[,1:8],2,sd),3),
                                  "SE" = round(colMeans(est.res[,9:16]),3),
                                  "CP" = round(colMeans(est.res[,17:24]),3))
  return(list(rho0.res = rho0.res, coeff.res = coeff.res))
}
