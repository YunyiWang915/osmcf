################ 1. Data generating function ################
# 1) function for clinical trial dataset
cdata_gen <- function(Nc, p, r, mu, lambda, nu, theta, beta, alpha){
  C.R = rexp(Nc, rate = r) # censoring time C ~ Exp(r)
  
  X1.R = rbinom(Nc, size = 1, prob = p) # covariates X1 ~ Ber(q), X2 ~ U(0,1)
  X2.R = runif(Nc)
  
  T.B = rweibull(Nc, shape = mu, scale = (1/lambda)*exp(-cbind(X1.R,X2.R)%*%beta/mu))  # T.B ~ Weibull(mu, lambda^mu*exp(beta'*X.R))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      )
  T.O = rweibull(Nc, shape = nu, scale = (1/theta)*exp(-cbind(X1.R,X2.R)%*%alpha/nu)) # T.O ~ Weibull(nu, theta^nu*exp(alpha'*X.R))
  
  T.R = pmin(T.B, T.O) # true event time
  Y.R = pmin(T.R, C.R) # observed event time
  
  Delta.R = rep(NA, Nc)
  Delta.R[Y.R == T.B] = 1 # cause of failure indicator
  Delta.R[Y.R == T.O] = 0
  delta.R = (T.R <= C.R)*1 # event indicator
  
  # data frame
  Y.R = Y.R/0.6 # observed event time normalized by the standard deviation of Y based on a large dataset
  obs_data = as.data.frame(cbind(Y.R, Delta.R, delta.R, X1.R, X2.R))  
  obs_data$id = seq(1,Nc,1)
  return(obs_data)
}

# 2) function for observational study dataset
odata_gen<-function(No, p, r, mu, lambda, nu, theta, beta, alpha){
  C = rexp(No, rate = r) # censoring time C ~ Exp(r)
  
  X1 = rbinom(No, size = 1, prob = p) # covariates X1 ~ Ber(q), X2 ~ U(0,1)
  X2 = runif(No)
  
  T.B = rweibull(No, shape = mu, scale = (1/lambda)*exp(-cbind(X1,X2)%*%beta/mu))  # T.B ~ Weibull(mu, lambda^mu*exp(beta'*X))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       )
  T.O = rweibull(No, shape = nu, scale = (1/theta)*exp(-cbind(X1,X2)%*%alpha/nu)) # T.O ~ Weibull(nu, theta^nu*exp(alpha'*X))
  
  Tt = pmin(T.B, T.O) # true event time
  Y = pmin(Tt, C) # observed event time
  
  Delta = rep(NA, No)
  Delta[Y == T.B] = 1 # cause of failure indicator
  Delta[Y == T.O] = 0
  delta = (Tt <= C)*1 # event indicator
  
  # data frame
  Y = Y/0.6 # observed event time normalized by the standard deviation of Y based on a large dataset
  obs_data = as.data.frame(cbind(Y, delta, X1, X2))  
  obs_data$id = seq(1,No,1)
  return(obs_data)
}


################ 2. Estimation function for gamma ################
gamma.fun = function(data, ind.wt){
  mydata = data[data$delta.R == 1,]
  Y.R = mydata$Y.R
  Y.R_0 = log(Y.R)
  Y.R_2 = Y.R^2
  
  if (ind.wt == 0){
    Y.R.mat = data.frame(cbind(Y.R_0,Y.R,Y.R_2, mydata[c("Delta.R", "X1.R", "X2.R")]))
    mod = glm((1-Delta.R) ~ Y.R_0 + Y.R + Y.R_2 + X1.R + X2.R,
              data = Y.R.mat, family = "binomial")
  }
  
  if (ind.wt == 1){
    Y.R.mat = data.frame(cbind(Y.R_0,Y.R,Y.R_2, mydata[c("Delta.R", "X1.R", "X2.R","rwt")]))
    mod = glm((1-Delta.R) ~ Y.R_0 + Y.R + Y.R_2 + X1.R + X2.R,
              data = Y.R.mat, family = "binomial", weights = rwt)
  }
  
  coeff = coef(summary(mod))[,1]
  coeff_se = coef(summary(mod))[,2]
  
  vcov = vcov(mod)
  
  return(list(coeff = coeff, coeff_se = coeff_se, vcov = vcov))
}


################ 3. Estimation function for beta ################
# negative log partial likelihood function
neglog_pl_beta = function(odat, beta, ind.rho_est, gammas.hat, events_id){
  Y = odat$Y
  Y.mat = cbind(rep(1,length(Y)), log(Y), Y, Y^2)
  
  gammas.true = c(-0.314, 0.5, 0, 0, 0.4, -0.5) # true value
  if (ind.rho_est == 1){gamma1 = gammas.hat[1:4]; gamma2 = gammas.hat[-c(1:4)]}
  if (ind.rho_est == 0){gamma1 = gammas.true[1:4]; gamma2 = gammas.true[-c(1:4)]}
  
  # Estimate rho0
  rho0.hat = exp(as.vector(Y.mat%*%gamma1))
  odat = cbind(odat,rho0.hat)
  
  log_parlike = 0
  for (l in 1:length(events_id)){
    id = events_id[l]
    Yi = odat$Y[id] # Yi
    Xi = c(odat$X1[id], odat$X2[id]) # Xi
    risk.id = which(odat$Y >= Yi) # risk set R(Yi)
    Xj = cbind(odat$X1[risk.id], odat$X2[risk.id])
    
    rho0.hat = odat$rho0.hat[id]
    
    log_parlike = log_parlike - Xi%*%beta - log(1+rho0.hat*exp(Xi%*%gamma2)) + log(sum(exp(Xj%*%beta)*(1+rho0.hat*exp(Xj%*%gamma2))))
  }
  return(log_parlike)
}


################ 4. Logit link function ################
g.logit <- function(xx){exp(xx) / (1 + exp(xx))}
gd.logit <- function(xx){exp(xx) / (1+exp(xx))^2} # derivative of logit link function


################ 5. Main function for estimating coefficients and their se ################
osmcf.pr <- function(Nc, No, p, r, mu, lambda, nu, theta, beta, alpha, par1, ind.rho_est){
  #1) Generate clinical trial dataset
  cdat = cdata_gen(Nc, p, r, mu, lambda, nu, theta, beta, alpha)
  
  #2) Estimate gamma's and their se
  mod_gamma = gamma.fun(data = cdat, ind.wt = 0)
  gammas.hat = mod_gamma$coeff
  gammas.se = mod_gamma$coeff_se
  gammas.cov = mod_gamma$vcov
  
  #3) Generate observational study dataset and order it by observed event time
  odat = odata_gen(No, p, r, mu, lambda, nu, theta, beta, alpha)
  odat = odat[order(odat$Y),] # order by event time
  events_id = which(odat$delta == 1) # observed events index
  
  #4) Estimate beta's
  # par1 = c(0.1,0.1)
  beta.hat = nlm(f = neglog_pl_beta, par1, odat = odat, ind.rho_est, gammas.hat = gammas.hat, events_id = events_id)$estimate # Use nlm() to minimize negative log partial likelihood
  
  #5) Compute beta's se using perturbation
  beta.hat.pr = data.frame(matrix(nrow = 200, ncol = 2)) 
  for (k in 1:200) {
    cdat.pr = cdat
    cdat.pr$rwt = rexp(Nc,1) # generating random weights using standard exponential distribution
    mod_gamma.pr = gamma.fun(data = cdat.pr, ind.wt = 1)
    gammas.hat.pr = mod_gamma.pr$coeff
    
    odat.pr = odat
    odat.pr$rwt = rexp(No,1)
    odat.pr = odat.pr[order(odat.pr$Y),]
    events_id.pr = which(odat.pr$delta == 1)
    
    beta.hat.pr[k,] = nlm(f = neglog_pl_beta, par1, odat = odat.pr, ind.rho_est, gammas.hat = gammas.hat.pr, events_id = events_id.pr)$estimate
  }
  
  beta.se = apply(beta.hat.pr,2,sd)
  coeff.est = c(gammas.hat, beta.hat, gammas.se, beta.se)
  return(list(coeff_est = coeff.est,
              gamma1_cov = gammas.cov[1:4,1:4]))
}
