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
gamma.fun = function(data){
  mydata = data[data$delta.R == 1,]
  Y.R = mydata$Y.R
  Y.R_0 = log(Y.R)
  Y.R_2 = Y.R^2
  Y.R.mat = data.frame(cbind(Y.R_0,Y.R,Y.R_2, mydata[c("Delta.R", "X1.R", "X2.R")]))
  
  mod = glm((1-Delta.R) ~ Y.R_0 + Y.R + Y.R_2 + X1.R + X2.R,
            data = Y.R.mat, family = "binomial")
  
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
  # odat = odat[is.finite(rowSums(odat)),] # remove the rows with +/-Inf and rows containing NA
  
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
osmcf <- function(Nc, No, p, r, mu, lambda, nu, theta, beta, alpha, par1, ind.rho_est){
  #1) Generate clinical trial dataset
  cdat = cdata_gen(Nc, p, r, mu, lambda, nu, theta, beta, alpha)
  
  #2) Estimate gamma's and their se
  mod_gamma = gamma.fun(data = cdat)
  gammas.hat = mod_gamma$coeff
  gammas.se = mod_gamma$coeff_se
  gammas.cov = mod_gamma$vcov
  # gamma1.hat = gammas.hat[1:4]
  # gamma1.se = gammas.se[1:4]
  # gamma2.hat = gammas.hat[-c(1:4)]
  # gamma2.se = mod_gamma$coeff_se[-c(1:4)]
  
  #3) Generate observational study dataset and order it by observed event time
  odat = odata_gen(No, p, r, mu, lambda, nu, theta, beta, alpha)
  odat = odat[order(odat$Y),] # order by event time
  events_id = which(odat$delta == 1) # observed events index
  # Y = odat$Y
  # X1 = odat$X1
  # X2 = odat$X2
  # delta = odat$delta
  
  # #4) Estimate rho0
  # Y.mat = cbind(rep(1,length(Y)), log(Y), Y, Y^2)
  # rho0.hat = exp(as.vector(Y.mat%*%gamma1.hat))
  # #rho0.true = exp(log(nu) + nu*log(theta) - log(mu) - mu*log(lambda) + (nu-mu)*log(Y))
  # 
  # #5) Order the dataset by observed event time
  # odat = cbind(odat,rho0.hat)
  # odat = odat[is.finite(rowSums(odat)),] # remove the rows with +/-Inf and rows containing NA
  # odat = odat[order(odat$Y),] # order by event time
  # events_id = which(odat$delta == 1) # observed events index
  
  #4) Estimate beta's
  # par1 = c(0.1,0.1)
  beta.hat = nlm(f = neglog_pl_beta, par1, odat = odat, ind.rho_est, gammas.hat = gammas.hat, events_id = events_id)$estimate # Use nlm() to minimize negative log partial likelihood
  
  #5) Compute beta's se
  A.inv = mod_gamma$vcov*Nc # variance-covariance matrix of (gamma.hat - gamma0), A^-1 
  Cbeta.hat.sum = matrix(0, nrow = length(beta.hat), ncol = length(beta.hat))
  Cgamma.hat.sum = matrix(0, ncol = length(gammas.hat), nrow = length(beta.hat))
  # D.hat.sum = matrix(0, nrow = length(beta.hat), ncol = length(beta.hat))
  for (i in 1:length(events_id)){
    id = events_id[i]
    Yi = odat$Y[id] # Yi
    Xi = c(odat$X1[id], odat$X2[id]) # Xi
    risk.id = which(odat$Y >= Yi) # risk set R(Yi)
    Xj = cbind(odat$X1[risk.id], odat$X2[risk.id])
    Yi.tilde = cbind(rep(1,length(Yi)), log(Yi), Yi, Yi^2)
    nRmat = cbind(matrix(rep(Yi.tilde, length(risk.id)), nrow = length(risk.id), byrow = TRUE), Xj)
    
    phiij = c(exp(Xj%*%beta.hat) * (1+g.logit(nRmat%*%gammas.hat)))
    psiij = c(exp(Xj%*%beta.hat) * gd.logit(nRmat%*%gammas.hat)) * (nRmat)
    
    uprime = matrix(0, nrow = length(beta.hat), ncol = length(beta.hat))
    for (j in 1:length(risk.id)) {
      uprime = uprime + Xj[j,]%*%t(Xj[j,]) * phiij[j]
    }
    Cbeta.hat.sum = Cbeta.hat.sum + (uprime * sum(phiij) - colSums(Xj*phiij) %*% t(colSums(Xj*phiij))) / (sum(phiij))^2
    Cgamma.hat.sum = Cgamma.hat.sum + (t(Xj)%*%psiij * sum(phiij) - colSums(Xj*phiij) %*% t(colSums(psiij))) / (sum(phiij))^2
    Uibeta_gamma = Xi - colSums(Xj*phiij)/sum(phiij)
    # D.hat.sum = D.hat.sum + Uibeta_gamma %*% t(Uibeta_gamma)
  }
  Cbeta.hat = Cbeta.hat.sum/No
  Cgamma.hat = Cgamma.hat.sum/No
  # D.hat = D.hat.sum/No
  # Cbeta.hat
  # Cgamma.hat
  # D.hat
  
  beta.vcov = solve(Cbeta.hat)/No + solve(Cbeta.hat) %*% Cgamma.hat %*% (A.inv/Nc) %*% t(Cgamma.hat) %*% solve(Cbeta.hat)
  beta.se = sqrt(diag(beta.vcov))
  # beta.se
  
  coeff.est = c(gammas.hat, beta.hat, gammas.se, beta.se)
  return(list(coeff_est = coeff.est,
              gamma1_cov = gammas.cov[1:4,1:4]))
}

#Nc = 500000; No = 500; p = 0.5; r = 0.1; mu = 1.5; lambda = 0.5; beta = c(-0.5,1); alpha = c(-0.1,0.5); nu = 2; theta = 0.5; par1 = c(0.1,0.1)
