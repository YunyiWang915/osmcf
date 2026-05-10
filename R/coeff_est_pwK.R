# -----------------------------
# Piecewise-constant PH time generator
# -----------------------------
#' Simulate proportional hazards event times with a piecewise-constant baseline hazard
#'
#' This internal function generates event times from a proportional hazards model
#' with a piecewise-constant baseline hazard using inverse transform sampling.
#'
#' @param n Number of subjects.
#' @param eta Numeric vector of linear predictors. Its length must be equal to `n`.
#' @param cuts Numeric vector of increasing cut-points defining the finite intervals.
#' @param lam Numeric vector of baseline hazard values. Its length must be
#'   `length(cuts) + 1`.
#'
#' @return A numeric vector of simulated event times.
#'
#' @keywords internal

rph_piecewise_UK <- function(n, eta, cuts, lam) {
  stopifnot(length(eta) == n)
  stopifnot(length(lam) == length(cuts) + 1)
  stopifnot(all(diff(cuts) > 0), all(lam > 0))

  U <- runif(n)
  E <- -log(U)
  z <- E / exp(eta)  # target baseline cumulative hazard H0(T)

  dt <- diff(c(0, cuts))                      # lengths of finite intervals
  H_end <- cumsum(lam[1:length(cuts)] * dt)   # cum baseline hazard at each cutoff

  T <- numeric(n)

  for (i in seq_len(n)) {
    zi <- z[i]
    k <- which(zi <= H_end)[1]

    if (is.na(k)) {
      # last interval
      H_prev <- if (length(H_end) == 0) 0 else H_end[length(H_end)]
      t_start <- if (length(cuts) == 0) 0 else cuts[length(cuts)]
      T[i] <- t_start + (zi - H_prev) / lam[length(lam)]
    } else {
      H_prev <- if (k == 1) 0 else H_end[k - 1]
      t_start <- if (k == 1) 0 else cuts[k - 1]
      T[i] <- t_start + (zi - H_prev) / lam[k]
    }
  }
  T
}



################ 1. Data generating function ################
#' Generate a clinical trial dataset
#'
#' Generates a clinical trial dataset with observed cause-of-death information
#' under piecewise-constant cause-specific proportional hazards models.
#'
#' @param Nc Number of subjects in the clinical trial dataset.
#' @param p Probability parameter for generating the binary covariate.
#' @param r Rate parameter for the exponential censoring distribution.
#' @param cutsB Numeric vector of cut-points for the disease-specific baseline hazard.
#' @param lamB Numeric vector of disease-specific baseline hazard values.
#' @param cutsO Numeric vector of cut-points for the other-cause baseline hazard.
#' @param lamO Numeric vector of other-cause baseline hazard values.
#' @param beta Numeric vector of regression coefficients for disease-specific death.
#' @param alpha Numeric vector of regression coefficients for other-cause death.
#'
#' @return A data frame with variables `Y.R`, `Delta.R`, `delta.R`, `X1.R`,
#'   `X2.R`, and `id`.
#'
#' @export

cdata_genK <- function(Nc, p, r,
                       cutsB, lamB, cutsO, lamO,
                       beta, alpha) {

  C.R  <- rexp(Nc, rate = r)
  X1.R <- rbinom(Nc, 1, p)
  X2.R <- runif(Nc)

  etaB <- as.numeric(cbind(X1.R, X2.R) %*% beta)
  etaO <- as.numeric(cbind(X1.R, X2.R) %*% alpha)

  T.B <- rph_piecewise_UK(Nc, etaB, cutsB, lamB)
  T.O <- rph_piecewise_UK(Nc, etaO, cutsO, lamO)

  T.R <- pmin(T.B, T.O)
  Y.R <- pmin(T.R, C.R)

  delta.R <- as.integer(T.R <= C.R)

  Delta.R <- rep(NA_integer_, Nc)
  Delta.R[delta.R == 1 & T.B <= T.O] <- 1L
  Delta.R[delta.R == 1 & T.O <  T.B] <- 0L

  data.frame(Y.R, Delta.R, delta.R, X1.R, X2.R, id = seq_len(Nc))
}


#' Generate an observational study dataset
#'
#' Generates an observational study dataset under piecewise-constant
#' cause-specific proportional hazards models. The latent cause-of-death
#' indicator is included for simulation evaluation.
#'
#' @param No Number of subjects in the observational study dataset.
#' @param p Probability parameter for generating the binary covariate.
#' @param r Rate parameter for the exponential censoring distribution.
#' @param cutsB Numeric vector of cut-points for the disease-specific baseline hazard.
#' @param lamB Numeric vector of disease-specific baseline hazard values.
#' @param cutsO Numeric vector of cut-points for the other-cause baseline hazard.
#' @param lamO Numeric vector of other-cause baseline hazard values.
#' @param beta Numeric vector of regression coefficients for disease-specific death.
#' @param alpha Numeric vector of regression coefficients for other-cause death.
#'
#' @return A data frame with variables `Y`, `Delta`, `delta`, `X1`, `X2`, and `id`.
#'
#' @export

odata_genK <- function(No, p, r,
                       cutsB, lamB, cutsO, lamO,
                       beta, alpha) {

  C  <- rexp(No, rate = r)
  X1 <- rbinom(No, 1, p)
  X2 <- runif(No)

  etaB <- as.numeric(cbind(X1, X2) %*% beta)
  etaO <- as.numeric(cbind(X1, X2) %*% alpha)

  T.B <- rph_piecewise_UK(No, etaB, cutsB, lamB)
  T.O <- rph_piecewise_UK(No, etaO, cutsO, lamO)

  Tt <- pmin(T.B, T.O)
  Y  <- pmin(Tt, C)

  delta <- as.integer(Tt <= C)

  Delta <- rep(NA_integer_, No)
  Delta[delta == 1 & T.B <= T.O] <- 1L
  Delta[delta == 1 & T.O <  T.B] <- 0L

  data.frame(Y, Delta, delta, X1, X2, id = seq_len(No))
}


################ 2. Estimation function for gamma ################
#' Estimate proportional cause ratio model parameters
#'
#' Estimates the parameters in the proportional cause ratio model using clinical
#' trial data with observed cause-of-death information.
#'
#' @param data Clinical trial dataset. The data should contain `Y.R`, `Delta.R`,
#'   `delta.R`, `X1.R`, and `X2.R`. If `ind.wt = 1`, the data should also contain
#'   a perturbation weight variable named `rwt`.
#' @param ind.wt Indicator for perturbation weighting. Use `0` for the original
#'   estimation and `1` for weighted estimation.
#'
#' @return A list with components:
#' \describe{
#'   \item{coeff}{Estimated proportional cause ratio model coefficients.}
#'   \item{coeff_se}{Estimated standard errors of the coefficients.}
#'   \item{vcov}{Estimated variance-covariance matrix.}
#' }
#'
#' @export


gamma_fun = function(data, ind.wt){
  mydata = data[data$delta.R == 1,]
  Y.R = mydata$Y.R
  Y.R_0 = log(Y.R)
  Y.R_2 = Y.R^2

  if (ind.wt == 0){
    gamma.mat = data.frame(cbind(Y.R_0,Y.R,Y.R_2, mydata[c("Delta.R", "X1.R", "X2.R")]))
    mod = glm((1-Delta.R) ~ Y.R_0 + Y.R + Y.R_2 + X1.R + X2.R,
              data = gamma.mat, family = "binomial")
  }
  if (ind.wt == 1){
    gamma.mat = data.frame(cbind(Y.R_0,Y.R,Y.R_2, mydata[c("Delta.R", "X1.R", "X2.R","rwt")]))
    mod = suppressWarnings(
      glm((1-Delta.R) ~ Y.R_0 + Y.R + Y.R_2 + X1.R + X2.R,
        data = gamma.mat,
        family = "binomial",
        weights = gamma.mat$rwt))
  }

  coeff = coef(summary(mod))[,1]
  coeff_se = coef(summary(mod))[,2]

  vcov = vcov(mod)

  return(list(coeff = coeff, coeff_se = coeff_se, vcov = vcov))
}


################ 3. Estimation function for beta ################
#' Negative log partial likelihood for beta estimation
#'
#' Internal objective function used in nonlinear optimization for estimating
#' disease-specific survival regression coefficients.
#'
#' @param odat Observational study dataset.
#' @param beta Numeric vector of disease-specific survival regression coefficients.
#' @param gammas.hat Estimated proportional cause ratio model parameters.
#' @param events_id Integer vector indicating rows corresponding to observed events.
#' @param ind.wt Indicator for perturbation weighting. Use `0` for the original
#'   estimation and `1` for weighted estimation.
#'
#' @return A numeric value of the negative log partial likelihood.
#'
#' @keywords internal

neglog_pl_beta = function(odat, beta, gammas.hat, events_id, ind.wt){
  Y = odat$Y
  Y.mat = cbind(rep(1,length(Y)), log(Y), Y, Y^2)

  gamma1 = gammas.hat[1:4]; gamma2 = gammas.hat[-c(1:4)]

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

    if (ind.wt == 0){
      log_parlike = log_parlike - Xi%*%beta - log(1+rho0.hat*exp(Xi%*%gamma2)) + log(sum(exp(Xj%*%beta)*(1+rho0.hat*exp(Xj%*%gamma2))))
    }
    if (ind.wt == 1){
      Wi = odat$rwt[id]
      log_parlike = log_parlike - Wi*(Xi%*%beta + log(1+rho0.hat*exp(Xi%*%gamma2)) - log(sum(exp(Xj%*%beta)*(1+rho0.hat*exp(Xj%*%gamma2)))))
    }
  }
  return(log_parlike)
}



################ 4. Main function for estimating coefficients and their se ################
#' Estimate disease-specific survival model parameters using a proportional cause ratio model
#'
#' Implements the simulation and estimation procedure for modeling
#' disease-specific survival in observational studies with missing cause-of-death
#' information. The function leverages cause-of-death information from clinical
#' trial data through a proportional cause ratio model.
#'
#' Clinical trial and observational study datasets are generated under
#' piecewise-constant cause-specific proportional hazards models. The
#' proportional cause ratio model parameters are estimated from the clinical
#' trial data, and the disease-specific survival regression coefficients are
#' estimated from the observational study data. Perturbation resampling is used
#' to estimate standard errors for the disease-specific survival regression
#' coefficients.
#'
#' @importFrom stats coef glm nlm rbinom rexp runif sd
#' @param Nc Number of subjects in the clinical trial dataset.
#' @param No Number of subjects in the observational study dataset.
#' @param p Probability parameter for generating the binary covariate.
#' @param r Rate parameter for the exponential censoring distribution.
#' @param cuts_rho Numeric vector of cut-points for the piecewise-constant
#'   baseline hazards.
#' @param lamB Numeric vector of disease-specific baseline hazard values.
#' @param lamO Numeric vector of other-cause baseline hazard values.
#' @param beta Numeric vector of true regression coefficients for
#'   disease-specific death.
#' @param alpha Numeric vector of true regression coefficients for other-cause death.
#' @param par1 Numeric vector of initial values for beta estimation.
#' @param B_perturb Number of perturbation resampling replicates used to estimate
#'   the standard errors of the disease-specific survival regression coefficients.
#'
#' @return A list with components:
#' \describe{
#'   \item{gamma_est}{Estimated proportional cause ratio model coefficients.}
#'   \item{beta_est}{Estimated disease-specific survival regression coefficients.}
#'   \item{gamma_se}{Estimated standard errors of the proportional cause ratio model coefficients.}
#'   \item{beta_se}{Estimated standard errors of the disease-specific survival regression coefficients.}
#'   \item{gamma1_cov}{Estimated variance-covariance matrix for the baseline time-dependent component of the proportional cause ratio model.}
#' }
#'
#' @export
#'
#' @examples
#' # Piecewise-constant baseline hazards
#' lamB <- seq(0.18, 0.08, length.out = 10)
#' lamO <- seq(0.05, 0.10, length.out = 10)
#'
#' # Cut-points defining 10 intervals
#' cuts_rho <- c(0.241, 0.506, 0.813, 1.164, 1.603,
#'               2.130, 2.843, 3.868, 5.701)
#'
#' set.seed(1)
#'
#' res <- dssPCR_pwK(Nc = 200, No = 300, p = 0.5, r = 0.12,
#'                   cuts_rho = cuts_rho, lamB = lamB, lamO = lamO,
#'
#'   # True covariate effects
#'   beta  = c(-0.5, 1),
#'   alpha = c(-0.1, 0.5),
#'
#'   # Initial values for beta estimation
#'   par1 = c(0, 0),
#'
#'   # Number of perturbation replicates
#'   B_perturb = 20
#' )
#'
#' res$gamma_est
#' res$beta_est
#' res$gamma_se
#' res$beta_se
#' res$gamma1_cov

dssPCR_pwK <- function(Nc, No, p, r, cuts_rho, lamB, lamO,
                      beta, alpha, par1, B_perturb = 200){
  #1) Generate clinical trial dataset
  cdat = cdata_genK(Nc, p, r, cutsB = cuts_rho, lamB, cutsO = cuts_rho, lamO, beta, alpha)

  #2) Estimate gamma's and their se
  mod_gamma = gamma_fun(data = cdat, ind.wt = 0)
  gammas.hat = mod_gamma$coeff
  gammas.se = mod_gamma$coeff_se
  gammas.cov = mod_gamma$vcov

  #3) Generate observational study dataset and order it by observed event time
  odat = odata_genK(No, p, r, cutsB = cuts_rho, lamB, cutsO = cuts_rho, lamO, beta, alpha)
  odat = odat[order(odat$Y),] # order by event time
  events_id = which(odat$delta == 1) # observed events index

  #4) Estimate beta's
  # par1 = c(0,0)
  beta.hat = nlm(f = neglog_pl_beta, par1, odat = odat, gammas.hat = gammas.hat,
                 events_id = events_id, ind.wt = 0)$estimate # Use nlm() to minimize negative log partial likelihood

  #5) Compute beta's se using perturbation
  beta.hat.pr = data.frame(matrix(nrow = B_perturb, ncol = 2))
  for (k in 1:B_perturb) {
    cdat.pr = cdat
    cdat.pr$rwt = rexp(Nc,1) # generating random weights using standard exponential distribution
    mod_gamma.pr = gamma_fun(data = cdat.pr, ind.wt = 1)
    gammas.hat.pr = mod_gamma.pr$coeff

    odat.pr = odat
    odat.pr$rwt = rexp(No,1)

    beta.hat.pr[k,] = nlm(f = neglog_pl_beta, par1, odat = odat.pr, gammas.hat = gammas.hat.pr,
                          events_id = events_id, ind.wt = 1)$estimate
  }

  beta.se = apply(beta.hat.pr,2,sd)
  coeff.est = c(gammas.hat, beta.hat, gammas.se, beta.se)
  return(list(
    gamma_est = gammas.hat,
    beta_est = beta.hat,
    gamma_se = gammas.se,
    beta_se = beta.se,
    gamma1_cov = gammas.cov[1:4,1:4]
  ))
}
