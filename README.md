# dssPCR

## Overview

`dssPCR` implements methods for modeling disease-specific survival (DSS)
in observational studies with missing cause-of-death information using a
proportional cause ratio model.

The package leverages cause-of-death information from clinical trial data to
estimate DSS regression parameters in observational studies where cause of
death is unavailable.

## Installation

```r
# install.packages("devtools")
devtools::install_github("YunyiWang915/dssPCR")
```



## Key Arguments

| Argument | Description |
|---|---|
| `Nc` | Number of subjects in the clinical trial dataset |
| `No` | Number of subjects in the observational study dataset |
| `p` | Probability parameter for generating the binary covariate |
| `r` | Rate parameter for the exponential censoring distribution |
| `cuts_rho` | Cut-points defining the piecewise-constant hazard intervals |
| `lamB` | Disease-specific baseline hazard values across intervals |
| `lamO` | Other-cause baseline hazard values across intervals |
| `beta` | True regression coefficients for disease-specific death |
| `alpha` | True regression coefficients for other-cause death |
| `par1` | Initial values used in beta optimization |
| `B_perturb` | Number of perturbation resampling replicates for standard error estimation |



## Quick Start: Full Estimation with `dssPCR_pwK()`

```r
library(dssPCR)

set.seed(1)
Nc <- 1000
No <- 1000
r <- 0.12

lamB <- seq(0.18, 0.08, length.out = 10)
lamO <- seq(0.05, 0.10, length.out = 10)

cuts_rho <- c(0.241, 0.506, 0.813, 1.164, 1.603, 2.130, 2.843, 3.868, 5.701)

res <- dssPCR_pwK(Nc = Nc, No = No, p = 0.5, r = r, cuts_rho = cuts_rho,
                  lamB = lamB, lamO = lamO,
                  beta = c(-0.5, 1), alpha = c(-0.1, 0.5),
                  par1 = c(0, 0), B_perturb = 200)

res$gamma_est
res$beta_est
res$gamma_se
res$beta_se
```

## Step-by-Step Workflow

The full procedure can also be run step by step to illustrate how the package works internally.

### Step 1: Specify simulation settings

```r
Nc <- 1000
No <- 1000
p <- 0.5
r <- 0.12

lamB <- seq(0.18, 0.08, length.out = 10)
lamO <- seq(0.05, 0.10, length.out = 10)

cuts_rho <- c(0.241, 0.506, 0.813, 1.164, 1.603, 2.130, 2.843, 3.868, 5.701)

beta  <- c(-0.5, 1)
alpha <- c(-0.1, 0.5)
```

### Step 2: Generate a clinical trial dataset

```r
set.seed(2)

cdat <- cdata_genK(Nc = Nc, p = p, r = r, 
                   cutsB = cuts_rho, lamB = lamB, 
                   cutsO = cuts_rho, lamO = lamO,
                   beta = beta, alpha = alpha)

head(cdat)
```

### Step 3: Estimate proportional cause ratio model parameters

```r
mod_gamma <- gamma_fun(data = cdat, ind.wt = 0)

gammas.hat <- mod_gamma$coeff
gammas.se  <- mod_gamma$coeff_se

gammas.hat
gammas.se
```

### Step 4: Generate an observational study dataset

```r
odat <- odata_genK(No = No, p = p, r = r, 
                   cutsB = cuts_rho, lamB = lamB,
                   cutsO = cuts_rho, lamO = lamO,
                   beta = beta, alpha = alpha)

head(odat)
```

### Step 5: Run the full estimation procedure

```r
res <- dssPCR_pwK(Nc = Nc, No = No, p = p, r = r, 
                  cuts_rho = cuts_rho, lamB = lamB, 
                  lamO = lamO, beta = beta, alpha = alpha,
                  par1 = c(0, 0), B_perturb = 200)

res$gamma_est
res$beta_est
res$gamma_se
res$beta_se
```

## Main Functions

| Function | Description |
|---|---|
| `dssPCR_pwK()` | Runs the full simulation and estimation procedure |
| `cdata_genK()` | Generates a clinical trial dataset with observed cause-of-death information |
| `odata_genK()` | Generates an observational study dataset with missing cause-of-death information |
| `gamma_fun()` | Estimates proportional cause ratio model parameters |

---

## Reference

Wang Y, Shen Y, Ning J. Modeling Disease-specific Survival in Observational Studies 
with Missing Cause of Death: Leveraging Information from Clinical Trial Data.
