# als-simulations.R 
#
# Simulation code for ALS simulations
#
# The goal is to produce simulations for 3 scenarios:
# - Y(1)_i=Y(0)_i for all i
# - Y(1)_i=Y(0)_i + t_i for all i (t_i normal, possibly)
# - Y(1)_i=Y(0)_i + f(X_i) for all i
#
# We'll compare the following balancing methods for causal effect estimation:
# - pscore matching vs 
# - pscore weighting vs
# -
#
# Outcome: ALSFRS-R (0-48 score)
#
#
# Covariates to "include" (i.e. make realistic simulations):
# sex, age, time of day, BMI, ... 
#-------------------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(SuperLearner)
library(kableExtra)
require(haven)
require(lme4)
require(tableone)
require(broom)
require(coin)
require(extraDistr)
#-------------------------------------------------------------------------------
# External functions
source('als-helper-functions.R')
#-------------------------------------------------------------------------------
# Fixed simulation parameters

# Covariates: age, bmi, months since diagnosis, baseline ALSFRSR
meanvec <- c(55,25,6, 32)
sigma <- diag(rep(1,4)) ; diag(sigma) <- c(4,2,1,16)
n <- 100

# sigma <- matrix(c(4, 0.1, -0.3,
#                   0.1, 2, 0,
#                   -0.3, 0, 1), nrow=3, byrow=F)

# Baseline pot. outcome DGP (last col is 1 b.c. baseline alsfrs-r)
betas <- c(1,-0.1, -0.1, -0.8, 1)

# "Treatment" assignment
alphas <- c(0 ,-0.2, -0.1, 0.4, 0.6, 0.5)
#-------------------------------------------------------------------------------
# Simulation scenario 1: sharp causal null holds

gammas <- rep(0, 5) ; tau <- 0

# covariates (way I'm thinking of doing this, eventually, is to approx dist of 
# covariates in the data we have with a MVN and use the estimated params to gen
# data)

df <- gen_data(n, # sample size
         meanvec, sigma, # covariate params
         alphas, # ps coefficients
         tau, betas, gammas)

X <- gen_covariates(n,meanvec, sigma)
A <- tmt_model(X,alphas)
pot_outs <- pot_outcome_model(X,tau,betas,gammas)

df <- data.frame(cbind(X,A,pot_outs))
df <- df %>% mutate(Y = ifelse(A==1,Y1,Y0))
#-------------------------------------------------------------------------------
# Simulation scenario 2: constant tmt effect

gammas <- rep(0, 5) ; tau <- 0

# covariates (way I'm thinking of doing this, eventually, is to approx dist of 
# covariates in the data we have with a MVN and use the estimated params to gen
# data)

X <- gen_covariates(n,meanvec, sigma)
A <- tmt_model(X,alphas)
pot_outs <- pot_outcome_model(X,tau,betas,gammas)

df <- data.frame(cbind(X,A,pot_outs))
df <- df %>% mutate(Y = ifelse(A==1,Y1,Y0)) # observed outcome
#-------------------------------------------------------------------------------
# Simulation scenario 3: heterogenous treatment effects

gammas <- rep(0, 5) ; tau <- 0

# covariates (way I'm thinking of doing this, eventually, is to approx dist of 
# covariates in the data we have with a MVN and use the estimated params to gen
# data)

X <- gen_covariates(n,meanvec, sigma)
A <- tmt_model(X,alphas)
pot_outs <- pot_outcome_model(X,tau,betas,gammas)

df <- data.frame(cbind(X,A,pot_outs))
df <- df %>% mutate(Y = ifelse(A==1,Y1,Y0))




