# als-helper-functions.R
#
# External functions for the ALS simulation

expit <- function(x) {
  return(exp(x) / (1+exp(x)))
}

gen_covariates <- function(n, meanvec, sigma) {
  
  # Simulate continuous convariates (assumed to be normal)
  X <- mvtnorm::rmvnorm(n,mean=meanvec,sigma = sigma)
  
  # Also simulate sex (assume to be independent)
  sex <- rbinom(n,size=1,prob=0.5)
  
  # Last column is (assumed to be) baseline ALSFRS-R. Round numbers
  X[,ncol(X)] <- round(X[,ncol(X)])
  
  X <- data.frame(sex,X)
  
  return(X)
}

tmt_model <- function(X,alphas) {
  #' Code for generating "treatment" assignment
  #'
  #'
  #'
  
  n <- nrow(X) ; X <- scale(X) ; X <- as.matrix(cbind(1,X))
  probs <- expit(X%*%alphas) 
  A <- rbinom(n,size=1,prob = probs)
  
  return(data.frame(A))
  
}

pot_outcome_model <- function(X, tau, betas, gammas) {
  #' Code for genetating potential outcomes
  #'
  #'
  #'
  #'
  n <- nrow(X) ; X <- as.matrix(X)

  # Means for potential outcomes  
  mu0 <- X%*%betas
  mu1 <- tau + X%*%gammas
  
  # Generate potential outcomes
  Y0 <- extraDistr::rtpois(n,mu0,a=0,b=48)
  Y1 <- Y0 + extraDistr::rtpois(n,mu1,a=0,b=48-Y0)
  
  return(data.frame(Y0,Y1))
  
}

gen_data <- function(n, # sample size
                     meanvec, sigma, # covariate params
                     alphas, # ps coefficients
                     tau, betas, gammas) { # outcome model params
  
  # Generate covariates
  X <- gen_covariates(n,meanvec, sigma)
  
  # Treatment model
  A <- tmt_model(X,alphas)
  
  # Potential outcomes
  pot_outs <- pot_outcome_model(X,tau,betas,gammas)
  
  # Form the df
  df <- data.frame(cbind(X,A,pot_outs))
  df <- df %>% mutate(Y = ifelse(A==1,Y1,Y0))
  
  return(df)
  
}

run_sim <- function(params) {
  
  # Extract elements from params
  list2env(pararms, envir = environment())
  
  # Get all vars that start with X
  vars <- grep("^(?!Y|A)", names(df), value = TRUE, perl = TRUE)
  formula_string <- paste("A ~", paste(vars, collapse=" + "))
  formula_obj <- as.formula(formula_string)
  
  for (s in 1:nsim) {
    
    # simulate data
    df <- gen_data(n, # sample size
                   meanvec, sigma, # covariate params
                   alphas, # ps coefficients
                   tau, betas, gammas)
    
    # extraxt coeffs for ps regression

    # ps matching
    psm <-MatchIt::matchit(formula_obj, data = df,
                  method = "nearest", estimand = "ATT")
    lm_psm <- lm(Y ~ A, data=psm,weights = 'weights')
    avg_comparisons(fit1, variables = "A",
                    vcov = ~subclass,
                    newdata = subset(md, A == 1),
                    wts = "weights")
    
    
    
    # exact matching
    exm <- MatchIt::matchit(A ~ ., data=df_psm, method='cardinality')
    
    
    
  }
}
