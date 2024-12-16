###############################################################################
# Load Packages
###############################################################################
library(maxLik)
library(rlist)
library(matrixcalc)
library(Rfast)
library(sasLM)

###############################################################################
# Helper Functions
###############################################################################

# Rotation matrix R for a given angle delta
R <- function(delta) {
  # Create a 2x2 rotation matrix
  R <- matrix(NA, 2, 2)
  R[1, 1] <- cos(delta)
  R[1, 2] <- -sin(delta)
  R[2, 1] <- sin(delta)
  R[2, 2] <- cos(delta)
  return(R)
}

# Compute the principal matrix square root using spectral decomposition
matroot <- function(mat) {
  ev <- eigen(mat)
  return(ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors))
}

###############################################################################
# Likelihood-Related Functions
###############################################################################

# Log-Likelihood function for t-distributed errors (xi)
loglike <- function(xi, df) {
  df1 <- as.vector(df)[1] 
  df2 <- as.vector(df)[2]
  
  xi1 <- xi[, 1]
  xi2 <- xi[, 2]
  
  c1 <- log(gamma((df1 + 1) / 2)) - log(gamma(df1 / 2)) - 0.5 * log(pi * (df1 - 2))
  c2 <- log(gamma((df2 + 1) / 2)) - log(gamma(df2 / 2)) - 0.5 * log(pi * (df2 - 2))
  
  l1 <- -(df1 + 1) / 2 * log(1 + xi1^2 / (df1 - 2))
  l2 <- -(df2 + 1) / 2 * log(1 + xi2^2 / (df2 - 2))
  
  l <- sum(c1 + l1 + c2 + l2)
  return(l)
}

# Score function u (derivative of log-likelihood w.r.t. delta indirectly)
score <- function(xi, df) {
  df1 <- as.vector(df)[1] 
  df2 <- as.vector(df)[2]
  
  if (any(class(xi) == "matrix")) {
    xi1 <- xi[, 1]
    xi2 <- xi[, 2]
  } else {
    xi1 <- xi[1]
    xi2 <- xi[2]
  }
  
  u <- (df2 + 1) * xi1 * xi2 / (df2 - 2 + xi2^2) - (df1 + 1) * xi1 * xi2 / (df1 - 2 + xi1^2)
  return(u)
}

# Derivative of u with respect to delta
d_u_delta_function <- function(xi, df) {
  df1 <- as.vector(df)[1] 
  df2 <- as.vector(df)[2]
  
  xi1 <- xi[, 1]
  xi2 <- xi[, 2]
  
  c1 <- (df2 + 1) * ((df2 - 2) * (xi2^2 - xi1^2) + xi2^4 + xi1^2 * xi2^2)
  c2 <- (df2 - 2 + xi2^2)^2
  c3 <- (df1 + 1) * ((df1 - 2) * (xi1^2 - xi2^2) + xi1^4 + xi1^2 * xi2^2)
  c4 <- (df1 - 2 + xi1^2)^2
  
  result <- (c1 / c2) + (c3 / c4)
  return(result)
}

# Derivative of delta w.r.t. omega 
d_delta_omega_function <- function(theta, d_u_delta, d_delta_omega) {
  omega <- as.vector(theta)[1]
  beta <- as.vector(theta)[2]
  kappa <- as.vector(theta)[3]

  return(1 - beta + (beta + kappa * d_u_delta) * d_delta_omega)
}

# Derivative of delta w.r.t. beta
d_delta_beta_function <- function(theta, delta, d_u_delta, d_delta_beta) {
  omega <- as.vector(theta)[1]
  beta <- as.vector(theta)[2]
  kappa <- as.vector(theta)[3]
  
  return((-1) * omega + delta + (beta + kappa * d_u_delta) * d_delta_beta)
}

# Derivative of delta w.r.t. kappa
d_delta_kappa_function <- function(theta, u, d_u_delta, d_delta_kappa) {
  omega <- as.vector(theta)[1]
  beta <- as.vector(theta)[2]
  kappa <- as.vector(theta)[3]
  
  return(u + (beta + kappa * d_u_delta) * d_delta_kappa)
}

# Delta function: updates delta based on the score-driven process
delta_function <- function(theta, delta, u) {
  omega <- as.vector(theta)[1]
  beta <- as.vector(theta)[2]
  kappa <- as.vector(theta)[3]
  
  return(omega * (1 - beta) + beta * delta + kappa * u)
}

###############################################################################
# Score-Driven Dynamic Rotation
###############################################################################

# Likelihood estimation function
lambda_est_max <- function(theta, et, df) {
  omega <- as.vector(theta)[1]
  beta <- as.vector(theta)[2]
  kappa <- as.vector(theta)[3]
  
  df1 <- as.vector(df)[1] 
  df2 <- as.vector(df)[2]
  
  xi_t <- matrix(0, NROW(et), 2)
  
  # Initialization
  delta <- omega
  xi <- t(t(R(delta)) %*% (et[1, ]))
  xi_t[1, ] <- xi
  u <- score(xi, df)
  
  i <- 2
  while (i <= NROW(et)) {
    delta <- delta_function(theta = c(omega, beta, kappa), delta = delta, u = u)
    xi <- t(t(R(delta)) %*% et[i, ])
    xi_t[i, ] <- xi
    u <- score(xi, df)
    i <- i + 1
  }
  
  loglik <- loglike(xi_t, df)
  return(loglik)
}

# Analytics function returns various intermediate quantities
analytics <- function(theta, et, df) {
  omega <- as.vector(theta)[1]
  beta <- as.vector(theta)[2]
  kappa <- as.vector(theta)[3]
  
  df1 <- as.vector(df)[1] 
  df2 <- as.vector(df)[2]
  
  delta_t <- matrix(0, NROW(et), 1)
  xi_t  <- matrix(0, NROW(et), 2)
  
  u_t  <- matrix(0, NROW(et), 1)        # Score (u)
  d_u_delta_t <- matrix(0, NROW(et), 1) # Derivative of u wrt. delta
  d_l_theta_t <- matrix(0, NROW(et), 3) # Derivative of log-likelihood wrt. theta
  d_l_delta_t <- matrix(0, NROW(et), 1) # Derivative of log-likelihood wrt. delta (possibly u)
  d_delta_theta_t <- matrix(0, NROW(et), 3) # Derivative of delta wrt. theta
  
  # Initialization
  delta_init <- omega
  d_delta_omega <- 0
  d_delta_beta <- 0
  d_delta_kappa <- 0
  d_u_delta <- 0
  u <- 0
  
  # Initial state
  delta <- delta_function(theta, delta_init, u)
  delta_t[1, ] <- delta
  
  xi <- t(t(R(delta)) %*% et[1, ])
  xi_t[1, ] <- xi
  
  u_t[1, ] <- score(xi, df)
  
  d_delta_omega <- d_delta_omega_function(theta, d_u_delta, d_delta_omega)
  d_delta_beta <- d_delta_beta_function(theta, delta = delta_init, d_u_delta, d_delta_beta)
  d_delta_kappa <- d_delta_kappa_function(theta, u, d_u_delta, d_delta_kappa)
  
  d_u_delta <- d_u_delta_function(xi, df)
  d_u_delta_t[1, ] <- d_u_delta
  
  d_delta_theta <- c(d_delta_omega, d_delta_beta, d_delta_kappa)
  d_delta_theta_t[1, ] <- d_delta_theta
  
  u <- u_t[1, ] 
  d_l_theta_t[1, ] <- u * d_delta_theta
  
  # Iteration
  i <- 2
  while (i <= NROW(et)) {
    delta <- delta_function(theta = c(omega, beta, kappa), delta = delta, u = u)
    delta_t[i, ] <- delta
    
    xi <- t(t(R(delta)) %*% et[i, ])
    xi_t[i, ] <- xi
    
    d_delta_omega <- d_delta_omega_function(theta, d_u_delta, d_delta_omega)
    d_delta_beta <- d_delta_beta_function(theta, delta = delta_t[i - 1, ], d_u_delta, d_delta_beta)
    d_delta_kappa <- d_delta_kappa_function(theta, u, d_u_delta, d_delta_kappa)
    
    d_u_delta <- d_u_delta_function(xi, df)
    d_u_delta_t[i, ] <- d_u_delta
    
    d_delta_theta <- c(d_delta_omega, d_delta_beta, d_delta_kappa)
    d_delta_theta_t[i, ] <- d_delta_theta
    
    u_t[i, ] <- score(xi, df)
    u <- u_t[i, ]
    d_l_theta_t[i, ] <- crossprod(u_t[i, ], d_delta_theta)
    
    i <- i + 1
  }
  
  return_list <- list(u_t, delta_t, d_u_delta_t, xi_t, d_delta_theta_t, d_l_theta_t)
  return(return_list)
}

# Gradient function returns the gradient of the log-likelihood
grad_function <- function(theta, et, df) {
  omega <- as.vector(theta)[1]
  beta <- as.vector(theta)[2]
  kappa <- as.vector(theta)[3]
  
  df1 <- as.vector(df)[1] 
  df2 <- as.vector(df)[2]
  
  delta_t <- matrix(0, NROW(et), 1)
  xi_t  <- matrix(0, NROW(et), 2)
  
  u_t  <- matrix(0, NROW(et), 1)      
  d_u_delta_t <- matrix(0, NROW(et), 1)
  d_l_theta_t <- matrix(0, NROW(et), 3)
  d_l_delta_t <- matrix(0, NROW(et), 1)
  d_delta_theta_t <- matrix(0, NROW(et), 3)
  
  # Initialization
  delta_init <- omega
  d_delta_omega <- 0
  d_delta_beta <- 0
  d_delta_kappa <- 0
  d_u_delta <- 0
  u <- 0
  
  delta <- delta_function(c(omega, beta, kappa), delta_init, u)
  delta_t[1, ] <- delta
  
  xi <- t(t(R(delta)) %*% et[1, ])
  xi_t[1, ] <- xi
  
  u_t[1, ] <- score(xi, df)
  
  d_delta_omega <- d_delta_omega_function(theta, d_u_delta, d_delta_omega)
  d_delta_beta <- d_delta_beta_function(theta, delta = delta_init, d_u_delta, d_delta_beta)
  d_delta_kappa <- d_delta_kappa_function(theta, u, d_u_delta, d_delta_kappa)
  
  d_u_delta <- d_u_delta_function(xi, df)
  d_u_delta_t[1, ] <- d_u_delta
  
  d_delta_theta <- c(d_delta_omega, d_delta_beta, d_delta_kappa)
  d_delta_theta_t[1, ] <- d_delta_theta
  
  u <- u_t[1, ]
  d_l_theta_t[1, ] <- u * d_delta_theta
  
  i <- 2
  while (i <= NROW(et)) {
    delta <- delta_function(theta = c(omega, beta, kappa), delta = delta, u = u)
    delta_t[i, ] <- delta
    
    xi <- t(t(R(delta)) %*% et[i, ])
    xi_t[i, ] <- xi
    
    d_delta_omega <- d_delta_omega_function(theta, d_u_delta, d_delta_omega)
    d_delta_beta <- d_delta_beta_function(theta, delta = delta_t[i - 1, ], d_u_delta, d_delta_beta)
    d_delta_kappa <- d_delta_kappa_function(theta, u, d_u_delta, d_delta_kappa)
    
    d_u_delta <- d_u_delta_function(xi, df)
    d_u_delta_t[i, ] <- d_u_delta
    
    d_delta_theta <- c(d_delta_omega, d_delta_beta, d_delta_kappa)
    d_delta_theta_t[i, ] <- d_delta_theta
    
    u_t[i, ] <- score(xi, df)
    u <- u_t[i, ]
    d_l_theta_t[i, ] <- u_t[i, ] * d_delta_theta
    
    i <- i + 1
  }
  
  return_list <- d_l_theta_t
  return(return_list)
}
