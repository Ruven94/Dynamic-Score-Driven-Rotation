###############################################################################
# Setup and Data Loading
###############################################################################
# Required Packages
library(quantmod)
library(PerformanceAnalytics)
library(BEKKs)

# Download data: Gold ETF (GLD) and MSCI World ETF (URTH) from Yahoo Finance
getSymbols("GLD", src = "yahoo", from = "2014-01-01", to = Sys.Date())
getSymbols("URTH", src = "yahoo", from = "2014-01-01", to = Sys.Date())

###############################################################################
# Data Preparation
###############################################################################
# Extract closing prices
gold_prices <- Cl(GLD)
msci_prices <- Cl(URTH)

# Compute daily log-returns
gold_returns <- dailyReturn(gold_prices, type = "log")
msci_returns <- dailyReturn(msci_prices, type = "log")

# De-mean the returns
gold_returns_demean <- gold_returns - mean(gold_returns)
msci_returns_demean <- msci_returns - mean(msci_returns)

# Combine into a matrix
data <- cbind(gold_returns_demean, msci_returns_demean)
data <- as.matrix(data)

###############################################################################
# BEKK Model Estimation
###############################################################################
# Specify an asymmetric BEKK model
bekk_s1 <- bekk_spec(model = list(type = "bekk", asymmetric = TRUE))
# Fit the BEKK model to the data
bekk_m1 <- bekk_fit(bekk_s1, data)
# Extract residuals (innovations)
et <- as.matrix(bekk_m1$e_t)

###############################################################################
# Dynamic Score-Driven Rotation Estimation
###############################################################################
# The function dynamic_score_driven_rot_estimation is assumed to be defined 
# and available in the environment. It takes 'et' and ranges for df1, df2.
# Here, we run it for df1 and df2 in {5,6}.
results <- dynamic_score_driven_rot_estimation(et,
                                               df1_range = seq(5, 6),
                                               df2_range = seq(5, 6))

###############################################################################
# Analyze Results
###############################################################################
# Extract parameters and compute xi and delta series using the analytics function
theta_est <- results[1, 3:5]
df_est    <- results[1, 1:2]

xi    <- analytics(theta = theta_est, et = et, df = df_est)[[4]]
delta <- analytics(theta = theta_est, et = et, df = df_est)[[2]]

# Plot the estimated delta series
plot(delta, type = "l", main = "Estimated Delta over Time", ylab = "Delta", xlab = "Time")