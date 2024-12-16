# Dynamic Score-Driven Rotation Estimation

This repository implements a dynamic score-driven rotation model using BEKK GARCH residuals and a t-distributed likelihood. The code provides:

1. **Helper Functions**: For matrix operations, likelihood calculations, scoring, and gradients.
2. **Dynamic Estimation**: A routine to estimate parameters using a grid search followed by optimization via `maxBHHH`.
3. **Empirical Example**: An application using daily return data from the Gold ETF (GLD) and the MSCI World ETF (URTH).

## File Structure

- **Helper.R**  
  Contains helper functions and utilities, including:
  - Rotation matrix function
  - Matrix square root function
  - Log-likelihood and score functions
  - Derivatives and analytical gradients for dynamic rotation estimation
  
- **Estimation.R**  
  Implements the `dynamic_score_driven_rot_estimation` function. Given a range of `df1` and `df2` values, this script:
  - Performs a grid search over `(omega, beta, kappa)` to find good initial parameters.
  - Uses `maxBHHH` (from `maxLik`) to optimize the parameters from the best initial candidates.
  - Returns sorted results based on the highest likelihood.

- **Empirical.R**  
  Shows a practical example:
  - Fetches historical data for GLD and URTH from Yahoo Finance.
  - Computes log-returns and fits an asymmetric BEKK model to obtain residuals.
  - Runs the dynamic score-driven rotation estimation on the residuals.
  - Extracts and plots the estimated delta series over time.

## Requirements

- **R Version**: Tested on R 4.0+
- **R Packages**:  
  - `quantmod` (for data retrieval and financial time series handling)  
  - `PerformanceAnalytics` (for return calculations)  
  - `BEKKs` (for fitting BEKK GARCH models)  
  - `maxLik` (for BHHH optimization)  
  - `rlist`, `matrixcalc`, `Rfast`, `sasLM` (for various mathematical and list utilities)
  
To install missing packages:
```r
install.packages(c("quantmod", "PerformanceAnalytics", "BEKKs", "maxLik", "rlist", "matrixcalc", "Rfast", "sasLM"))
