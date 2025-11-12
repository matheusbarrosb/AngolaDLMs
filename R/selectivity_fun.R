library(tidyverse)

# 1. SELECTIVITY FUNCTIONS ----------------------------
# Double-normal (Gaussian) selectivity - dome-shaped
dome_selectivity = function(L, SL50, SL_sd, SR50, SR_sd) {
  sel_asc     = pnorm(L, mean = SL50, sd = SL_sd)
  sel_desc    = 1 - pnorm(L, mean = SR50, sd = SR_sd)
  selectivity = sel_asc * sel_desc
  return(selectivity)
}

# Double-logistic selectivity - dome-shaped
double_logistic = function(L, SL50, SL95, SR50, SR95) {
  sel_asc     = 1 / (1 + exp(-log(19) * (L - SL50) / (SL95 - SL50)))
  sel_desc    = 1 / (1 + exp(-log(19) * (L - SR50) / (SR95 - SR50)))
  selectivity = sel_asc * (1 - sel_desc)
  return(selectivity)
}

# 2. NON-PARAMETRIC SELECTIVITY ESTIMATION ------------
estimate_selectivity_nonparametric = function(data,
                                               species_name,
                                               bin_width = 1,
                                               smooth_span = 0.3,
                                               n_attempts = 10) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("Estimating selectivity for:", species_name, "\n")
  cat(strrep("=", 70), "\n\n")
  
  years = sort(unique(data$Ano))
  results_list = list()
  
  for (yr in years) {
    
    cat("Processing year:", yr, "\n")
    
    # Filter data for this year
    data_year = data %>% filter(Ano == yr)
    
    # Expand data
    expanded_data = data_year %>% uncount(N)
    
    n_samples <- nrow(expanded_data)
    cat("  Sample size:", n_samples, "\n")
    
    if (n_samples < 50) {
      cat("  WARNING: Small sample size, skipping year\n\n")
      next
    }
    
    # Create length bins
    min_length = floor(min(expanded_data$LT))
    max_length = ceiling(max(expanded_data$LT))
    length_bins = seq(min_length, max_length + bin_width, by = bin_width)
    
    # Observed frequency
    obs_hist = hist(expanded_data$LT, breaks = length_bins, plot = FALSE)
    observed_freq = obs_hist$counts
    bin_mids = obs_hist$mids
    
    mode_length = bin_mids[which.max(observed_freq)]
    q25 = quantile(expanded_data$LT, 0.25)
    q75 = quantile(expanded_data$LT, 0.75)
    
    left_side_idx = bin_mids <= mode_length
    left_freq = observed_freq[left_side_idx]
    left_mids = bin_mids[left_side_idx]
    
    # Fit smooth curve to left side and extrapolate
    if (sum(left_side_idx) > 5) {
      # Use loess smoothing
      smooth_left = loess(left_freq ~ left_mids, span = smooth_span)
      
      # Predict for all lengths using the ascending pattern
      # Mirror the left side to create expected right side
      expected_freq = numeric(length(bin_mids))
      
      for (i in seq_along(bin_mids)) {
        # Distance from mode
        dist_from_mode = abs(bin_mids[i] - mode_length)
        mirror_length = mode_length - dist_from_mode
        
        if (mirror_length >= min(left_mids) & mirror_length <= max(left_mids)) {
          expected_freq[i] = predict(smooth_left, newdata = data.frame(left_mids = mirror_length))
        } else {
          expected_freq[i] = observed_freq[i]  # Use observed for very small lengths
        }
      }
      
      expected_freq[expected_freq < 0] = 0.1
      expected_freq = pmax(expected_freq, 0.1)  # Avoid zeros
      
    } else {
      # Fallback: use log-normal as expected distribution
      log_mean = log(mode_length)
      log_sd = sd(log(expanded_data$LT))
      expected_freq = dlnorm(bin_mids, log_mean, log_sd) * sum(observed_freq)
      expected_freq = pmax(expected_freq, 0.1)
    }
    
    # --- NEGATIVE LOG-LIKELIHOOD FUNCTIONS ---
    
    # For dome selectivity
    nll_dome = function(params) {
      SL50  = params[1]
      SL_sd = params[2]
      SR50  = params[3]
      SR_sd = params[4]
      
      # Calculate selectivity
      sel = dome_selectivity(bin_mids, SL50, SL_sd, SR50, SR_sd)
      
      # Expected observed = expected population * selectivity
      predicted = expected_freq * sel
      predicted = pmax(predicted, 1e-10)
      
      # Negative log-likelihood (Poisson)
      nll = -sum(dpois(observed_freq, lambda = predicted, log = TRUE))
      
      if (!is.finite(nll)) return(1e10)
      return(nll)
    }
    
    # For double-logistic selectivity
    nll_logistic = function(params) {
      SL50 = params[1]
      SL95 = params[2]
      SR50 = params[3]
      SR95 = params[4]
      
      sel = double_logistic(bin_mids, SL50, SL95, SR50, SR95)
      
      predicted = expected_freq * sel
      predicted = pmax(predicted, 1e-10)
      
      nll = -sum(dpois(observed_freq, lambda = predicted, log = TRUE))
      
      if (!is.finite(nll)) return(1e10)
      return(nll)
    }
    
    # --- FIT DOME SELECTIVITY ---
    cat("  Fitting dome (double-normal) selectivity...\n")
    
    best_dome = NULL
    best_dome_nll = Inf
    
    for (attempt in 1:n_attempts) {
      
      start_dome <- c(
        SL50  = q25 + rnorm(1, 0, 3),
        SL_sd = abs(rnorm(1, 4, 1)),
        SR50  = q75 + rnorm(1, 0, 3),
        SR_sd = abs(rnorm(1, 4, 1))
      )
      
      lower_dome = c(min_length, 0.5, mode_length, 0.5)
      upper_dome = c(mode_length, 15, max_length, 15)
      
      tryCatch({
        fit_dome = optim(
          par     = start_dome,
          fn      = nll_dome,
          method  = "L-BFGS-B",
          lower   = lower_dome,
          upper   = upper_dome,
          control = list(maxit = 1000)
        )
        
        if (fit_dome$value < best_dome_nll) {
          best_dome = fit_dome
          best_dome_nll = fit_dome$value
        }
        
      }, error = function(e) {
        cat("    Dome attempt", attempt, "failed\n")
      })
    }
    
    # --- FIT DOUBLE-LOGISTIC SELECTIVITY ---
    cat("  Fitting double-logistic selectivity...\n")
    
    best_logistic = NULL
    best_logistic_nll = Inf
    
    for (attempt in 1:n_attempts) {
      
      start_SL50 = q25 + rnorm(1, 0, 3)
      start_SR50 = q75 + rnorm(1, 0, 3)
      
      start_logistic = c(
        SL50 = start_SL50,
        SL95 = start_SL50 + abs(rnorm(1, 4, 1)),
        SR50 = start_SR50,
        SR95 = start_SR50 + abs(rnorm(1, 4, 1))
      )
      
      lower_logistic = c(min_length, min_length + 1, 
                          mode_length, mode_length + 1)
      upper_logistic = c(mode_length, mode_length + 10,
                          max_length, max_length + 10)
      
      tryCatch({
        fit_logistic = optim(
          par     = start_logistic,
          fn      = nll_logistic,
          method  = "L-BFGS-B",
          lower   = lower_logistic,
          upper   = upper_logistic,
          control = list(maxit = 1000)
        )
        
        if (fit_logistic$value < best_logistic_nll) {
          best_logistic <- fit_logistic
          best_logistic_nll <- fit_logistic$value
        }
        
      }, error = function(e) {
        cat("    Logistic attempt", attempt, "failed\n")
      })
    }
    
    # --- STORE RESULTS ---
    if (!is.null(best_dome) && !is.null(best_logistic)) {
      
      # Calculate AIC
      k_dome = 4  
      k_logistic = 4
      
      aic_dome = 2 * k_dome + 2 * best_dome$value
      aic_logistic = 2 * k_logistic + 2 * best_logistic$value
      
      delta_aic = aic_logistic - aic_dome
      
      cat("  Dome NLL:", round(best_dome_nll, 2), 
          " | AIC:", round(aic_dome, 2), "\n")
      cat("  Logistic NLL:", round(best_logistic_nll, 2), 
          " | AIC:", round(aic_logistic, 2), "\n")
      cat("  Delta AIC:", round(delta_aic, 2), 
          "(negative = dome better)\n\n")
      
      results_list[[as.character(yr)]] <- list(
        year = yr,
        n_samples = n_samples,
        observed_freq = observed_freq,
        expected_freq = expected_freq,
        length_bins = length_bins,
        bin_mids = bin_mids,
        dome = list(
          params = best_dome$par,
          nll = best_dome$value,
          aic = aic_dome,
          convergence = best_dome$convergence
        ),
        logistic = list(
          params = best_logistic$par,
          nll = best_logistic$value,
          aic = aic_logistic,
          convergence = best_logistic$convergence
        ),
        delta_aic = delta_aic,
        best_model = ifelse(delta_aic < 0, "dome", "logistic")
      )
    }
  }
  
  return(results_list)
}

# 3. ALTERNATIVE: FIT DIRECTLY TO DECLINE RATE

# This method fits selectivity by assuming the right tail decline is due to selectivity
# Not mortality (since we don't know mortality)

estimate_selectivity_from_decline = function(data,
                                              species_name,
                                              bin_width = 1,
                                              n_attempts = 10) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("Estimating selectivity from decline rate:", species_name, "\n")
  cat(strrep("=", 70), "\n\n")
  
  years = sort(unique(data$Ano))
  results_list = list()
  
  for (yr in years) {
    
    cat("Processing year:", yr, "\n")
    
    data_year = data %>% filter(Ano == yr)
    expanded_data = data_year %>% uncount(N)
    
    n_samples = nrow(expanded_data)
    cat("  Sample size:", n_samples, "\n")
    
    if (n_samples < 50) {
      cat("  WARNING: Small sample size, skipping\n\n")
      next
    }
    
    # Create frequency distribution
    min_length = floor(min(expanded_data$LT))
    max_length = ceiling(max(expanded_data$LT))
    length_bins = seq(min_length, max_length + bin_width, by = bin_width)
    
    obs_hist = hist(expanded_data$LT, breaks = length_bins, plot = FALSE)
    observed_freq = obs_hist$counts
    bin_mids = obs_hist$mids
    
    # Find mode (peak)
    mode_idx = which.max(observed_freq)
    mode_length = bin_mids[mode_idx]
    
    # Identify ascending and descending portions
    ascending_idx  = 1:mode_idx
    descending_idx = mode_idx:length(bin_mids)
    
    # Key idea: fit selectivity to the RATIO of observed to expected
    # Expected = constant (flat) or slowly declining
    
    # Simple approach: assume flat expected frequency at mode level
    expected_flat = max(observed_freq)
    
    # Calculate implied selectivity
    implied_sel = observed_freq / expected_flat
    implied_sel = pmin(implied_sel, 1)  # Cap at 1
    implied_sel = pmax(implied_sel, 0.01)  # Floor at 0.01
    
    # Now fit parametric selectivity curves to this implied selectivity
    
    # --- DOME SELECTIVITY ---
    nll_dome_direct = function(params) {
      SL50  = params[1]
      SL_sd = params[2]
      SR50  = params[3]
      SR_sd = params[4]
      
      sel_fitted = dome_selectivity(bin_mids, SL50, SL_sd, SR50, SR_sd)
      
      # Sum of squared errors between implied and fitted selectivity
      sse = sum((implied_sel - sel_fitted)^2)
      
      return(sse)
    }
    
    # --- DOUBLE-LOGISTIC ---
    nll_logistic_direct <- function(params) {
      SL50 = params[1]
      SL95 = params[2]
      SR50 = params[3]
      SR95 = params[4]
      
      sel_fitted = double_logistic(bin_mids, SL50, SL95, SR50, SR95)
      
      sse <- sum((implied_sel - sel_fitted)^2)
      
      return(sse)
    }
    
    # Fit both models
    q25 = quantile(expanded_data$LT, 0.25)
    q75 = quantile(expanded_data$LT, 0.75)
    
    best_dome = NULL
    best_dome_sse = Inf
    
    for (attempt in 1:n_attempts) {
      start_dome <- c(
        SL50 = q25 + rnorm(1, 0, 2),
        SL_sd = abs(rnorm(1, 3, 1)),
        SR50 = q75 + rnorm(1, 0, 2),
        SR_sd = abs(rnorm(1, 3, 1))
      )
      
      tryCatch({
        fit = optim(start_dome, nll_dome_direct, method = "L-BFGS-B",
                     lower = c(min_length, 0.5, mode_length, 0.5),
                     upper = c(mode_length, 10, max_length, 10))
        
        if (fit$value < best_dome_sse) {
          best_dome = fit
          best_dome_sse = fit$value
        }
      }, error = function(e) {})
    }
    
    best_logistic = NULL
    best_logistic_sse = Inf
    
    for (attempt in 1:n_attempts) {
      start_SL50 = q25 + rnorm(1, 0, 2)
      start_SR50 = q75 + rnorm(1, 0, 2)
      
      start_logistic = c(
        SL50 = start_SL50,
        SL95 = start_SL50 + abs(rnorm(1, 3, 0.5)),
        SR50 = start_SR50,
        SR95 = start_SR50 + abs(rnorm(1, 3, 0.5))
      )
      
      tryCatch({
        fit <- optim(start_logistic, nll_logistic_direct, method = "L-BFGS-B",
                     lower = c(min_length, min_length + 0.5, mode_length, mode_length + 0.5),
                     upper = c(mode_length, mode_length + 10, max_length, max_length + 10))
        
        if (fit$value < best_logistic_sse) {
          best_logistic <- fit
          best_logistic_sse <- fit$value
        }
      }, error = function(e) {})
    }
    
    if (!is.null(best_dome) && !is.null(best_logistic)) {
      
      # Use AIC based on SSE (treating as normal errors)
      n = length(observed_freq)
      aic_dome = n * log(best_dome_sse / n) + 2 * 4
      aic_logistic = n * log(best_logistic_sse / n) + 2 * 4
      delta_aic = aic_logistic - aic_dome
      
      cat("  Dome SSE:", round(best_dome_sse, 4), 
          " | AIC:", round(aic_dome, 2), "\n")
      cat("  Logistic SSE:", round(best_logistic_sse, 4), 
          " | AIC:", round(aic_logistic, 2), "\n")
      cat("  Delta AIC:", round(delta_aic, 2), "\n\n")
      
      results_list[[as.character(yr)]] <- list(
        year = yr,
        n_samples = n_samples,
        observed_freq = observed_freq,
        implied_selectivity = implied_sel,
        length_bins = length_bins,
        bin_mids = bin_mids,
        dome = list(
          params = best_dome$par,
          sse = best_dome$value,
          aic = aic_dome,
          convergence = best_dome$convergence
        ),
        logistic = list(
          params = best_logistic$par,
          sse = best_logistic$value,
          aic = aic_logistic,
          convergence = best_logistic$convergence
        ),
        delta_aic = delta_aic,
        best_model = ifelse(delta_aic < 0, "dome", "logistic")
      )
    }
  }
  
  return(results_list)
}


# 4. CORRECTION FUNCTION
correct_for_selectivity = function(data, selectivity_results, 
                                    model_type = "best",
                                    min_selectivity = 0.05) {
  
  corrected_data_list <- list()
  
  for (yr_name in names(selectivity_results)) {
    
    yr = as.numeric(yr_name)
    res = selectivity_results[[yr_name]]
    
    # Determine which model
    if (model_type == "best") {
      use_model = res$best_model
    } else {
      use_model = model_type
    }
    
    # Get parameters
    if (use_model == "dome") {
      params = res$dome$params
      SL50  = params[1]
      SL_sd = params[2]
      SR50  = params[3]
      SR_sd = params[4]
    } else {
      params = res$logistic$params
      SL50 = params[1]
      SL95 = params[2]
      SR50 = params[3]
      SR95 = params[4]
    }
    
    # Get data for this year
    data_year <- data %>% filter(Ano == yr)
    
    # Calculate selectivity
    if (use_model == "dome") {
      selectivity = dome_selectivity(data_year$LT, SL50, SL_sd, SR50, SR_sd)
    } else {
      selectivity = double_logistic(data_year$LT, SL50, SL95, SR50, SR95)
    }
    
    # Correct frequencies
    selectivity_corrected = pmax(selectivity, min_selectivity)
    corrected_freq = data_year$N / selectivity_corrected
    
    # Create corrected dataset
    corrected_year = data_year %>%
      mutate(
        N_original  = N,
        selectivity = selectivity,
        N           = corrected_freq,
        corrected   = TRUE,
        model_used  = use_model
      )
    
    corrected_data_list[[yr_name]] = corrected_year
  }
  
  corrected_data <- bind_rows(corrected_data_list)
  
  return(corrected_data)
}


# 5. PLOTTING FUNCTIONS
plot_selectivity_results <- function(results, species_name, method = "nonparametric") {
  
  n_years = length(results)
  n_cols  = min(3, n_years)
  n_rows  = ceiling(n_years / n_cols)
  
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
  
  for (yr_name in names(results)) {
    
    res = results[[yr_name]]
    yr  = res$year
    bin_mids = res$bin_mids
    
    # Calculate fitted selectivity curves
    dome_params = res$dome$params
    sel_dome    = dome_selectivity(bin_mids, dome_params[1], dome_params[2],
                                 dome_params[3], dome_params[4])
    
    log_params   = res$logistic$params
    sel_logistic = double_logistic(bin_mids, log_params[1], log_params[2],
                                    log_params[3], log_params[4])
    
    # Plot observed data
    plot(bin_mids, res$observed_freq, type = "h", lwd = 3, col = "gray70",
         xlab = "Length (cm)", ylab = "Frequency",
         main = paste0(species_name, " - ", yr, 
                       "\nBest: ", res$best_model, 
                       " (ΔAIC = ", round(res$delta_aic, 1), ")"),
         ylim = c(0, max(res$observed_freq) * 1.2))
    
    # Add selectivity curves (scaled)
    sel_dome_scaled = sel_dome * max(res$observed_freq)
    sel_log_scaled  = sel_logistic * max(res$observed_freq)
    
    lines(bin_mids, sel_dome_scaled, col = "blue", lwd = 2, lty = 2)
    lines(bin_mids, sel_log_scaled, col = "red", lwd = 2, lty = 2)
    
    # Add expected frequency if available
    if (!is.null(res$expected_freq)) {
      lines(bin_mids, res$expected_freq, col = "green", lwd = 1, lty = 3)
    }
    
    legend("topright", 
           legend = c("Observed", "Dome sel", "Logistic sel", "Expected"),
           col = c("gray70", "blue", "red", "green"),
           lty = c(1, 2, 2, 3), lwd = c(3, 2, 2, 1),
           bty = "n", cex = 0.7)
  }
  
  par(mfrow = c(1, 1))
}

# Summary table
summarize_selectivity = function(results) {
  
  summary_df = map_df(results, function(res) {
    tibble(
      year = res$year,
      n_samples = res$n_samples,
      dome_SL50 = round(res$dome$params[1], 2),
      dome_SR50 = round(res$dome$params[3], 2),
      dome_sd   = round(res$dome$params[2], 2),  
      dome_AIC  = round(res$dome$aic, 2),
      logistic_SL50 = round(res$logistic$params[1], 2),
      logistic_SR50 = round(res$logistic$params[3], 2),
      logistic_AIC  = round(res$logistic$aic, 2),
      delta_AIC     = round(res$delta_aic, 2),
      best_model    = res$best_model
    )
  })
  
  return(summary_df)
}

# ==============================================================================
# 6. USAGE EXAMPLE
# ==============================================================================

# METHOD 1: Non-parametric (mirrors left side)
# results <- estimate_selectivity_nonparametric(
#   data = your_data,
#   species_name = "Sparisoma cretense",
#   bin_width = 1,
#   smooth_span = 0.3,  # Smoothing parameter (0.2-0.5)
#   n_attempts = 10
# )

# METHOD 2: Direct from decline (simpler, but more assumptions)
# results <- estimate_selectivity_from_decline(
#   data = your_data,
#   species_name = "Sparisoma cretense",
#   bin_width = 1,
#   n_attempts = 10
# )

# View results
# summary_table <- summarize_selectivity(results)
# print(summary_table)

# Plot
# plot_selectivity_results(results, "Sparisoma cretense")

# Correct data
# corrected_data <- correct_for_selectivity(data, results, model_type = "best")