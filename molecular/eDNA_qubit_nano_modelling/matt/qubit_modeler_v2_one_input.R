# ==============================================================================
# QUBIT CALIBRATION: Single-file end-to-end script
# - Uses samples WITH Qubit values for calibration
# - Predicts MISSING Qubit values from NanoDrop data
# - LOD handling, robust fit with retries, model selection, CV, diagnostics
# ==============================================================================

# ---- Packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(robustbase)
  library(janitor)
})

# ---- Load single data file ----
cat("Loading data...\n")
data_all <- read_csv("ND_QB_comparison_csv.csv", show_col_types = FALSE) %>% 
  clean_names()

# Check for required columns
required_cols <- c("sample_name", "nano_ng_ul", "a260_a280", "a260_a230", "real_qubit")
missing_cols <- setdiff(required_cols, names(data_all))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Standardize column names
data_all <- data_all %>%
  rename(
    sample_id = sample_name,
    ND_conc = nano_ng_ul,
    a260_280 = a260_a280,
    a260_230 = a260_a230,
    qubit_conc = real_qubit
  )

cat(sprintf("Total samples loaded: %d\n", nrow(data_all)))

# ---- Config ----
CONFIG <- list(
  log_floor    = 0.1,
  winsor_probs = c(0.05, 0.95),
  select_rule  = list(
    cv_rmse_win     = 0.02,
    scale_tie_break = TRUE
  ),
  LOD_Q = 0.05
)

# ---- Split into calibration (has Qubit) and prediction (missing Qubit) ----
cal <- data_all %>% filter(!is.na(qubit_conc))
nd  <- data_all %>% filter(is.na(qubit_conc))

cat(sprintf("Calibration samples (with Qubit): %d\n", nrow(cal)))
cat(sprintf("Prediction samples (missing Qubit): %d\n", nrow(nd)))

if (nrow(cal) < 10) {
  stop("Insufficient calibration samples (need at least 10 with Qubit values)")
}

# ---- Handle below-LOD in calibration ----
cal <- cal %>%
  mutate(is_below_lod = qubit_conc <= CONFIG$LOD_Q)
cal_train <- cal %>% filter(!is_below_lod)

cat(sprintf("Calibration samples used: %d | Excluded below LOD: %d\n",
            nrow(cal_train), sum(cal$is_below_lod)))

# ---- Feature engineering function ----
engineer_features <- function(df, winsor_limits = NULL, is_calibration = TRUE) {
  if (is_calibration) {
    df <- df %>%
      mutate(
        logND = log10(pmax(ND_conc, CONFIG$log_floor)),
        logQ  = log10(pmax(qubit_conc, CONFIG$log_floor)),
        d280  = a260_280 - 2,
        d230  = a260_230 - 2
      )
  } else {
    df <- df %>%
      mutate(
        logND = log10(pmax(ND_conc, CONFIG$log_floor)),
        d280  = a260_280 - 2,
        d230  = a260_230 - 2
      )
  }
  
  df <- df %>% distinct(sample_id, .keep_all = TRUE)
  
  if (!is.null(winsor_limits)) {
    df <- df %>%
      mutate(
        d230_w = pmin(pmax(d230, winsor_limits$d230[1]), winsor_limits$d230[2]),
        d280_w = pmin(pmax(d280, winsor_limits$d280[1]), winsor_limits$d280[2])
      )
    return(df)
  } else {
    limits <- list(
      d230 = quantile(df$d230, CONFIG$winsor_probs, na.rm = TRUE),
      d280 = quantile(df$d280, CONFIG$winsor_probs, na.rm = TRUE)
    )
    df <- df %>%
      mutate(
        d230_w = pmin(pmax(d230, limits$d230[1]), limits$d230[2]),
        d280_w = pmin(pmax(d280, limits$d280[1]), limits$d280[2])
      )
    return(list(data = df, limits = limits))
  }
}

# ---- Engineer features for calibration ----
cat("\nEngineering features...\n")
cal_result    <- engineer_features(cal_train, is_calibration = TRUE)
cal_train     <- cal_result$data
winsor_limits <- cal_result$limits

cat(sprintf("Winsor limits  d230:[%.3f, %.3f]  d280:[%.3f, %.3f]\n",
            winsor_limits$d230[1], winsor_limits$d230[2],
            winsor_limits$d280[1], winsor_limits$d280[2]))

# ---- Remove any rows with missing engineered features ----
cal_train <- cal_train %>% filter(complete.cases(logQ, logND, d230_w, d280_w))
cat(sprintf("Calibration rows after QC: %d\n", nrow(cal_train)))

# ---- Robust fitting function with retries ----
safe_lmrob <- function(formula, data) {
  ctrl1 <- lmrob.control(
    setting = "KS2014",
    maxit.scale = 2000,
    max.it = 400,
    refine.tol = 1e-7,
    trace = FALSE
  )
  fit <- try(lmrob(formula, data = data, control = ctrl1), silent = TRUE)
  if (!inherits(fit, "try-error") && is.finite(fit$scale)) return(fit)
  
  ctrl2 <- lmrob.control(
    setting = "KS2014",
    init = "S",
    maxit.scale = 4000,
    max.it = 800,
    refine.tol = 1e-6,
    trace = FALSE
  )
  fit <- try(lmrob(formula, data = data, control = ctrl2), silent = TRUE)
  if (!inherits(fit, "try-error") && is.finite(fit$scale)) return(fit)
  
  ctrl3 <- lmrob.control(
    setting = "KS2014",
    psi = "huber",
    maxit.scale = 4000,
    max.it = 800,
    refine.tol = 1e-6,
    trace = FALSE
  )
  fit <- lmrob(formula, data = data, control = ctrl3)
  fit
}

# ---- Fit models ----
cat("\n=== Fitting robust models ===\n")
m1_full   <- safe_lmrob(logQ ~ logND + d230_w + d280_w, data = cal_train)
m1_simple <- safe_lmrob(logQ ~ logND + d280_w, data = cal_train)

cat(sprintf("Full model:   R²=%.4f  scale=%.4f\n",
            summary(m1_full)$r.squared, m1_full$scale))
cat(sprintf("Simple model: R²=%.4f  scale=%.4f\n",
            summary(m1_simple)$r.squared, m1_simple$scale))

# ---- K-fold Cross-Validation ----
set.seed(1)
K <- min(5, nrow(cal_train))
fold_id <- sample(rep(1:K, length.out = nrow(cal_train)))

cv_err <- function(mod_formula) {
  cv <- lapply(1:K, function(k){
    trn <- cal_train[fold_id != k, , drop = FALSE]
    tst <- cal_train[fold_id == k, , drop = FALSE]
    fit <- safe_lmrob(mod_formula, trn)
    data.frame(obs = tst$logQ, pred = predict(fit, newdata = tst))
  }) %>% bind_rows()
  
  rmse <- sqrt(mean((cv$obs - cv$pred)^2))
  medfold <- median(10^abs(cv$obs - cv$pred))
  q75fold <- unname(quantile(10^abs(cv$obs - cv$pred), 0.75))
  c(rmse = rmse, med_fold = medfold, q75_fold = q75fold)
}

cv_full   <- cv_err(logQ ~ logND + d230_w + d280_w)
cv_simple <- cv_err(logQ ~ logND + d280_w)

cat(sprintf("\nCV Results (RMSE | Median× | 75th×):\n"))
cat(sprintf("  Full:   %.3f | %.2fx | %.2fx\n",
            cv_full["rmse"], cv_full["med_fold"], cv_full["q75_fold"]))
cat(sprintf("  Simple: %.3f | %.2fx | %.2fx\n",
            cv_simple["rmse"], cv_simple["med_fold"], cv_simple["q75_fold"]))

# ---- Model selection ----
delta <- CONFIG$select_rule$cv_rmse_win
choose_simple <- function() {
  if (abs(cv_full["rmse"] - cv_simple["rmse"]) <= delta) {
    if (CONFIG$select_rule$scale_tie_break) {
      if (m1_simple$scale < m1_full$scale - 1e-9) return(TRUE)
      if (abs(m1_simple$scale - m1_full$scale) <= 1e-6) return(TRUE)
      return(FALSE)
    } else {
      return(TRUE)
    }
  } else {
    return(cv_simple["rmse"] < cv_full["rmse"])
  }
}

use_simple <- choose_simple()
best_model <- if (use_simple) m1_simple else m1_full
best_tag   <- if (use_simple) "simple" else "full"

cat(sprintf("\n=== CHOSEN MODEL: %s ===\n", best_tag))
print(summary(best_model))

# ---- Save calibration diagnostics ----
cal_train$.w <- best_model$rweights
weights_out <- cal_train %>%
  arrange(.w) %>%
  select(sample_id, qubit_conc, ND_conc, a260_280, a260_230, d230_w, d280_w, .w)

write_csv(weights_out, sprintf("calibration_weights_%s.csv", best_tag))

# ---- Predict for samples with missing Qubit values ----
if (nrow(nd) > 0) {
  cat(sprintf("\n=== Predicting %d samples with missing Qubit values ===\n", nrow(nd)))
  
  nd_prep <- engineer_features(nd, winsor_limits = winsor_limits, is_calibration = FALSE)
  
  # Check for samples outside calibration range
  nd_prep <- nd_prep %>%
    mutate(
      outside_range = (
        ND_conc < min(cal_train$ND_conc) | ND_conc > max(cal_train$ND_conc) |
          a260_280 < min(cal_train$a260_280) | a260_280 > max(cal_train$a260_280) |
          a260_230 < min(cal_train$a260_230) | a260_230 > max(cal_train$a260_230)
      )
    )
  
  n_outside <- sum(nd_prep$outside_range, na.rm = TRUE)
  if (n_outside > 0) {
    warning(sprintf("%d samples are outside calibration range - predictions may be unreliable", n_outside))
  }
  
  # Generate predictions with prediction intervals
  pi_nd <- as.data.frame(predict(best_model, newdata = nd_prep, interval = "prediction"))
  names(pi_nd) <- c("logQ_pred", "logQ_lwr", "logQ_upr")
  
  nd_preds <- nd_prep %>%
    bind_cols(pi_nd) %>%
    mutate(
      qubit_pred = pmax(10^logQ_pred, 0),
      qubit_lb   = pmax(10^logQ_lwr, 0),
      qubit_ub   = pmax(10^logQ_upr, 0),
      CI_width   = qubit_ub - qubit_lb,
      ND_to_Qb_correction_factor = ifelse(ND_conc > 0, qubit_pred / ND_conc, NA_real_)
    ) %>%
    select(sample_id, ND_conc, a260_280, a260_230, 
           qubit_pred, qubit_lb, qubit_ub, CI_width,
           ND_to_Qb_correction_factor, outside_range)
  
  write_csv(nd_preds, sprintf("predicted_qubit_values_%s.csv", best_tag))
  cat(sprintf("Predictions saved to: predicted_qubit_values_%s.csv\n", best_tag))
  
  # Summary statistics
  cat("\nPrediction Summary:\n")
  cat(sprintf("  Mean predicted Qubit: %.2f ng/µL\n", mean(nd_preds$qubit_pred, na.rm = TRUE)))
  cat(sprintf("  Median predicted Qubit: %.2f ng/µL\n", median(nd_preds$qubit_pred, na.rm = TRUE)))
  cat(sprintf("  Range: [%.2f, %.2f] ng/µL\n", 
              min(nd_preds$qubit_pred, na.rm = TRUE),
              max(nd_preds$qubit_pred, na.rm = TRUE)))
  cat(sprintf("  Samples outside calibration range: %d\n", n_outside))
  
} else {
  cat("\nNo samples with missing Qubit values - no predictions needed.\n")
}

# ---- Combine calibration + predictions into complete dataset ----
cal_complete <- cal_train %>%
  mutate(
    qubit_pred = qubit_conc,
    qubit_lb = NA_real_,
    qubit_ub = NA_real_,
    CI_width = NA_real_,
    ND_to_Qb_correction_factor = ifelse(ND_conc > 0, qubit_conc / ND_conc, NA_real_),
    outside_range = FALSE,
    data_type = "calibration"
  ) %>%
  select(sample_id, ND_conc, a260_280, a260_230, 
         qubit_conc, qubit_pred, qubit_lb, qubit_ub, CI_width,
         ND_to_Qb_correction_factor, outside_range, data_type)

if (nrow(nd) > 0) {
  pred_complete <- nd_preds %>%
    mutate(
      qubit_conc = NA_real_,
      data_type = "predicted"
    ) %>%
    select(sample_id, ND_conc, a260_280, a260_230,
           qubit_conc, qubit_pred, qubit_lb, qubit_ub, CI_width,
           ND_to_Qb_correction_factor, outside_range, data_type)
  
  complete_dataset <- bind_rows(cal_complete, pred_complete)
} else {
  complete_dataset <- cal_complete
}

write_csv(complete_dataset, sprintf("complete_dataset_with_predictions_%s.csv", best_tag))
cat(sprintf("\nComplete dataset saved to: complete_dataset_with_predictions_%s.csv\n", best_tag))

# ---- Diagnostic plots ----
cal_diag <- cal_train %>%
  mutate(
    .fitted = predict(best_model, newdata = cal_train),
    .resid  = logQ - .fitted,
    Q_obs   = qubit_conc,
    Q_pred  = 10^.fitted
  )

p1 <- ggplot(cal_diag, aes(x = Q_obs, y = Q_pred, color = .w)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_gradient(low = "red", high = "blue", name = "Robust\nWeight") +
  coord_equal() +
  labs(
    x = "Observed Qubit (ng/µL)",
    y = "Predicted Qubit (ng/µL)",
    title = sprintf("Calibration: Observed vs Predicted (%s model)", best_tag),
    subtitle = sprintf("R² = %.3f | n = %d", summary(best_model)$r.squared, nrow(cal_diag))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
p1

p2 <- ggplot(cal_diag, aes(x = .fitted, y = .resid, color = .w)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_gradient(low = "red", high = "blue", name = "Robust\nWeight") +
  labs(
    x = "Fitted log10(Qubit)",
    y = "Residual (log10 scale)",
    title = sprintf("Residual Plot (%s model)", best_tag)
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())


print(p1)
print(p2)

ggsave(sprintf("calibration_obs_vs_pred_%s.png", best_tag), p1, 
       width = 7, height = 5.5, dpi = 300)
ggsave(sprintf("calibration_residuals_%s.png", best_tag), p2, 
       width = 7, height = 5.5, dpi = 300)

# ---- Performance metrics ----
res_log <- cal_diag$.resid
GMFE  <- 10^mean(abs(res_log), na.rm = TRUE)
MedFE <- 10^median(abs(res_log), na.rm = TRUE)

cat("\n=== Model Performance ===\n")
cat(sprintf("In-sample Geometric Mean Fold Error: %.2fx\n", GMFE))
cat(sprintf("In-sample Median Fold Error: %.2fx\n", MedFE))
cat(sprintf("CV RMSE (log10): %.3f\n", ifelse(use_simple, cv_simple["rmse"], cv_full["rmse"])))
cat(sprintf("CV Median Fold Error: %.2fx\n", ifelse(use_simple, cv_simple["med_fold"], cv_full["med_fold"])))

# ---- Domain of validity ----
domain <- cal_train %>%
  summarise(
    ND_min = min(ND_conc, na.rm = TRUE),
    ND_max = max(ND_conc, na.rm = TRUE),
    Q_min = min(qubit_conc, na.rm = TRUE),
    Q_max = max(qubit_conc, na.rm = TRUE),
    a280_min = min(a260_280, na.rm = TRUE),
    a280_max = max(a260_280, na.rm = TRUE),
    a230_min = min(a260_230, na.rm = TRUE),
    a230_max = max(a260_230, na.rm = TRUE)
  )

cat("\n=== Calibration Domain (do not extrapolate beyond) ===\n")
print(domain)
write_csv(domain, sprintf("calibration_domain_%s.csv", best_tag))

cat("\n=== Analysis Complete ===\n")

