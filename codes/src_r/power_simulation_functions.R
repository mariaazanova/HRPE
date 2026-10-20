library(brms)
options(brms.backend = "cmdstanr")
library(posterior)
library(future)
library(future.apply)

set_seed <- 423
set.seed(set_seed)

.pick_draw_id <- function(fit, mode = c("single_draw", "closest_to_mean"), seed = NULL) {
  mode <- match.arg(mode)
  dr <- posterior::as_draws_df(fit)
  b_cols <- grep("^b_", colnames(dr), value = TRUE)
  if (mode == "single_draw") {
    if (!is.null(seed)) set.seed(seed)
    return(sample.int(nrow(dr), 1))
  } else {
    b_means <- colMeans(dr[, b_cols, drop = FALSE])
    # Euclidean distance on fixed effects only (stable, fast)
    d <- as.matrix(dr[, b_cols, drop = FALSE])
    dd <- sqrt(rowSums((t(t(d) - b_means))^2))
    return(which.min(dd))
  }
}

.build_newdata <- function(fit, pilot_df, n_sub, n_trials, id_var = "subID") {
  # Model formula and variables
  ff <- fit$formula$formula
  resp <- all.vars(ff[[2]])[1]
  rhs_vars <- setdiff(all.vars(ff), resp)
  
  # Keep only variables that exist in pilot_df
  rhs_vars <- rhs_vars[rhs_vars %in% names(pilot_df)]
  need_vars <- unique(c(rhs_vars, id_var))
  need_vars <- need_vars[need_vars %in% names(pilot_df)]
  
  # Complete cases on required predictors
  base_df <- pilot_df[, need_vars, drop = FALSE]
  base_df <- base_df[stats::complete.cases(base_df), , drop = FALSE]
  stopifnot(nrow(base_df) > 0)
  
  # Preserve factor levels exactly as in pilot data
  for (v in names(base_df)) {
    if (is.factor(pilot_df[[v]])) {
      base_df[[v]] <- factor(as.character(base_df[[v]]),
                             levels = levels(pilot_df[[v]]))
    }
  }
  
  # Create new subject IDs
  new_ids <- sprintf("sim_sub_%03d", seq_len(n_sub))
  
  # Sample n_trials rows per subject from base_df (predictors only; outcomes ignored)
  by_sub <- lapply(seq_len(n_sub), function(i) {
    idx <- sample.int(nrow(base_df), n_trials, replace = TRUE)
    nd  <- base_df[idx, , drop = FALSE]
    # assign NEW levels to grouping factor
    nd[[id_var]] <- factor(new_ids[i])
    nd
  })
  
  nd_all <- do.call(rbind, by_sub)
  
  # Ensure id_var is a factor and allows new levels at predict-time
  nd_all[[id_var]] <- factor(nd_all[[id_var]])
  
  list(newdata = nd_all, resp = resp)
}

simulate_power_BF_pp <- function(
    brms_fit, pilot_df,
    n_sim = 100, n_sub = 100, n_trials = 300,
    refit_iter = 2000, chains = 2, cores = 2,
    BF_threshold = 10,
    hypothesis_direction = c("equal", "greater", "less"),
    param_source = c("single_draw", "closest_to_mean"),
    base_seed = 1234,
    diagnostics_gate = TRUE,
    adapt_delta = 0.95, max_treedepth = 12,
    contrasts_spec = NULL,
    include_intercept = FALSE,
    rhat_threshold = 1.01,
    gating_scope = c("global", "per_term"),
    eps = 1e-12
) {
  stopifnot(inherits(brms_fit, "brmsfit"))
  hypothesis_direction <- match.arg(hypothesis_direction)
  param_source <- match.arg(param_source)
  gating_scope <- match.arg(gating_scope)
  
  old_contr <- getOption("contrasts")
  if (!is.null(contrasts_spec)) {
    options(contrasts = contrasts_spec)
    on.exit(options(contrasts = old_contr), add = TRUE)
  }
  
  # Fixed effects to test (optionally drop Intercept)
  fixef_all   <- rownames(fixef(brms_fit))
  fixef_names <- if (include_intercept) fixef_all else setdiff(fixef_all, "Intercept")
  K <- length(fixef_names)
  if (K == 0) stop("No fixed effects to test (did you exclude the Intercept?).")
  
  # Parallel
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = cores)
  
  # Collect per-sim results (no gating here)
  sim_list <- future.apply::future_lapply(
    seq_len(n_sim),
    function(i_sim) {
      set.seed(base_seed + i_sim)
      
      nd_info <- .build_newdata(brms_fit, pilot_df, n_sub = n_sub, n_trials = n_trials, id_var = "subID")
      newdata <- nd_info$newdata
      resp    <- nd_info$resp
      
      draw_id <- .pick_draw_id(brms_fit, mode = param_source, seed = base_seed + 100000 + i_sim)
      
      has_draw_ids_arg <- "draw_ids" %in% names(formals(brms::posterior_predict))
      if (has_draw_ids_arg) {
        ysim <- brms::posterior_predict(
          brms_fit, newdata = newdata, allow_new_levels = TRUE,
          sample_new_levels = "gaussian", re_formula = NULL,
          draws = 1, draw_ids = draw_id
        )
      } else {
        ysim <- brms::posterior_predict(
          brms_fit, newdata = newdata, allow_new_levels = TRUE,
          sample_new_levels = "gaussian", re_formula = NULL,
          draws = 1
        )
      }
      newdata[[resp]] <- as.numeric(ysim[1, ])
      
      sim_fit <- update(
        brms_fit,
        newdata = newdata,
        chains = chains, iter = refit_iter, cores = cores, refresh = 0,
        control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
      )
      
      # Diagnostics
      rhat_vals <- tryCatch(brms::rhat(sim_fit), error = function(e) NA_real_)
      # Global max R-hat across *all* parameters (conservative)
      rhat_max_all <- tryCatch({
        if (is.matrix(rhat_vals)) max(rhat_vals[, 1], na.rm = TRUE) else max(rhat_vals, na.rm = TRUE)
      }, error = function(e) NA_real_)
      
      # Per-term R-hat for fixed effects
      if (is.matrix(rhat_vals)) rhat_vals <- rhat_vals[, 1]
      b_rhat <- rhat_vals[paste0("b_", fixef_names)]
      # if names missing, keep NA
      
      # Divergences (global flag)
      np <- tryCatch(brms::nuts_params(sim_fit), error = function(e) NULL)
      div_flag <- FALSE
      if (!is.null(np)) {
        divs <- tryCatch(do.call(rbind, np), error = function(e) NULL)
        if (!is.null(divs) && "divergent__" %in% names(divs)) {
          div_flag <- sum(divs$divergent__ > 0, na.rm = TRUE) > 0
        }
      }
      
      # Evidence ratios for requested terms
      BF_raw <- numeric(K)
      for (j in seq_len(K)) {
        term <- fixef_names[j]
        hstr <- switch(
          hypothesis_direction,
          "equal"   = paste0(term, " = 0"),
          "greater" = paste0(term, " > 0"),
          "less"    = paste0(term, " < 0")
        )
        hres <- hypothesis(sim_fit, hstr)
        BF_raw[j] <- as.numeric(hres$hypothesis$Evid.Ratio)
      }
      
      # Sanitize tiny/invalid negatives
      BF_raw[is.nan(BF_raw)] <- NA_real_
      BF_raw[is.infinite(BF_raw) & BF_raw < 0] <- NA_real_
      BF_raw[BF_raw < 0 & BF_raw > -eps] <- 0
      BF_raw[BF_raw < 0] <- NA_real_
      
      names(BF_raw) <- fixef_names
      
      # Final orientation ("> BF_threshold" ⇒ evidence for effect/H1)
      BF_final <- if (hypothesis_direction == "equal") 1 / BF_raw else BF_raw
      BF_final[is.nan(BF_final)] <- NA_real_
      BF_final[is.infinite(BF_final) & BF_final < 0] <- NA_real_
      BF_final[BF_final < 0 & BF_final > -eps] <- 0
      BF_final[BF_final < 0] <- NA_real_
      names(BF_final) <- fixef_names
      
      list(
        bf_raw     = BF_raw,
        bf_final   = BF_final,
        b_rhat     = b_rhat,
        rhat_max   = rhat_max_all,
        divergent  = div_flag
      )
    },
    future.seed = TRUE
  )
  
  # Bind "ALL" (no gating) matrices
  BF_raw_all   <- do.call(rbind, lapply(sim_list, `[[`, "bf_raw"))
  BF_final_all <- do.call(rbind, lapply(sim_list, `[[`, "bf_final"))
  Rhat_b_mat   <- do.call(rbind, lapply(sim_list, `[[`, "b_rhat"))
  rhat_max_vec <- vapply(sim_list, `[[`, numeric(1), "rhat_max")
  div_vec      <- vapply(sim_list, `[[`, logical(1), "divergent")
  
  colnames(BF_raw_all)   <- fixef_names
  colnames(BF_final_all) <- fixef_names
  colnames(Rhat_b_mat)   <- fixef_names
  rownames(BF_raw_all)   <- rownames(BF_final_all) <- rownames(Rhat_b_mat) <- NULL
  
  # ---- Apply diagnostics gating to produce the usual (gated) outputs ----
  BF_raw_gated   <- BF_raw_all
  BF_final_gated <- BF_final_all
  
  if (diagnostics_gate) {
    if (gating_scope == "global") {
      bad_row <- div_vec | !is.finite(rhat_max_vec) | (rhat_max_vec >= rhat_threshold)
      BF_raw_gated[bad_row, ]   <- NA_real_
      BF_final_gated[bad_row, ] <- NA_real_
    } else { # per_term
      # mask cells whose *own* fixed-effect R-hat is high; fallback to global if all NAs
      per_term_bad <- !is.na(Rhat_b_mat) & Rhat_b_mat >= rhat_threshold
      all_na_row   <- apply(Rhat_b_mat, 1, function(x) all(is.na(x)))
      # Per-cell mask
      BF_raw_gated[per_term_bad]   <- NA_real_
      BF_final_gated[per_term_bad] <- NA_real_
      # Fallback global mask for rows with no per-term R-hat info
      if (any(all_na_row)) {
        fallback_bad <- (div_vec | !is.finite(rhat_max_vec) | (rhat_max_vec >= rhat_threshold)) & all_na_row
        BF_raw_gated[fallback_bad, ]   <- NA_real_
        BF_final_gated[fallback_bad, ] <- NA_real_
      }
      # Divergences always nuke the whole row (conservative)
      if (any(div_vec)) {
        BF_raw_gated[div_vec, ]   <- NA_real_
        BF_final_gated[div_vec, ] <- NA_real_
      }
    }
  }
  
  # Validity message: rows fully OK after gating
  fully_valid <- rowSums(is.na(BF_final_gated)) == 0
  message(sprintf("Valid sims (all terms OK): %d / %d (%.1f%%).",
                  sum(fully_valid), n_sim, 100 * mean(fully_valid)))
  
  # Summaries from *gated* FINAL matrix (matches manual checks)
  evidence_effect <- setNames(colMeans(BF_final_gated >  BF_threshold,    na.rm = TRUE), fixef_names)
  evidence_null   <- setNames(colMeans(BF_final_gated < (1 / BF_threshold), na.rm = TRUE), fixef_names)
  
  list(
    # GATED 
    BF_matrix_raw    = as.data.frame(BF_raw_gated,   check.names = FALSE),
    BF_matrix_final  = as.data.frame(BF_final_gated, check.names = FALSE),
    evidence_effect  = as.data.frame(as.list(evidence_effect), check.names = FALSE),
    evidence_null    = as.data.frame(as.list(evidence_null),   check.names = FALSE),
    
    # NEW: ALL (no gating at all)
    BF_matrix_raw_all   = as.data.frame(BF_raw_all,   check.names = FALSE),
    BF_matrix_final_all = as.data.frame(BF_final_all, check.names = FALSE),
    
    # NEW: Diagnostics matrices/vectors 
    Rhat_matrix_fixed = as.data.frame(Rhat_b_mat, check.names = FALSE),  # R-hat for tested fixed effects
    Rhat_max_all      = rhat_max_vec,                                     # global max R-hat per sim
    Divergences       = div_vec,
    
    settings = list(
      n_sim = n_sim, n_sub = n_sub, n_trials = n_trials,
      BF_threshold = BF_threshold,
      hypothesis_direction = hypothesis_direction,
      param_source = param_source,
      diagnostics_gate = diagnostics_gate,
      adapt_delta = adapt_delta, max_treedepth = max_treedepth,
      include_intercept = include_intercept,
      rhat_threshold = rhat_threshold,
      gating_scope = gating_scope,
      eps = eps
    )
  )
}

evidence_from_bf_rhat <- function(
    BF_final_mat,                 # matrix/data.frame, "> BF_thr" => supports effect/H1
    Rhat_mat = NULL,              # per-term R-hat (same shape as BF_final_mat), optional for global
    BF_threshold = 10,
    rhat_threshold = 1.01,
    gating_scope = c("none", "global", "per_term"),
    rhat_max_vec = NULL,          # optional: global max R-hat per sim (vector length nrow)
    divergences = NULL            # optional: logical vector length nrow; TRUE => drop row
) {
  gating_scope <- match.arg(gating_scope)
  B <- as.matrix(BF_final_mat)
  nr <- nrow(B); nc <- ncol(B)
  
  if (!is.null(divergences)) {
    B[divergences, ] <- NA_real_
  }
  
  if (gating_scope == "per_term") {
    stopifnot(!is.null(Rhat_mat))
    R <- as.matrix(Rhat_mat)
    stopifnot(nrow(R) == nr && ncol(R) == nc)
    bad <- !is.na(R) & (R >= rhat_threshold)
    B[bad] <- NA_real_
  } else if (gating_scope == "global") {
    if (!is.null(rhat_max_vec)) {
      bad_row <- !is.finite(rhat_max_vec) | (rhat_max_vec >= rhat_threshold)
    } else if (!is.null(Rhat_mat)) {
      R <- as.matrix(Rhat_mat)
      row_max <- apply(R, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
      bad_row <- !is.finite(row_max) | (row_max >= rhat_threshold)
    } else {
      stop("For global gating, provide either rhat_max_vec or Rhat_mat.")
    }
    B[bad_row, ] <- NA_real_
  } # else "none": no extra masking
  
  evidence_effect <- colMeans(B >  BF_threshold,    na.rm = TRUE)
  evidence_null   <- colMeans(B < (1 / BF_threshold), na.rm = TRUE)
  
  list(
    evidence_effect = as.data.frame(as.list(evidence_effect), check.names = FALSE),
    evidence_null   = as.data.frame(as.list(evidence_null),   check.names = FALSE),
    valid_rows      = rowSums(is.na(B)) == 0
  )
}

