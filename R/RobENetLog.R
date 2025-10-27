robust_elasticnet_logistic = function(X, y,
                                      alpha_seq = seq(0.1, 0.9, length.out=41),
                                      h_prop_csteps = 0.75,
                                      num_csteps = 20,
                                      nfolds = 5,
                                      standardize =  TRUE,
                                      alpha_wt_seq = seq(0, 1, length.out=41),
                                      initial_subset_cutoff_factor = 3,
                                      reweight_cutoff_factor = 3
                                      )
    {

    initial_subset = init_subset(X = X,
                                standardize = standardize,
                                initial_subset_cutoff_factor = initial_subset_cutoff_factor)


    h_n = length(initial_subset$subset_init)
    #h_n = floor(h_prop_csteps * n)
    if (h_n < 2) stop("h_prop_csteps is too small, resulting in less than 2 samples for initial subset.")
    #subset_init = order(score_x)[1:h_n]

    y_subset_init = y[initial_subset$subset_init]
    if (length(unique(y_subset_init)) < 2) {
        stop(paste("Initial subset for alpha (glmnet) tuning has only one class. Samples:", length(y_subset_init),
                   "Unique Y:", paste(unique(y_subset_init), collapse=","),
                   "Smallest score_x:", paste(head(sort(initial_subset$score_x)), collapse=", ")))
    }

    nfolds_init_cv = nfolds
    if (length(y_subset_init) < nfolds) {
        warning(paste("Number of samples in initial subset (",length(y_subset_init),") is less than nfolds (",nfolds,"). Reducing nfolds_init_cv to ", length(y_subset_init)))
        nfolds_init_cv = max(2, length(y_subset_init))
    }

    pre_tune = .pre_tuning(alpha_seq = alpha_seq,
                           y_subset_init = y_subset_init,
                           nfolds_init_cv = nfolds_init_cv,
                           subset_init = initial_subset$subset_init,
                           X_scaled = initial_subset$X_scaled,
                           subset_X = initial_subset$subset_X)


    overall_best_model = NULL
    overall_best_subset_c = NULL
    min_total_deviance_on_full_data = Inf
    tuned_alpha_wt = ifelse(length(alpha_wt_seq) > 0, alpha_wt_seq[1], 0.5)

    # initialize list to store convergence history ###
    convergence_history_list <- list()

    for (current_alpha_wt in alpha_wt_seq) {
        if (length(unique(y[initial_subset$subset_init])) < 2) {
            warning(sprintf("Skipping alpha_wt=%.2f: Initial subset has only one class.", current_alpha_wt))
            next
        }
        current_cstep_model = tryCatch(
            glmnet(initial_subset$subset_X, y[initial_subset$subset_init], alpha = pre_tune$best_alpha_glmnet,
                   lambda = pre_tune$best_lambda_for_csteps, family = "binomial", standardize=FALSE),
            error = function(e){
                warning(sprintf("Initial glmnet fit failed for C-steps with alpha_wt=%.2f: %s", current_alpha_wt, e$message))
                return(NULL)
            }
        )
        if(is.null(current_cstep_model)) next

        # Initialize history for this alpha_wt ###
        history_this_awt <- data.frame()

        subset_prev_c = NULL
        temp_subset_c = initial_subset$subset_init
        for (i_cstep in 1:num_csteps) {

            #calculate and store deviance on current subset ###
            preds_on_subset <- predict(current_cstep_model, newx = initial_subset$X_scaled[temp_subset_c, , drop=FALSE], type="response", s=pre_tune$best_lambda_for_csteps)
            preds_on_subset <- pmax(pmin(as.vector(preds_on_subset), 1 - 1e-8), 1e-8)
            y_temp_subset_c_for_dev <- y[temp_subset_c]
            #deviance_on_subset <- -2 * sum(y_temp_subset_c_for_dev * log(preds_on_subset) + (1-y_temp_subset_c_for_dev)*log(1-preds_on_subset), na.rm=TRUE)
            deviance_on_subset = sum(deviance_res(y_temp_subset_c_for_dev,preds_on_subset),na.rm = TRUE)

            history_this_awt <- rbind(history_this_awt, data.frame(
                alpha_wt = current_alpha_wt,
                step = i_cstep,
                subset_deviance = deviance_on_subset
            ))

            pred_probs_full = tryCatch(
                predict(current_cstep_model, newx = initial_subset$X_scaled, type = "response", s = pre_tune$best_lambda_for_csteps),
                error = function(e) { warning(sprintf("Predict failed in C-step %d for alpha_wt %.2f", i_cstep, current_alpha_wt)); NULL}
            )
            if(is.null(pred_probs_full)) break

            pred_probs_full = pmax(pmin(as.vector(pred_probs_full), 1 - 1e-8), 1e-8)

            #dev_resid_values_full = -2 * (y * log(pred_probs_full) + (1 - y) * log(1 - pred_probs_full))
            dev_resid_values_full = deviance_res(y,pred_probs_full)


            scaled_score_x = min_max_scale(score_x)

            scaled_dev_resid = min_max_scale(abs(dev_resid_values_full))

            score_combined = current_alpha_wt * scaled_score_x + (1 - current_alpha_wt) * scaled_dev_resid


            temp_subset_c = order(score_combined)[1:floor(h_prop_csteps * n)]#[1:h_n] # it is dynamic currently, try with fixed 75% as suggested in c-step classical methods

            if (!is.null(subset_prev_c) && all(sort(temp_subset_c) == sort(subset_prev_c))) {
                break
            } # check for convergence
            subset_prev_c = temp_subset_c

            y_temp_subset_c = y[temp_subset_c]
            if (length(temp_subset_c) < 2 || length(unique(y_temp_subset_c)) < 2) {
                warning(sprintf("Subset for C-step for alpha_wt=%.2f became too small or univariant at C-step %d. Breaking C-steps.", current_alpha_wt, i_cstep))
                break
            }
            current_cstep_model = tryCatch(
                glmnet(initial_subset$X_scaled[temp_subset_c, ], y_temp_subset_c, alpha = pre_tune$best_alpha_glmnet,
                       lambda = pre_tune$best_lambda_for_csteps, family = "binomial", standardize=FALSE),
                error = function(e){
                    warning(sprintf("glmnet fit failed during C-step %d for alpha_wt=%.2f: %s", i_cstep, current_alpha_wt, e$message))
                    return(NULL)
                }
            )
            if(is.null(current_cstep_model)) break
        } # end of i_cstep loop

        #store the convergence history for this alpha_wt ###
        convergence_history_list[[as.character(current_alpha_wt)]] <- history_this_awt

        if(is.null(current_cstep_model)) {
            warning(sprintf("C-steps failed to produce a final model for alpha_wt = %.2f. Skipping.", current_alpha_wt))
            next
        }

        final_pred_probs_for_this_alpha_wt = tryCatch(
            predict(current_cstep_model, newx = initial_subset$X_scaled, type = "response", s = pre_tune$best_lambda_for_csteps),
            error = function(e) { warning(sprintf("Final predict failed for alpha_wt %.2f", current_alpha_wt)); NULL}
        )
        if(is.null(final_pred_probs_for_this_alpha_wt)) next

        final_pred_probs_for_this_alpha_wt = pmax(pmin(as.vector(final_pred_probs_for_this_alpha_wt), 1 - 1e-8), 1e-8)

        #current_total_deviance_full = -2 * sum(y * log(final_pred_probs_for_this_alpha_wt) + (1 - y) * log(1 - final_pred_probs_for_this_alpha_wt), na.rm = TRUE)
        current_total_deviance_full = sum(deviance_res(y,final_pred_probs_for_this_alpha_wt),na.rm = TRUE)

        if (is.finite(current_total_deviance_full) && current_total_deviance_full < min_total_deviance_on_full_data) {
            min_total_deviance_on_full_data = current_total_deviance_full
            tuned_alpha_wt = current_alpha_wt
            overall_best_model = current_cstep_model
            overall_best_subset_c = temp_subset_c
        }
    } # end of current_alpha_wt loop

    if (is.null(overall_best_model)) {
        # Fallback: if no alpha_wt worked, try to use the initial model from alpha tuning
        warning("Failed to find any suitable model after tuning alpha_wt. Trying to use initial model if possible.")
        # Check if initial model based on best_alpha_glmnet can be refit or used.
        # This part would require a more robust fallback strategy.
        # For now, we will error out if no model is found.
        if(exists("pre_tune$best_alpha_glmnet") && exists("pre_tune$best_lambda_for_csteps") && length(unique(y[initial_subset$subset_init])) >= 2){
            overall_best_model <- tryCatch(
                glmnet(initial_subset$subset_X, y[initial_subset$subset_init], alpha = pre_tune$best_alpha_glmnet,
                       lambda =  pre_tune$best_lambda_for_csteps, family = "binomial", standardize=FALSE),
                error = function(e) NULL
            )
            if(!is.null(overall_best_model)){
                tuned_alpha_wt <- NA # Indicate that alpha_wt tuning failed
                overall_best_subset_c <-  initial_subset$subset_init
                warning("Fell back to model based on initial subset (no C-steps or alpha_wt tuning).")

                # Recalculate deviance for this fallback model
                final_pred_probs_fallback = predict(overall_best_model, newx = initial_subset$X_scaled, type = "response", s = pre_tune$best_lambda_for_csteps)
                final_pred_probs_fallback = pmax(pmin(as.vector(final_pred_probs_fallback), 1 - 1e-8), 1e-8)
                #min_total_deviance_on_full_data = -2 * sum(y * log(final_pred_probs_fallback) + (1 - y) * log(1 - final_pred_probs_fallback), na.rm = TRUE)
                min_total_deviance_on_full_data = sum(deviance_res(y,final_pred_probs_fallback))

            } else {
                stop("Failed to find any suitable model after tuning alpha_wt, and fallback failed.")
            }
        } else {
            stop("Failed to find any suitable model after tuning alpha_wt, and initial parameters for fallback are missing.")
        }
    }


    # reweighting step
    reweight_threshold_quantile <- 0.95 # e.g., consider top 2.5% of residuals as outliers for reweighting
    # This is often chosen based on chi-squared quantiles for residuals.
    # For logistic regression, deviance residuals don't have a simple chi-squared.
    # A common choice is a cutoff like 2 or 2.5 for standardized residuals,
    # or a quantile of the robustly estimated residuals.

    # Final predictions and scores are based on the 'overall_best_model' before rewiightsing
    final_pred_probs_before_reweight = predict(overall_best_model, newx = initial_subset$X_scaled, type = "response", s = pre_tune$best_lambda_for_csteps)
    final_pred_probs_before_reweight = pmax(pmin(as.vector(final_pred_probs_before_reweight), 1 - 1e-8), 1e-8)
    raw_beta = overall_best_model$beta
    raw_model = overall_best_model

    # pearson resuidals

    #
    p_robust_from_predict <- predict(overall_best_model, newx = initial_subset$X_scaled, type = "response", s = pre_tune$best_lambda_for_csteps)
    p_robust <- pmax(pmin(as.vector(p_robust_from_predict), 1 - 1e-8), 1e-8)

    final_dev_resids_before_reweight = ( (y - p_robust) / sqrt(p_robust * (1-p_robust) + 1e-8) ) # Added epsilon for stability
    #Alternative for residuals: Pearson residuals, standardized Pearson, or standardized deviance residuals
    # For simplicity, let's use the absolute deviance residuals.
    #print(final_dev_resids_before_reweight);print((final_dev_resids_before_reweight)^2)
    #abs_dev_resids <- abs(final_dev_resids_before_reweight) # Deviance residuals are non-negative
    #print(abs_dev_resids)
    # Determine the cutoff for "good" observations based on these residuals
    #    This is a crucial step. The paper likely uses a chi-squared quantile for
    #    standardized residuals in linear models. For logistic, it's less direct.
    #    A common robust approach is to use a quantile of these residuals.
    #    For example, if we expect (1-h_prop_csteps) outliers, we can use that.
    #    Or a fixed quantile like 0.975 for a 2.5% outlier threshold.
    #cutoff_resids <- quantile(abs_dev_resids, probs = reweight_threshold_quantile, na.rm = TRUE)

    abs_residuals <- abs(as.vector(final_dev_resids_before_reweight)) # Ensure it's a vector

    median_abs_res <- median(abs_residuals, na.rm = TRUE)
    mad_abs_res <- mad(abs_residuals, center = median_abs_res, constant = 1.4826, na.rm = TRUE)

    if (is.na(mad_abs_res) || mad_abs_res < .Machine$double.eps) { # If MAD is zero
        warning("MAD of absolute residuals is zero. Falling back to quantile for reweighting.")
        cutoff_resids_dynamic <- quantile(abs_residuals, probs = 0.95, na.rm = TRUE)
    } else {
        # Outliers are those e.g., > median + 2.5 * MAD (or 3 * MAD)
        # This is analogous to z-score > 2.5 or 3 for normally distributed data
        reweight_cutoff_factor = reweight_cutoff_factor
        cutoff_resids_dynamic <- median_abs_res + reweight_cutoff_factor * mad_abs_res
        cat(paste("Dynamic cutoff (MAD based):", cutoff_resids_dynamic, "\n"))
    }

    observation_weights_reweight <- ifelse(abs_residuals <= cutoff_resids_dynamic, 1, 0)
    num_reweighted_good_obs <- sum(observation_weights_reweight)
    cat(paste("Dynamic Reweighting: num_good_obs =", num_reweighted_good_obs, "\n"))

    y_good_obs <- y[observation_weights_reweight == 1] # Subset of y for good observations

    # Determine nfolds for reweighted CV adaptively
    nfolds_reweight_cv <- nfolds # Start with the main nfolds parameter

    if (num_reweighted_good_obs > 0 && length(unique(y_good_obs)) >= 2) {
        min_class_count_reweight <- min(table(factor(y_good_obs, levels=c(0,1))))

        # Adjust nfolds if the number of good observations or class counts are too small
        if (num_reweighted_good_obs < nfolds_reweight_cv || min_class_count_reweight < nfolds_reweight_cv) {
            nfolds_reweight_cv = max(2, min_class_count_reweight) # At least 2 folds, or based on min class
            if (num_reweighted_good_obs < nfolds_reweight_cv && num_reweighted_good_obs > 0) {
                nfolds_reweight_cv = num_reweighted_good_obs # Fallback to LOOCV for the good obs
            }
        }
        # Ensure nfolds_reweight_cv is at least 2 if there are at least 2 good observations
        if (num_reweighted_good_obs >= 2 && nfolds_reweight_cv < 2) {
            nfolds_reweight_cv <- 2
        }
        # cv.glmnet requires nfolds >= 3 unless N < 3 (then nfolds=N for LOOCV).
        # More robustly: if num_reweighted_good_obs < 3, nfolds_reweight_cv becomes num_reweighted_good_obs for LOOCV.
        if (num_reweighted_good_obs < 3 && num_reweighted_good_obs > 0) {
            nfolds_reweight_cv <- num_reweighted_good_obs
        }


    } else {
        # Not enough data or classes to even attempt determining nfolds_reweight_cv
        nfolds_reweight_cv = 0 # Will lead to can_reweight = FALSE
    }


    # --- Check if reweighting is feasible ---
    can_reweight <- TRUE
    warning_message <- ""

    if (num_reweighted_good_obs < 2) { # Need at least 2 observations for any fit
        can_reweight <- FALSE
        warning_message <- "Reweighting step: Less than 2 'good' observations."
    } else if (length(unique(y_good_obs)) < 2) { # Need both classes
        can_reweight <- FALSE
        warning_message <- "Reweighting step: Only one class present in the 'good' subset."
    } else if (nfolds_reweight_cv < 2 && num_reweighted_good_obs >=2) { # This should ideally not be hit if logic above is right
        # but as a safeguard if nfolds_reweight_cv calculation leads to <2
        can_reweight <- FALSE
        warning_message <- paste("Reweighting step: Determined nfolds_reweight_cv =", nfolds_reweight_cv, "which is too small for CV. Need at least 2.")
    }
    # Note: The condition `num_reweighted_good_obs < ncol(X_scaled)` is NOT used here.

    if (!can_reweight) {
        warning(paste(warning_message, "Skipping reweighting, using previous best model."))
        # overall_best_model from C-steps remains the final model.
        # (No further code in this 'if' block is needed for the reweighting part itself)
    } else {
        # Proceed with reweighting using cv.glmnet and glmnet
        cat(paste("Proceeding with reweighting. num_good_obs =", num_reweighted_good_obs, "nfolds_reweight_cv =", nfolds_reweight_cv, "\n"))

        foldid_reweight <- stratified_fold(y_good_obs, nfolds_reweight_cv) # Pass y_good_obs

        cv_reweighted_fit <- tryCatch(
            cv.glmnet(initial_subset$X_scaled[observation_weights_reweight == 1, , drop=FALSE], # Subset X_scaled
                      y_good_obs,                                                # Use y_good_obs
                      alpha = pre_tune$best_alpha_glmnet,
                      family = "binomial",
                      nfolds = nfolds_reweight_cv, # Use the adaptively determined nfolds
                      foldid = foldid_reweight,
                      type.measure = "deviance",   # Or "auc"
                      standardize = FALSE),        # X_scaled is already scaled
            error = function(e) {
                warning(paste("cv.glmnet for reweighting failed:", e$message, ". Using previous best model."))
                return(NULL)
            }
        )

        if (!is.null(cv_reweighted_fit)) {
            reweighted_lambda <- cv_reweighted_fit$lambda.1se
            if(is.null(reweighted_lambda) || !is.finite(reweighted_lambda)) {
                reweighted_lambda <- cv_reweighted_fit$lambda.min
            }
            if(is.null(reweighted_lambda) || !is.finite(reweighted_lambda)){ # If still problematic
                warning("Could not determine a valid lambda from reweighted CV. Using previous best model.")
                # overall_best_model remains unchanged
            } else {
                reweighted_model_final <- tryCatch(
                    glmnet(initial_subset$X_scaled[observation_weights_reweight == 1, , drop=FALSE],
                           y_good_obs,
                           alpha = pre_tune$best_alpha_glmnet,
                           lambda = reweighted_lambda,
                           family = "binomial",
                           standardize = FALSE),
                    error = function(e) {
                        warning(paste("Final glmnet for reweighting failed:", e$message, ". Using previous best model."))
                        return(NULL)
                    }
                )

                if (!is.null(reweighted_model_final)) {
                    cat("Reweighting successful. Updating final model.\n")
                    # IMPORTANT: Update the variables that define the final model
                    overall_best_model <- reweighted_model_final
                    best_lambda_for_csteps <- reweighted_lambda # Update this to reflect the lambda used by the new overall_best_model
                    overall_best_subset_c <- which(observation_weights_reweight == 1)
                }
                # If reweighted_model_final is NULL, overall_best_model remains the one from C-steps
            }
        }
        # If cv_reweighted_fit is NULL, overall_best_model remains the one from C-steps
    }

    #final calculations
    final_pred_probs_on_full_data = predict(overall_best_model, newx = initial_subset$X_scaled, type = "response",s = best_lambda_for_csteps)
    final_pred_probs_on_full_data = pmax(pmin(as.vector(final_pred_probs_on_full_data), 1 - 1e-8), 1e-8)
    #final_dev_resid_values = -2 * (y * log(final_pred_probs_on_full_data) + (1 - y) * log(1 - final_pred_probs_on_full_data))
    final_dev_resid_values = deviance_res(y,final_pred_probs_on_full_data)

    #scaled_score_x_final = scale_norm(score_x_update)
    scaled_score_x_final = min_max_scale(initial_subset$score_x)
    scaled_dev_resid_final = min_max_scale(abs(final_dev_resid_values))

    final_combined_score = tuned_alpha_wt * scaled_score_x_final + (1 - tuned_alpha_wt) * scaled_dev_resid_final
    final_x_score_component = tuned_alpha_wt * scaled_score_x_final
    final_y_score_component = (1 - tuned_alpha_wt) * scaled_dev_resid_final

    final_coefs_sparse = coef(overall_best_model, s = best_lambda_for_csteps)
    final_coefs_beta_vector = as.vector(overall_best_model$beta) # Beta coefficients only
    names(final_coefs_beta_vector) = rownames(overall_best_model$beta)
    print(paste("wt: ",tuned_alpha_wt))
    print(paste("etune: ",best_alpha_glmnet))
    #print(final_combined_score);print(final_x_score_component);print(final_y_score_component)

    #add diagnostic data to the return list ###
    convergence_history <- do.call(rbind, convergence_history_list)
    rownames(convergence_history) <- NULL

    return(list(model = overall_best_model,
                best_alpha_glmnet = best_alpha_glmnet,
                best_lambda = best_lambda_for_csteps,
                tuned_alpha_wt = tuned_alpha_wt,
                subset_used_for_final_model = overall_best_subset_c,
                coefs_with_intercept = final_coefs_sparse,
                coefs_beta_only = final_coefs_beta_vector,
                prob_on_full_data = final_pred_probs_on_full_data,
                combined_score_final = final_combined_score,
                total_deviance_on_full_data = min_total_deviance_on_full_data,
                dev_resid_values_final = final_dev_resid_values,
                x_score_component_final = final_x_score_component,
                y_score_component_final = final_y_score_component,
                person_res = final_dev_resids_before_reweight,
                reweight = observation_weights_reweight,
                mad_cutoff = cutoff_resids_dynamic,
                raw_beta = raw_beta,
                raw_model = raw_model,
                convergence_history = convergence_history,
                cv_grid_initial = cv_grid_initial,
                scaling_params_list = scaling_params_list

    ))
}
