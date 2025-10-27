#' @title Perform initial hyperparameter tuning on a subset
#'
#' @description
#' This internal function performs cross-validated hyperparameter tuning for
#' an Elastic Net logistic regression model. The tuning is done on an initial
#' subset of non-outlying observations to ensure robustness. The function
#' identifies the best `alpha` and `lambda` for the given subset.
#'
#' @param alpha_seq A numeric vector of `alpha` values to be tested during the
#'   cross-validation.
#' @param y_subset_init The numeric vector of response values (0 or 1) corresponding
#'   to the initial non-outlier subset.
#' @param nfolds_init_cv An integer specifying the number of folds for cross-validation
#'   on the initial subset.
#' @param subset_init An integer vector of row indices of the initial non-outlier subset.
#' @param X_scaled A numeric matrix of the predictor variables, already scaled.
#'   This is typically the `X_scaled` output from `init_subset`.
#' @param subset_X A numeric matrix of predictor variables corresponding to the
#'   initial non-outlier subset.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item \strong{best_alpha_glmnet:} The `alpha` value that resulted in the
#'     best (minimum) mean cross-validated deviance.
#'   \item \strong{best_lambda_for_csteps:} The `lambda.1se` from the best `glmnet`
#'     model, or `lambda.min` if `lambda.1se` is not available.
#'   \item \strong{cv_grid_initial:} A data frame containing the full cross-validation
#'     grid of deviance values for all tested `alpha` and `lambda` combinations.
#' }
#'
#' @details
#' The function uses `glmnet::cv.glmnet()` to perform a grid search over the
#' specified `alpha_seq`. Stratified k-fold cross-validation is used to ensure
#' a proportional representation of classes in each fold. The function is designed
#' to handle cases where `glmnet` fails for certain `alpha` values by continuing
#' the search. If no valid model is found, an error is thrown. The resulting
#' best `alpha` and corresponding `lambda` (or `lambda.min`) are returned for
#' subsequent steps.
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom magrittr %>%
#' @importFrom stats median mad
#' @noRd
#'
.pre_tuning = function(alpha_seq = alpha_seq,
                       y_subset_init = y_subset_init,
                       nfolds_init_cv = nfolds_init_cv,
                       subset_init = initial_subset$subset_init,
                       X_scaled = initial_subset$X_scaled,
                       subset_X = initial_subset$subset_X) {

    best_alpha_glmnet = alpha_seq[1]
    best_cvm_init = Inf
    best_lambda_for_csteps = NULL # Initialize

    # strarification for proportionate representation of classes in the fold as of response variable
    fold_id_init = stratified_fold(y_subset_init, nfolds_init_cv)

    # dataframe to store full cv grid
    cv_grid_initial_list <- list()

    for (a_glmnet in alpha_seq) {
        cv_results = tryCatch(
            glmnet::cv.glmnet(subset_X, y_subset_init, alpha = a_glmnet,
                              family = "binomial", nfolds = nfolds_init_cv, type.measure = "deviance", foldid = fold_id_init,
                              standardize = FALSE), # X is already scaled
            error = function(e) {
                warning(sprintf("cv.glmnet failed for alpha_glmnet=%.2f: %s", a_glmnet, e$message))
                return(NULL)
            })

        ### store the results for the heatmap ###
        if (!is.null(cv_results)) {
            cv_grid_initial_list[[as.character(a_glmnet)]] <- data.frame(
                alpha = a_glmnet,
                lambda = cv_results$lambda,
                cvm = cv_results$cvm # Mean CV error (deviance)
            )
        }

        if (!is.null(cv_results) && !is.null(cv_results$cvm) && min(cv_results$cvm, na.rm=TRUE) < best_cvm_init) {
            best_cvm_init = min(cv_results$cvm, na.rm=TRUE)
            best_alpha_glmnet = a_glmnet
            best_lambda_for_csteps = cv_results$lambda.1se
            if(is.null(best_lambda_for_csteps) || !is.finite(best_lambda_for_csteps)) {
                best_lambda_for_csteps = cv_results$lambda.min
            }
        }
    } # end of alpha and lambda tuning,a_glmnet loop


    if (best_cvm_init == Inf || is.null(best_lambda_for_csteps)) {
        stop("Failed to tune initial glmnet alpha/lambda. cv.glmnet might have failed for all alpha values or lambda selection failed.")
    }

    # combine the stored CV grid ###
    cv_grid_initial <- do.call(rbind, cv_grid_initial_list)
    rownames(cv_grid_initial) <- NULL


    return(list(best_alpha_glmnet = best_alpha_glmnet,
                best_lambda_for_csteps = best_lambda_for_csteps,
                cv_grid_initial = cv_grid_initial ))
}


