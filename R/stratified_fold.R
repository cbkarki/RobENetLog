#' Create stratified k-folds for cross-validation
#'
#' This function generates a vector of fold IDs for stratified k-fold
#' cross-validation, ensuring that each fold has approximately the same
#' proportion of the response variable's classes as the full dataset.
#' It is a wrapper for `caret::createFolds()` with added robust handling of
#' edge cases.
#'
#' @param y_vec A numeric or character vector of the response variable. It will
#'   be coerced to a factor internally.
#' @param k_folds An integer specifying the number of folds for cross-validation.
#'
#' @return An integer vector of fold IDs, with the same length as `y_vec`.
#'   Each observation is assigned a number from 1 to `k_folds`.
#'
#' @details
#' The function includes several checks to handle potential issues:
#' \itemize{
#'   \item If the number of observations is less than `k_folds`, it reduces `k_folds`
#'   to the number of observations and issues a warning.
#'   \item If only one class is present in the response variable, it falls back
#'   to simple random sampling and issues a warning.
#'   \item If the smallest class count is less than `k_folds`, it reduces `k_folds`
#'   to the smallest class count (or a minimum of 2) and issues a warning.
#'   \item If the call to `caret::createFolds()` fails for any reason, it also
#'   falls back to simple random sampling.
#' }
#'
#' @examples
#' # Example 1: Standard stratified folding
#' set.seed(123)
#' y <- c(rep(1, 40), rep(2, 60))
#' folds <- stratified_fold(y, k_folds = 10)
#' folds
#' table(y, folds)
#'
#' # Example 2: Edge case with small number of observations
#' y_small <- c(1, 2, 1, 2)
#' folds_small <- stratified_fold(y_small, k_folds = 5)
#' folds_small
#' table(y_small, folds_small)
#'
#' # Example 3: Edge case with only one class
#' y_single_class <- rep(1, 20)
#' folds_single <- stratified_fold(y_single_class, k_folds = 5)
#' folds_single
#' table(y_single_class, folds_single)
#'
#' @importFrom caret createFolds
#' @export
stratified_fold <- function(y_vec, k_folds) {
    n_obs <- length(y_vec)
    if (n_obs < k_folds) {
        warning("Number of observations is less than k_folds. Setting k_folds to n_obs.")
        k_folds <- n_obs
    }
    y_factor <- as.factor(y_vec)
    if (length(levels(y_factor)) < 2) {
        warning("Only one class present for stratified_fold. Using simple random folds.")
        return(sample(rep(1:k_folds, length.out = n_obs)))
    }
    min_class_count <- min(table(y_factor))
    if(min_class_count < k_folds && min_class_count > 0){
        k_folds <- max(2, min_class_count) # Ensure k_folds is at least 2
        warning(paste("Smallest class count is less than k_folds. Reducing k_folds to", k_folds))
    } else if (min_class_count == 0) {
        warning("One class has zero observations. Using simple random folds.")
        return(sample(rep(1:k_folds, length.out = n_obs)))
    }

    foldid <- integer(n_obs)
    folds_list <- tryCatch(
        caret::createFolds(y_factor, k = k_folds, list = TRUE, returnTrain = FALSE),
        error = function(e) {
            warning("caret::createFolds failed: ", e$message, ". Using simple random folds.")
            NULL
        }
    )
    if (is.null(folds_list)) {
        return(sample(rep(1:k_folds, length.out = n_obs)))
    }
    for (i in 1:length(folds_list)) {
        foldid[folds_list[[i]]] <- i
    }
    return(foldid)
}
