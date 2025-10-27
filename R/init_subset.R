#' @title Create an initial subset based on robust outlyingness
#'
#' @description
#' This internal function identifies an initial subset of observations that
#' are considered "non-outliers" based on a robust measure of outlyingness.
#' This subset can be used for robust hyperparameter tuning in models like
#' Elastic Net.
#'
#' @param X A numeric matrix or object coercible to a numeric matrix.
#' @param standardize Logical. If `TRUE` (the default), the input matrix `X`
#'   is first scaled using a robust method (median and MAD).
#' @param initial_subset_cutoff_factor A numeric value used to determine the
#'   outlyingness cutoff. The cutoff is defined as,
#'  \deqn{\text{median}(score_x) + \text{initial\_subset\_cutoff\_factor} \times \text{mad}(score_x)}
#'  where, \deqn{score_{x_{i}} = \sum_{j=1}^p |X_{scaled_{ij}}|}
#'  \deqn{ X_{scaled_{ij}} = \frac{X_{ij}- \text{Med}(X_{j})}{\text{MAD}(X_j)} }
#'  \deqn{\text{MAD} = \text{median}(|x_i - \text{median}(x)|)}
#'
#' @return A list with five elements:
#' \itemize{
#'   \item \strong{subset_init:} An integer vector of row indices corresponding
#'     to the observations in the initial non-outlier subset.
#'   \item \strong{cutoff_score_x_dynamic:} The dynamic cutoff value used to
#'     determine the initial subset.
#'   \item \strong{score_x:} A numeric vector of outlyingness scores for each observation.
#'   \item \strong{Median:} A numeric vector containing the median of each column of `X`,
#'     used for centering.
#'   \item \strong{MAD:} A numeric vector containing the median absolute
#'     deviation (MAD) of each column of `X`, used for scaling.
#'   \item \strong{X_scaled:} Median centered and MAD scaled of the input matrix if `standardize = TRUE`, if `FALSE` same input matrix is returned.
#'  \item \strong{subset_X:} Initial robust subset of `X_scaled` based on `cutoff_score_x_dynamic`. `subset_X` will be the same as `X_scaled`
#'      if no `x_score` exceed the `cutoff_score_x_dynamic`.
#'
#' }
#'
#' @details
#' The function calculates an outlyingness score for each observation based on
#' the scaled absolute deviations from the median across all columns.
#' The method relies on the `robust_scale()` function, which must be available
#' internally within the package.
#'
#' @importFrom stats median mad rnorm
#' @examples
#' # Assuming robust_scale is an internal function available
#' # within the package.
#' # Generate some sample data
#' set.seed(123)
#' X_data <- matrix(rnorm(100), ncol = 5)
#'
#' # Run the function with default parameters
#' subset_result <- init_subset(X_data)
#'
#' # Inspect the results
#' str(subset_result)
#'
#' # Example with standardization disabled
#' subset_result_no_scale <- init_subset(X_data, standardize = FALSE)
#' str(subset_result_no_scale)
#'
#' # Example with a different cutoff factor
#' subset_result_strict <- init_subset(X_data, initial_subset_cutoff_factor = 2)
#' str(subset_result_strict)
#'
#' @export
init_subset = function(X,
                       standardize=TRUE,
                       initial_subset_cutoff_factor=3
) {
    n = nrow(X)
    p_cols = ncol(X)

    # matrix conversion
    if(!is.matrix(X)) {
        X = as.matrix(X)
    }

    # robust scaling
    scale_result = robust_scale(X)

    # if standardize is true, the robust scaling and outlyingness, if false no scaling is done
    if (standardize) {
        X_scaled = scale_result$X_scaled
        scaling_params_list <- list(centers = scale_result$Median, scales = scale_result$MAD)
    } else {
        X_scaled = as.matrix(X)
        scaling_params_list <- list(centers = rep(0, p_cols), scales = rep(1, p_cols))
    }
    if(!is.matrix(X_scaled)) X_scaled = as.matrix(X_scaled)

    # absolute values or each value, hence obtaining a positive score for each observation
    mad_outlyingness = apply(scale_result$X_scaled, 2, function(col) abs(col))

    # rowsum will provide score for each observation
    score_x = rowSums(mad_outlyingness, na.rm = TRUE)

    # percentile cutoff based on median based outlyingness
    # median of x_scores
    median_score_x = median(score_x, na.rm = TRUE)

    # mad of x_scores
    mad_score_x = mad(score_x, center = median_score_x, constant = 1.4826, na.rm = TRUE)

    # defining cutoff
    cutoff_score_x_dynamic = median_score_x + initial_subset_cutoff_factor * mad_score_x

    # initial subset, which will be used for hyperparameter tunning like alpha and lambda for Elastic Net
    subset_init = which(score_x <= cutoff_score_x_dynamic)

    return(list(subset_init = subset_init,
                cutoff_score_x_dynamic = cutoff_score_x_dynamic,
                score_x = score_x,
                Median = scaling_params_list$centers,
                MAD = scaling_params_list$scales,
                X_scaled = X_scaled,
                subset_X = X_scaled[subset_init,]
                )
           )
}
