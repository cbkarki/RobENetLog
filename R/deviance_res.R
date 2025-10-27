#' Calculate deviance residuals for binomial models
#'
#' This function computes the deviance residuals, a measure of lack of fit for
#' each observation in a binomial (e.g., logistic regression) model.
#'
#' @param y A numeric vector of true binary outcomes (0 or 1).
#' @param pred_probs A numeric vector of predicted probabilities (values
#'   between 0 and 1) from the model.
#'
#' @return A numeric vector of deviance residuals, with one value for each
#'   observation.
#'
#' @details
#' The deviance residual for each observation is calculated using the formula:
#' \deqn{
#'   r_{i} = -2 * (y_i * \log(\hat{p}_i) + (1 - y_i) * \log(1 - \hat{p}_i))
#' }
#' where \eqn{y_i} is the true outcome and \eqn{\hat{p}_i} is the predicted
#' probability. A small value for the residual indicates a good fit for that
#' observation. The function takes the square root of the expression to align
#' with the standard definition of a deviance residual.
#'
#' This function includes a small safeguard to handle cases where predicted
#' probabilities are exactly 0 or 1, preventing errors with `log(0)`.
#'
#' @examples
#' # Example with some true outcomes and predicted probabilities
#' true_y <- c(1, 0, 1, 1, 0)
#' pred_p <- c(0.9, 0.1, 0.6, 0.8, 0.4)
#' deviance_res(true_y, pred_p)
#'
#' # Example with perfect prediction for one observation
#' perfect_pred <- c(1, 0, 1, 1, 0)
#' pred_p_perfect <- c(1, 0.1, 0.6, 0.8, 0.4)
#' deviance_res(true_y, pred_p_perfect)
#'

deviance_res <- function(y, pred_probs) {
    if (length(y) != length(pred_probs)) {
        stop("y and pred_probs must have the same length.")
    }

    # Safeguard against log(0) and log(1) for numerical stability
    pred_probs <- pmax(.Machine$double.eps, pmin(1 - .Machine$double.eps, pred_probs))

    # Calculate the deviance residual
    residuals <- -2 * (y * log(pred_probs) + (1 - y) * log(1 - pred_probs))
    return(residuals)
}
