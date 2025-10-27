#' Robustly scale a matrix or vector
#'
#' This function performs robust scaling on a numeric matrix or vector by
#' centering the data using the median and scaling by the median absolute deviation (MAD).
#' This method is more resistant to outliers than traditional scaling with the
#' mean and standard deviation.
#'
#' @param X_mat A numeric matrix or object coercible to a numeric matrix.
#'
#' @return A list containing three elements:
#' \itemize{
#'   \item \strong{X_scaled:} A matrix with the scaled data.
#'   \item \strong{Median:} A numeric vector of the median for each column, used for centering.
#'   \item \strong{MAD:} A numeric vector of the MAD for each column, used for scaling.
#' }
#' If a column is constant or has all `NA`s after removing `NA`s, the scale is treated as zero, and the corresponding scaled column will be all zeros.
#'
#' @details The function iterates through each column of the input matrix.
#' For each column, it computes the median and the MAD.
#' Each data point \eqn{x_i} is then transformed using the formula \deqn{(x_i - \text{median}) / \text{MAD}},
#' The constant 1.4826 is used to make the MAD comparable to the standard deviation under the assumption of a normal distribution.
#' If a column's MAD is effectively zero, the scaled values for that column are set to zero to prevent division-by-zero errors.
#'
#' @examples
#' # Example with a standard matrix
#' mat <- matrix(c(1, 2, 3, 10, 20, 30, 100, 200, 300), nrow = 3)
#' robust_scale(mat)
#'
#' # Example with an outlier
#' mat_outlier <- matrix(c(1, 2, 3, 10, 20, 1000, 100, 200, 300), nrow = 3)
#' robust_scale(mat_outlier)
#'
#' # Example with a single column that is constant
#' const_col <- matrix(c(5, 5, 5, 10, 20, 30), ncol = 2)
#' robust_scale(const_col)
#'
#' # Example with a data frame
#' df <- data.frame(a = c(1, 2, 3), b = c(10, 20, 30))
#' robust_scale(df)
#'
#' @seealso \code{\link[stats]{median}}, \code{\link[stats]{mad}}
#'
#' @export
robust_scale = function(X_mat) {
    # Function code remains the same
    if (!is.matrix(X_mat)) {
        X_mat_input <- as.matrix(X_mat)
        if (!is.numeric(X_mat_input)) {
            stop("Input X_mat cannot be coerced to a numeric matrix.")
        }
    } else {
        X_mat_input <- X_mat
    }

    num_cols <- ncol(X_mat_input)
    num_rows <- nrow(X_mat_input)

    # Initialize objects to store results
    X_scaled_mat <- matrix(NA_real_, nrow = num_rows, ncol = num_cols)
    centers <- numeric(num_cols)
    scales <- numeric(num_cols)

    # Preserve column names if they exist
    if (!is.null(colnames(X_mat_input))) {
        colnames(X_scaled_mat) <- colnames(X_mat_input)
        names(centers) <- colnames(X_mat_input)
        names(scales) <- colnames(X_mat_input)
    }


    if (num_cols == 0 || num_rows == 0) { # Handle empty matrix case
        warning("Input matrix is empty.")
        return(list(X_scaled = X_scaled_mat, center = centers, scale = scales))
    }

    for (j in 1:num_cols) {
        col_data <- X_mat_input[, j]

        # Calculate median (center)
        m <- median(col_data, na.rm = TRUE)
        centers[j] <- m

        # Calculate MAD (scale)
        # Using constant = 1.4826 makes MAD comparable to SD for normal data
        mad_val <- mad(col_data, center = m, constant = 1.4826, na.rm = TRUE)

        if (is.na(mad_val) || mad_val < .Machine$double.eps) {
            # If MAD is zero or NA (e.g., column is constant or all NAs after median)
            # Set scale to 1 to avoid division by zero; effectively means no scaling, only centering.
            # Scaled values become 0 if they were equal to the median, or their deviation from median if scale is 1.
            # Your original code returned rep(0, length(col)). This implies (col - m) / mad_val results in 0s
            # if mad_val is effectively infinite or if we force result to 0.
            scales[j] <- ifelse(is.na(mad_val), NA_real_, mad_val) # Store the problematic MAD (or NA)
            X_scaled_mat[, j] <- 0
        } else {
            scales[j] <- mad_val
            X_scaled_mat[, j] <- (col_data - m) / mad_val
        }
    }

    return(list(X_scaled = X_scaled_mat, Median = centers, MAD = scales))
}
