#' Scales a numeric vector to a [0, 1] range.
#'
#' @description
#' This function takes a vector and performs min-max normalization
#'
#' @param x numeric vector.
#'
#' @return Returns the scaled vector, all values in the 0-1 range. If all the values in the vector are same it will return 0.5 which is the mid value
#'
#' @examples
#' # example 1
#' vec = sample(1:100,20)
#' min_max_scale(vec)
#'
#' # example 2
#' vec = rep(5, 10)
#' min_max_scale(vec)
#'
#' @export
min_max_scale = function(x) {
    if (length(unique(x)) == 1) {
        return(rep(0.5, length(x))) # Neutral value if constant
    }
    (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE) + .Machine$double.eps)
}


