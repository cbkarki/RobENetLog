#' Plot C-step convergence history
#'
#' Generates two ggplot objects showing the convergence history of
#' C-steps for different alpha weights, and the convergence path for the
#' optimal alpha weight. The function uses `ggplot2` for plotting, `viridis`
#' for color scales, and `patchwork` for arranging the two plots.
#'
#' @param convergence_history A data frame containing the convergence history.
#'   It is expected to have, at a minimum, the following columns:
#'   \itemize{
#'     \item `step`: The C-step number.
#'     \item `subset_deviance`: The deviance on the current subset.
#'     \item `alpha_wt`: The alpha weight used for that step.
#'   }
#'
#' @return The function prints a composite plot to the active graphics device
#'   but returns `NULL` invisibly. The composite plot is made of two subplots
#'   arranged vertically:
#'   \enumerate{
#'     \item **C-Step Convergence Plot:** Shows all C-step convergence lines,
#'           with each line representing a different alpha weight.
#'     \item **Optimum C-Step Convergence Plot:** Shows only the C-step
#'           convergence line for the alpha weight that resulted in the minimum
#'           subset deviance.
#'   }
#'
#' @import ggplot2
#' @import viridis
#' @import tidyr
#' @import patchwork
#' @noRd
#'
#' @examples
#' # Generate dummy data for demonstration
#' convergence_history_dummy <- expand.grid(step = 1:10, alpha_wt = seq(0.1, 1.0, by = 0.1))
#' convergence_history_dummy$subset_deviance <- with(
#'   convergence_history_dummy,
#'   (alpha_wt - 0.5)^2 * 100 + step * 5 + rnorm(nrow(convergence_history_dummy), sd = 2)
#' )
#'
#' # Call the function with the dummy data
#' plot_c_step_conv(convergence_history_dummy)

plot_c_step_conv = function(convergence_history) {
    # plot convergence: c-steps while tunning alpha wt
    min_alpha_wt = convergence_history[which.min(convergence_history$subset_deviance),]$alpha_wt

    # plot 1: convergence C-step
    p1 =
        ggplot(convergence_history, aes(x = step, y = subset_deviance, group = alpha_wt, color = alpha_wt)) +
        geom_line() +
        geom_point(size = 1) +
        scale_color_viridis_c(name = "alpha_wt") +
        labs(
            title = paste("C-Step Convergence Plot"),
            subtitle = "Each line represents a different alpha weight",
            x = "C-Step Number",
            y = "Deviance on Current Subset"
        ) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))


    # plot 2: the optimized path of c-step convergence
    p2 =
        ggplot(convergence_history %>% filter(alpha_wt == min_alpha_wt  ), aes(x = step, y = subset_deviance, group = alpha_wt, color = alpha_wt)) +
        geom_line() +
        geom_point(size = 1) +
        scale_color_viridis_c(name = "alpha_wt") +
        labs(
            title = paste("Optimum C-Step Convergence Plot"),
            subtitle = paste("The line represents optimum alpha_wt"),
            x = "C-Step Number",
            y = "Deviance on Current Subset"
        ) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

    # arrange plots
    print(p1 / p2)
}
