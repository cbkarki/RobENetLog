#' Plot initial hyperparameter tuning results
#'
#' This function generates a composite plot displaying the results of an
#' initial grid search for hyperparameters (alpha and lambda). The first plot
#' is a heatmap showing cross-validation deviance (`cvm`) across the grid,
#' and the second plot shows the minimum deviance path for the optimal alpha
#' value.
#'
#' @param cv_grid_initial A data frame containing the results of a grid search
#'   for hyperparameter tuning. It must contain the following columns:
#'   \itemize{
#'     \item `alpha`: The alpha value used in the model.
#'     \item `lambda`: The lambda value used in the model.
#'     \item `cvm`: The cross-validation mean deviance for the given alpha and lambda.
#'   }
#'
#' @return The function prints a composite plot to the active graphics device
#'   but returns `NULL` invisibly. The composite plot is made of two subplots
#'   arranged vertically:
#'   \enumerate{
#'     \item **Initial Hyperparameter Tuning:** A heatmap displaying the `cvm`
#'           values across the tested alpha and log(lambda) grid. The point
#'           with the minimum `cvm` is highlighted.
#'     \item **Optimum Initial Hyperparameter Tuning:** A line plot showing
#'           the `cvm` values as a function of log(lambda) for the optimal
#'           alpha value found in the grid search.
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
#' cv_grid_initial_dummy <- expand.grid(
#'   alpha = c(0.1, 0.5, 1.0),
#'   lambda = 10^seq(-3, 0, by = 0.5)
#' )
#' cv_grid_initial_dummy$cvm <- with(
#'   cv_grid_initial_dummy,
#'   (alpha - 0.5)^2 * 100 + log(lambda) * 5 + rnorm(nrow(cv_grid_initial_dummy))
#' )
#'
#' # Call the function with the dummy data
#' plot_al_la_heat(cv_grid_initial_dummy)
plot_al_la_heat = function(cv_grid_initial) {

    # plot 3:  heat map of grid search alpha nd lambda for initial subset
    min_point <- cv_grid_initial[which.min(cv_grid_initial$cvm), ]

    # curve for min deviance pair of alpha and lambda for which has the minimum cvm

    dev_min = cv_grid_initial %>% filter(alpha == min_point$alpha)

    p3 =
        ggplot(cv_grid_initial, aes(x = log(lambda), y = alpha, color = cvm)) +
        # Use geom_point to plot each calculated value
        geom_point(size = 3, alpha = 0.9) +
        # Use a color scale where lower is better
        scale_color_viridis_c(name = "CV Deviance", direction = -1) +
        # Highlight the minimum point
        geom_point(data = min_point, aes(x = log(lambda), y = alpha),
                   color = "red", size = 5, shape = 1, stroke = 1.5) + # Circle around the best point
        labs(
            title = paste("Initial Hyperparameter Tuning"),
            subtitle = paste0("Minimum CV Deviance at alpha=", round(min_point$alpha, 2),
                              ", log(lambda)=", round(log(min_point$lambda), 2)),
            x = "log(Lambda)",
            y = "Alpha"
        ) +
        scale_x_reverse() +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5)
        )

    # plot 4: minimum deviance for set from the grid search of alpha and lambda
    p4 =
        ggplot(dev_min, aes(x = log(lambda), y = cvm)) +
        # Use geom_point to plot each calculated value
        geom_point(size = 2, alpha = 0.8)  +
        # Use a color scale where lower is better
        scale_color_viridis_c(name = "CV Deviance", direction = -1) +
        # Highlight the minimum point
        # geom_point(data = min_point, aes(x = log(lambda), y = cvm),
        #            color = "red", size = 5, shape = 1, stroke = 1.5) + # Circle around the best point
        geom_vline(xintercept = log(min_point$lambda),linetype = 2,color="red")+
        labs(
            title = paste("Optimum Initial Hyperparameter Tuning"),
            subtitle = paste0("Minimum CV Deviance=", round(min(dev_min$cvm),2), "at alpha=", round(min(min_point$alpha), 2),
                              ", log(lambda)=", round(min(log(min_point$lambda)), 2)),
            x = "log(lambda)",
            y = "CV Devaince"
        ) +
        scale_x_reverse() +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5)
        )
    print(p3 / p4)
}
