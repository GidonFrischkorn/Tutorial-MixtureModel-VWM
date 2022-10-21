#' @title Means +/- SE
#'
#' @description This function computes the mean of the variable x plus the deviations
#' of one standard error below and above the mean. This is mainly used for plotting.
#' %% this is how to comment in roxygen2 docs, by the way.
#'
#' @usage mean_se(x)
#' %% Especially useful if you there are non-important pre-defined arguments to your function like non_sense in this example
#'
#' @param x a vector of data points for which the mean +/- se should be computed
#' @param mult the multiplier to use for computing the standard error
#'
#' @return a data frame containing three values: mean_x, mean_x - se, and mean_x + se
#'
#' @details 
#'
#' @examples
#'
#' x = rnorm(100, mean = 10, sd =5)
#' df_meanSE <- mean_se(x)
#' 
#' @author Ven Popov
#'
#' @export 
#' 
mean_se2 <- function (x, mult = 1.96) 
{
  x <- stats::na.omit(x)
  se <- mult * sqrt(stats::var(x)/length(x))
  mean <- mean(x)
  data.frame(y = mean, ymin = mean - se, ymax = mean + se)
}