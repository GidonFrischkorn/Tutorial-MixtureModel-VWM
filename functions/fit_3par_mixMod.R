#' @title Fit the 3-parameter mixture model to data
#'
#' @description Using the likelihood function of the 3-parameter mixture model,
#' this functions fits and returns parameter estimates of the 3-parameter mixture model
#' of 6-distributions: 1 for correct, 4 for non-target items, and 1 guessing.
#'
#' @usage fit_3par_mixMod(dat)
#' %% Especially useful if you there are non-important pre-defined arguments to your function like non_sense in this example
#'
#' @param dat a data frame containing the values for the target item (y) and 5 non-target items
#' @param init_values initial values to start estimation. 
#'
#' @return The function returns the estimated parameters as well as the negative loglikelihood of the converged estimation
#'
#' @details 
#'
#' @examples
#' 
#' @author Ven Popov
#'
#' @export 
#' 
fit_3par_mixMod <- function(dat, init_values=list(p_correct=3, p_other=1, sigma=9)) {
  require(stats4)
  LL_resp <- LL_3par_mixMod(dat) 
  # debug(LL_resp)
  fit <- mle(LL_resp, start = init_values)
  coef <- data.frame(t(fit@coef))
  p_correct <- exp(coef$p_correct)/(exp(coef$p_correct)+exp(coef$p_other)+exp(0))
  p_other <- exp(coef$p_other)/(exp(coef$p_correct)+exp(coef$p_other)+exp(0))
  coef$p_correct <- p_correct
  coef$p_other <- p_other
  coef$p_guess <- 1-coef$p_correct-coef$p_other
  coef$negll <- summary(fit)@m2logL
  return(round(coef,3))
}