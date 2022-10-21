#' @title Likelihood function of the 3-parameter mixture model (Bays et al., 2009)
#'
#' @description This function provides the likelihood of the 3-parameter mixture model
#' first introduced by Bays et al. (2009). It is used to estimate parameter of 
#' the model using maximum likelihood estimation. The function return the likelihood
#' of a set of parameters given the specified data (dat)
#' %% this is how to comment in roxygen2 docs, by the way.
#'
#' @usage LL_3par_mixMod(dat)
#' %% Especially useful if you there are non-important pre-defined arguments to your function like non_sense in this example
#'
#' @param dat a data frame containing the values for the target item (y) and 5 non-target items 
#'
#' @return the negative summed loglikelihood for the data to belong to the target distribution, any of the non-target distributions, or the guessing distribution
#'
#' @details 
#'
#' @examples
#' 
#' @author Ven Popov
#'
#' @export 
#' 
LL_3par_mixMod <- function(dat) {
  # transform x from degrees to radians
  rad = dat[c('y','V1','V2','V4','V5')]
  function(p_correct=2, p_other=1, sigma=9) {
    # trasnform laten prob
    p_c = exp(p_correct)/(exp(p_correct)+exp(p_other)+exp(0))
    p_o = exp(p_other)/(exp(p_correct)+exp(p_other)+exp(0))
    p_g = exp(0)/(exp(p_correct)+exp(p_other)+exp(0))
    # transform the normal sd into radians kappa for circular vonmises concentration parameter
    rad_sigma = sigma * pi /180
    kappa = (1/rad_sigma) ** 2
    l_norm <- brms::dvon_mises(rad$y, mu=0, kappa=kappa)
    l_norm1 <- brms::dvon_mises(rad$y, mu=rad$V1, kappa=kappa)
    l_norm2 <- brms::dvon_mises(rad$y, mu=rad$V2, kappa=kappa)
    l_norm4 <- brms::dvon_mises(rad$y, mu=rad$V4, kappa=kappa)
    l_norm5 <- brms::dvon_mises(rad$y, mu=rad$V5, kappa=kappa)
    l_unif <- brms::dvon_mises(rad$y, mu=0, kappa=0)
    likelihood <- p_c*l_norm + p_o/4*(l_norm1+l_norm2+l_norm4+l_norm5) + p_g*l_unif
    -sum(log(likelihood))
  }
}