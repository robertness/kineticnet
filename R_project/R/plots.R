#' Posterior histogram
#'
#' Generate a plot of the posterior of a parameter from a Stan fit.
#' @param param name of a parameter
#' @param fit Stanfit object
#' @param ground_truth value of the 'true' value, defaults to NULL

post_hist <- function(param,  fit, ground_truth = NULL){
  vals <- rstan::extract(fit, param)[[1]]
  rng <- range(vals)
  if(!is.null(ground_truth)){
    rng <- range(c(ground_truth, vals))
  }
  hist(vals, xlim = rng, main = paste(param, 'posterior'), col = 'grey')
  if(!is.null(ground_truth)){
    abline(v = params[[param]], col = 'red')
  }
}

#' Posterior prediction histogram
#'
#' Generate a histogram of the posterior predictive distribution and overlay it on the emperical distribution.
#' @param var_name
#' @param fit Stanfit object
#' @param obs observed values of the target variable
post_pred <- function(var_name, fit, obs){
  var_name_sim <- paste(var_name, "sim", sep = "_")
  sim = colMeans(rstan::extract(fit, var_name_sim)[[1]])
  N <- length(obs) # should be the same as the length of sim.
  pd <- data.frame(
    val = c(obs, sim),
    distribution = c(rep('empirical', N), rep('predicted', N))
  )
  ggplot(pd, aes(x = val, fill = distribution)) +
    geom_histogram(alpha = 0.5) +
    ggtitle(paste0(var_name, ": empirical vs post. predicted"))
}
