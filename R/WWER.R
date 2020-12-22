# TODO(brielin): Add a function combining welch_test and wwer for a single
# pair of phenotypes (instead of just all pairs below).

#' Simple implementation of a welch test.
#'
#' This tests the null hypothesis abs(beta1) = abs(beta2) against the two
#' alternatives abs(beta1) > abs(beta2) and abs(beta1) < abs(beta2).
#'
#' @param beta1 Float or vector of floats, mean of the first sample.
#' @param beta2 Float or vector of floats, mean of the second sample.
#' @param se1 Float or vector of floats, SD of estimate of mu1.
#' @param se2 Float or vector of floats, SD of estimate of mu2.
#' @param welch_thresh Float, p_value threshold for significance.
welch_test <- function(beta1, se1, beta2, se2, welch_thresh = 0.05) {
  n1 <- 1 / (se1**2)
  n2 <- 1 / (se2**2)
  t_val <- (abs(beta2) - abs(beta1)) / sqrt(se1^2 + se2^2)
  nu <- (se1^2 + se2^2)^2 / (se1^4 / (n1 - 1) + se2^4 / (n2 - 1))
  res_12 <- stats::pt(t_val, round(nu))
  res_21 <- 1-res_12
  sig_12 <- res_12 < welch_thresh
  sig_21 <- res_21 < welch_thresh
  res <- sig_12 - sig_21
  res[is.na(t_val)] <- 0
  return(list("pass" = res, "t" = t_val, "nu" = nu) )
}

#' Fits a regression of outcome on exposure using SE and weights as weights.
#' 
#' @param b_exp Vector of floats. SNP effects on exposure.
#' @param b_out Vector of floats. SNP effects on outcome.
#' @param se_exp Vector of floats. Standard error of SNP effects on exposure.
#' @param se_out Vector of floats. Standard error of SNP effects on outcome.
#' @param weights Vector of floats. Welch weights to use for regression, will be
#'   combined with the SE of the outcome to determine final weights.
#' @export
WWER <- function(b_exp, b_out, se_exp, se_out, weights){
  lm_res <- summary(stats::lm(
    sign(b_exp)*b_out ~ abs(b_exp),
    weights = weights/(mean(weights)*(se_out**2))))
  beta_hat <- lm_res$coefficients[2, 1]
  beta_se <- lm_res$coefficients[2, 2]/min(lm_res$sigma, 1)
  beta_p <- 2 * (1 - stats::pnorm(abs(beta_hat / beta_se)))
  return(
    list("beta.hat" = beta_hat, "beta.se" = beta_se, "beta.p.value" = beta_p))
}