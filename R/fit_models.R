#' Selects SNPs for inclusion in MR by comparing per-variance effect sizes.
#'
#' Notes: sumstats passed to this function must be computed on the per-variance
#' scale.
#'
#' @param sumstats List with elements "beta_hat", "se_hat", both M x D matrices.
#' @param snps_to_use List or NULL. A list named by phenotypes where each list
#'   entry is a list of SNPs that can be used for that phenotype. Usually
#'   the result of clumping to avoid correlated SNPs.
#' @param p_thresh Float, p-value threshold to use for SNP inclusion.
#' @param exclusive Bool. True to only use SNPs significant for one phenotype
#'   but *not* the other.
#' @param weight Bool. True to store welch-test weights for regression.
#' @param filter Double of NULL. If not NULL, filter variants with welch
#'   statistic less than filter.
#' @param verbose Bool. If true, print phenotype label during iteration.
#' @export
select_snps <- function(sumstats, snps_to_use = NULL, p_thresh = 5e-8,
                        exclusive = FALSE, weight = TRUE, filter = 1.6, verbose = FALSE) {
  z_scores <- as.matrix(abs(sumstats$beta_hat / sumstats$se_hat))
  p_vals <- 2 * (1 - stats::pnorm(z_scores))
  sig_p_vals <- dplyr::as_tibble(p_vals < p_thresh)

  selected_snps <- list()
  phenos <- colnames(sumstats$beta_hat)
  D <- length(phenos)
  snps <- rownames(sumstats$beta_hat)

  for(index in 1:D){
    pheno1 <- phenos[index]
    if(verbose){
      print(pheno1)
    }
    mask1 <- rep(TRUE, length(snps))
    if (!is.null(snps_to_use)) {
      p1_snps <- get(pheno1, snps_to_use)
      mask1 <- (snps %in% p1_snps)
    }
    sig1 <- dplyr::pull(sig_p_vals, pheno1)
    candidate1 <- sig1 & mask1
    candidate1[is.na(candidate1)] <- FALSE
    selected_snps[[pheno1]]$names <- snps[candidate1]
    for(pheno2 in phenos[index:D]){
      mask2 <- rep(TRUE, length(snps))
      if (!is.null(snps_to_use)) {
        p2_snps <- get(pheno2, snps_to_use)
        mask2 <- (snps %in% p2_snps)
      }
      sig2 <- dplyr::pull(sig_p_vals, pheno2)
      candidate2 <- sig2 & mask2
      candidate2[is.na(candidate2)] <- FALSE
      selected_snps[[pheno2]]$names <- snps[candidate2]

      keep <- candidate1 | candidate2
      b1 <- sumstats$beta_hat[keep, pheno1]
      b2 <- sumstats$beta_hat[keep, pheno2]
      s1 <- sumstats$se_hat[keep, pheno1]
      s2 <- sumstats$se_hat[keep, pheno2]

      welch_res <- welch_test(b1, s1, b2, s2)
      weights_12 <- -welch_res$t[candidate1[keep]]
      weights_21 <- welch_res$t[candidate2[keep]]

      if(!is.null(filter)){
        weights_12 <- ifelse(weights_12 > filter, weights_12, 0.0)
        weights_21 <- ifelse(weights_21 > filter, weights_21, 0.0)
      }
      if(isFALSE(weight)){
        weights_12 <- (weights_12 != 0)
        weights_21 <- (weights_21 != 0)
      }
      if(isTRUE(exclusive)){
        weights_12 <- weights_12 * as.numeric(!sig2[candidate1])
        weights_21 <- weights_21 * as.numeric(!sig2[candidate1])
      }

      selected_snps[[pheno1]][[pheno2]] <- weights_12
      selected_snps[[pheno2]][[pheno1]] <- weights_21
    }
  }
  return(selected_snps)
}

#' Calculates matrix of total causal effects using a specified method.
#'
#' @importFrom dplyr %>%
#'
#' @param sumstats List representing summary statistics from the second
#'   dataset. Must include entries "beta_hat" and "se_hat".
#' @param selected_snps A list of lists with names of each equal to the
#'   phenotype names. The inner lists are boolean vectors of length equal to
#'   the number of SNPs and TRUE indicating to use that SNP.
#' @param mr_method String, one of c("mean", "raps"). Method to use for TCE
#'   estimate between every pair.
#' @param min_instruments Integer. Return NA if there are less than
#'   this many instruments for a pair of phenotypes.
#' @param verbose Bpplean. True to print progress.
#' @param ... Additional parameters to pass to mr_method.
#' @export
fit_tce <- function(sumstats, selected_snps, mr_method = "egger_w",
                    min_instruments = 5, verbose = FALSE, ...) {
  mr_method_func <- switch(mr_method,
                           ps = function(b_exp, b_out, se_exp, se_out, weight){
                             suppressWarnings(mr.raps::mr.raps(
                               b_exp = b_exp, b_out = b_out, se_exp = se_exp, se_out = se_out))
                           },
                           aps = function(b_exp, b_out, se_exp, se_out, weight) {
                             suppressWarnings(mr.raps::mr.raps(
                               b_exp = b_exp, b_out = b_out, se_exp = se_exp, se_out = se_out,
                               over.dispersion = TRUE))
                           },
                           raps = function(b_exp, b_out, se_exp, se_out, weight) {
                             suppressWarnings(mr.raps::mr.raps(
                               b_exp = b_exp, b_out = b_out, se_exp = se_exp, se_out = se_out,
                               over.dispersion = TRUE, loss.function = "huber"))
                           },
                           egger_p = function(b_exp, b_out, se_exp, se_out, ...) {
                             input <- MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out, byse = se_out)
                             egger_res <- MendelianRandomization::mr_egger(input, robust = FALSE, penalized = TRUE)
                             return(list("beta.hat" = egger_res$Estimate, "beta.se" = egger_res$StdError.Est, "beta.p.value" = egger_res$Pvalue.Est))
                           },
                           egger = function(b_exp, b_out, se_exp, se_out, ...) {
                             input <- MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out, byse = se_out)
                             egger_res <- MendelianRandomization::mr_egger(input, robust = FALSE, penalized = FALSE)
                             return(list("beta.hat" = egger_res$Estimate, "beta.se" = egger_res$StdError.Est, "beta.p.value" = egger_res$Pvalue.Est))
                           },
                           mbe = function(b_exp, b_out, se_exp, se_out, ...){
                             input <- MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out, byse = se_out)
                             mbe_res <- MendelianRandomization::mr_mbe(input, seed = NA, stderror = "simple")
                             return(list("beta.hat" = mbe_res$Estimate, "beta.se" = mbe_res$StdError, "beta.p.value" = mbe_res$Pvalue))
                           },
                           median = function(b_exp, b_out, se_exp, se_out, ...){
                             input <- MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out, byse = se_out)
                             median_res <- MendelianRandomization::mr_median(input)
                             return(list("beta.hat" = median_res$Estimate, "beta.se" = median_res$StdError, "beta.p.value" = median_res$Pvalue))
                           },
                           ivw = function(b_exp, b_out, se_exp, se_out, ...){
                             input <- MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out, byse = se_out)
                             ivw_res <- MendelianRandomization::mr_ivw(input)
                             return(list("beta.hat" = ivw_res$Estimate, "beta.se" = ivw_res$StdError, "beta.p.value" = ivw_res$Pvalue))
                           },
                           mr_presso = function(b_exp, b_out, se_exp, se_out, ...){
                             input <- data.frame("b_exp" = b_exp, "b_out" = b_out, "se_exp" = se_exp, "se_out" = se_out)
                             mr_presso_res <- MRPRESSO::mr_presso(data = input, BetaOutcome = "b_out", BetaExposure = "b_exp", SdOutcome = "se_out", SdExposure = "se_exp",
                                                                  OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 1000)$`Main MR results`
                             result_index = 2
                             if(is.na(mr_presso_res$`Causal Estimate`[[result_index]])){
                               result_index = 1
                             }
                             return(list("beta.hat" = mr_presso_res$`Causal Estimate`[[result_index]], "beta.se" = mr_presso_res$Sd[[result_index]], "beta.p.value" = mr_presso_res$`P-value`[[result_index]]))
                           },
                           mr_mix = function(b_exp, b_out, se_exp, se_out, ...){
                             # TODO(brielin): This seems to flip the result?? Double check this.
                             mrmix_res <- MRMix::MRMix(b_exp, b_out, se_exp, se_out)
                             return(list("beta.hat" = mrmix_res$theta, "beta.se" = mrmix_res$SE_theta, "beta.p.value" = mrmix_res$pvalue_theta))
                           },
                           # TODO(brielin): current implementation (probably) won't work on real data
                           # because the global SNP matrix is not LD pruned (just per-phenotype).
                           # The SE is also asuming the posterior is normal which is probably wrong.
                           cause = function(X, variants){
                             params <- suppressWarnings(cause::est_cause_params(X, X$snp))
                             cause_res <- suppressWarnings(cause::cause(X=X, variants = variants, param_ests = params, force = TRUE))
                             summary_cause <- summary(cause_res)
                             quants <- summary_cause$quants[[2]]
                             beta.hat <- quants[1, "gamma"]
                             beta.se <- (quants[3, "gamma"] - quants[2, "gamma"])/(2*1.96)
                             beta.p.value <- summary_cause$p
                             return(list("beta.hat" = beta.hat, "beta.se" = beta.se, "beta.p.value" = beta.p.value))
                           },
                           egger_w = WWER
  )

  # all.equal returns a STRING if they dont have the same length??
  if (!isTRUE(all.equal(names(sumstats$beta_hat), names(selected_snps)))) {
    common_phenotypes <- intersect(
      colnames(sumstats$beta_hat), names(selected_snps))
    sumstats$beta_hat <- dplyr::select(
      sumstats$beta_hat, dplyr::one_of(common_phenotypes))
    sumstats$se_hat <- dplyr::select(
      sumstats$se_hat, dplyr::one_of(common_phenotypes))
    selected_snps <- purrr::map(selected_snps[common_phenotypes], function(x) {
      x[c(common_phenotypes, "names")]
    })
  }

  run_tce_row <- function(snps_to_use, exp) {
    if (verbose) {
      print(exp)
    }
    beta_sub <- sumstats$beta_hat[snps_to_use$names, ]
    se_sub <- sumstats$se_hat[snps_to_use$names, ]
    beta_exp <- beta_sub[, exp]
    se_exp <- se_sub[, exp]
    run_tce_entry <- function(beta_out, se_out, out){
      mask_or_weight = get(out, snps_to_use)
      snp_mask <- (mask_or_weight > 0) & !is.na(mask_or_weight) & !is.na(beta_out) & !is.na(beta_exp)
      n_instruments <- sum(snp_mask)
      if ((n_instruments < min_instruments) | (exp == out)){
        list("R" = NA, "SE" = NA, "N" = n_instruments, "p" = NA)
      } else {
        tryCatch(
          {
            if(mr_method != "cause"){
              mr_res <- mr_method_func(
                b_exp = beta_exp[snp_mask],
                b_out = beta_out[snp_mask],
                se_exp = se_exp[snp_mask],
                se_out = se_out[snp_mask],
                weight = mask_or_weight[snp_mask],
                ...
              )
            } else {
              X <- tibble::as_tibble(tibble::rownames_to_column(sumstats$beta_hat[, c(exp, out)], var = "snp")) %>%
                dplyr::rename(beta_hat_1 = exp, beta_hat_2 = out)
              X$seb1 <- sumstats$se_hat[,exp]
              X$seb2 <- sumstats$se_hat[,out]
              X$A1 <- "A"
              X$A2 <- "G"
              X <- dplyr::filter(X, !is.na(beta_hat_1), !is.na(beta_hat_2))
              variants <- dplyr::filter(X, snp %in% snps_to_use$names[snp_mask])$snp
              mr_res <- mr_method_func(X, variants)
            }

            return(list("R" = mr_res$beta.hat, "SE" = mr_res$beta.se,
                        "N" = n_instruments, "p" = mr_res$beta.p.value))
          },
          error = function(cond) {
            message(c("Error when processing ", exp, " ", out))
            message(cond)
            list("R" = NA, "SE" = NA, "N" = n_instruments, "p" = NA)
          }
        )
      }
    }
    return(purrr::pmap_dfr(
      list(beta_sub, se_sub, colnames(beta_sub)), run_tce_entry, .id = "out"))
  }
  tce_res <- purrr::imap_dfr(selected_snps, run_tce_row, .id = "exp")

  R_tce <- tce_res %>%
    tidyr::pivot_wider(
      names_from = "out", values_from = "R", id_cols = "exp") %>%
    tibble::column_to_rownames("exp")
  SE_tce <- tce_res %>%
    tidyr::pivot_wider(
      names_from = "out", values_from = "SE", id_cols = "exp") %>%
    tibble::column_to_rownames("exp")
  N_obs <- tce_res %>%
    tidyr::pivot_wider(
      names_from = "out", values_from = "N", id_cols = "exp") %>%
    tibble::column_to_rownames("exp")
  p_val <- tce_res %>%
    tidyr::pivot_wider(
      names_from = "out", values_from = "p", id_cols = "exp") %>%
    tibble::column_to_rownames("exp")
  diag(R_tce) <- 1.0
  diag(SE_tce) <- 0.0
  diag(N_obs) <- 0.0
  diag(p_val) <- 1.0

  return(list("R_tce" = as.matrix(R_tce), "SE_tce" = as.matrix(SE_tce),
              "N_obs" = as.matrix(N_obs), "p_val" = as.matrix(p_val)))
}


#' Helper function to do basic filtering of the TCE matrix.
#'
#' Large values of R, entries with a high SE, and row/columns with many nans
#' can be removed.
#'
#' @param R_tce Matrix or data.frame. Estimates of TCE.
#' @param SE_tce Matrix or data.frame. Standard errors of the entries in R_tce.
#' @param max_R Float. Set all entries where `abs(R_tce) > max_R` to `NA`.
#' @param max_SE Float. Set all entries whwere `SE > max_SE` tp `NA`.
#' @param max_nan_perc Float. Remove columns and rows that are more than
#'   `max_nan_perc` NAs.
#' @export
filter_tce <- function(R_tce, SE_tce, max_R = 1, max_SE = 0.5, max_nan_perc = 0.5) {
  R_tce[is.nan(SE_tce)] <- NA
  SE_tce[is.nan(SE_tce)] <- NA

  R_too_large <- abs(R_tce) > max_R
  R_tce[R_too_large] <- NA
  SE_tce[R_too_large] <- NA
  SE_too_large <- SE_tce > max_SE
  R_tce[SE_too_large] <- NA
  SE_tce[SE_too_large] <- NA

  row_nan_perc <- rowMeans(is.na(R_tce))
  col_nan_perc <- colMeans(is.na(R_tce))
  max_row_nan = max(row_nan_perc)
  max_col_nan = max(col_nan_perc)
  while((max_row_nan > max_nan_perc) | (max_col_nan > max_nan_perc)){
    if(max_row_nan >= max_col_nan){
      which_max_row_nan <- which.max(row_nan_perc)
      R_tce <- R_tce[-which_max_row_nan, -which_max_row_nan]
      SE_tce <- SE_tce[-which_max_row_nan, -which_max_row_nan]
    } else{
      which_max_col_nan <- which.max(col_nan_perc)
      R_tce <- R_tce[-which_max_col_nan, -which_max_col_nan]
      SE_tce <- SE_tce[-which_max_col_nan, -which_max_col_nan]
    }
    row_nan_perc <- rowMeans(is.na(R_tce))
    col_nan_perc <- colMeans(is.na(R_tce))
    max_row_nan = max(row_nan_perc)
    max_col_nan = max(col_nan_perc)
  }
  return(list("R_tce" = R_tce, "SE_tce" = SE_tce))
}
