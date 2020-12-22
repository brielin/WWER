#' Generates sumstats from our bi-directional Mendelian randomization model
#'
#' @param Na Number of samples in each study of phenotype A.
#' @param Nb Number of samples in each study of phenotype B.
#' @param M Total number of SNPs to simulate, including null SNPs.
#' @param sa Variance of the sampling distribution of SNPs effecting A.
#' @param sb Variance of the sampling distribution of SNPs effecting B.
#' @param su Variance of the sampling distribution of SNPs effecting U.
#' @param nu Causal effect of U on A.
#' @param eta Causal effect of U on B.
#' @param gamma Causal effect of A on B.
#' @param delta Causal effect of B on A.
#' @param p Proportion of non-null SNPs.
#' @param q Proportion of non-null SNPs effecting U.
#' @param r Proportion of non-null SNPs effecting A.
#' @param s Proportion of non-null SNPs effecting B.
#' @param af_ref_file String giving the location of a gzipped tsv with a column
#'   called `MAF`, giving a list of reference allele frequencies to use in
#'   simulation. Usually just an HM3 LD-score file for a particular chromosome.
generate_sumstats_bimr_varbeta <- function(Na, Nb, M, sa, sb, su, nu, eta, gamma, delta, p, q, r, s, af_ref_file = NULL){
  # Total effect shorthand.
  R2 = (1/(1-gamma*delta))**2
  Rab2 = gamma**2*R2
  Rba2 = delta**2*R2
  Rua2 = (nu + eta*delta)**2*R2
  Rub2 = (eta + nu*gamma)**2*R2
  gu = M*q*p*su
  ga = M*r*p*sa
  gb = M*s*p*sb

  Qa = 1 - Rua2 - gb*Rba2 - ga*R2
  Qb = 1 - Rub2 - ga*Rab2 - gb*R2

  ea = (R2*Qa - Rba2*Qb)/(R2**2 - Rab2*Rba2)
  eb = (R2*Qb - Rab2*Qa)/(R2**2 - Rab2*Rba2)

  if(ea < 0 | eb < 0){
    stop("Induced environmental variance is less than 0!")
  }

  hu = gu
  ha = hu*Rua2 + ga*R2 + gb*Rba2
  hb = hu*Rub2 + gb*R2 + ga*Rab2

  # TODO(brielin): This should go in a test.
  # 1 - Rua2 = gb*Rba2 + ga*R2 + eb*Rba2 + ea*R2
  # 1 - Rub2 = ga*Rab2 + gb*R2 + ea*Rab2 + eb*R2

  # Indicators for snp->phenotype effects.
  Zu = stats::rbinom(M, 1, p*q)
  Za = stats::rbinom(M, 1, p*r)
  Zb = stats::rbinom(M, 1, p*s)

  # snp->phenotype effect sizes.
  beta_a = stats::rnorm(M, 0, sqrt(sa))
  beta_b = stats::rnorm(M, 0, sqrt(sb))
  beta_u = stats::rnorm(M, 0, sqrt(su))

  # Network effects.
  ba = beta_a*Za*sqrt(R2) + beta_b*Zb*sqrt(Rba2) + beta_u*Zu*sqrt(Rua2)
  bb = beta_b*Zb*sqrt(R2) + beta_a*Za*sqrt(Rab2) + beta_u*Zu*sqrt(Rub2)
  se_beta_a = rep(1/sqrt(Na), M)
  se_beta_b = rep(1/sqrt(Nb), M)
  f = NULL

  # Correct for AF if not null
  if(!is.null(af_ref_file)){
    afs <- as.double(readr::read_tsv(gzfile(af_ref_file))$MAF)
    f = sample(afs, size=M, replace = T)
    af_sd = sqrt(2*f*(1-f))
    ba = ba/af_sd
    bb = bb/af_sd
    se_beta_a <- se_beta_a/af_sd
    se_beta_b <- se_beta_b/af_sd
  }

  # Simulate summary statistics
  ba_hat_1 = stats::rnorm(M, mean = ba, sd = se_beta_a)
  bb_hat_1 = stats::rnorm(M, mean = bb, sd = se_beta_b)
  ba_hat_2 = stats::rnorm(M, mean = ba, sd = se_beta_a)
  bb_hat_2 = stats::rnorm(M, mean = bb, sd = se_beta_b)

  rg = stats::cor(ba, bb)

  beta_hat_1 = cbind(ba_hat_1, bb_hat_1)
  beta_hat_2 = cbind(ba_hat_2, bb_hat_2)
  se_hat = cbind(se_beta_a, se_beta_b)
  beta = cbind(ba, bb)
  Z = cbind(Za, Zb, Zu)
  colnames(beta_hat_1) <- c("A", "B")
  rownames(beta_hat_1) <- paste0("rs", 1:nrow(beta_hat_1))
  colnames(beta_hat_2) <- c("A", "B")
  rownames(beta_hat_2) <- paste0("rs", 1:nrow(beta_hat_2))
  colnames(se_hat) <- c("A", "B")
  rownames(se_hat) <- paste0("rs", 1:nrow(se_hat))
  colnames(beta) <- c("A", "B")
  rownames(beta) <- paste0("rs", 1:nrow(beta))
  colnames(Z) <- c("A", "B", "U")
  rownames(Z) <- paste0("rs", 1:nrow(Z))
  return(list("sumstats_1" = list("beta_hat" = as.data.frame(beta_hat_1), "se_hat" = as.data.frame(se_hat)),
              "sumstats_2" = list("beta_hat" = as.data.frame(beta_hat_2), "se_hat" = as.data.frame(se_hat)),
              "beta" = as.data.frame(beta), "Z" = as.data.frame(Z), "f" = f,
              "hu" = hu, "ha" = ha, "hb" = hb,
              "ga" = ga, "gb" = gb, "ea" = ea, "eb" = eb,
              "sa" = sa, "sb" = sb, "su" = su, "rg" = rg,
              "Rab" = sqrt(Rab2), "Rba" = sqrt(Rba2), "Rua" = sqrt(Rua2), "Rub" = sqrt(Rub2)))
}


#' Generates sumstats from our bi-directional Mendelian randomization model
#'
#' @param Na Number of samples in each study of phenotype A.
#' @param Nb Number of samples in each study of phenotype B.
#' @param M Total number of SNPs to simulate, including null SNPs.
#' @param ha Heritability of phenotype A.
#' @param hb Heritability of phenotype B.
#' @param hu Heritability of phenotype U.
#' @param nu Causal effect of U on A.
#' @param eta Causal effect of U on B.
#' @param gamma Causal effect of A on B.
#' @param delta Causal effect of B on A.
#' @param p Proportion of non-null SNPs.
#' @param q Proportion of non-null SNPs effecting U.
#' @param r Proportion of non-null SNPs effecting A.
#' @param s Proportion of non-null SNPs effecting B.
#' @param af_ref_file String giving the location of a gzipped tsv with a column
#'   called `MAF`, giving a list of reference allele frequencies to use in
#'   simulation. Usually just an HM3 LD-score file for a particular chromosome.
#' @export
generate_sumstats_bimr <- function(Na, Nb, M, ha, hb, hu, nu, eta, gamma, delta, p, q, r, s, af_ref_file = NULL){
  vars <- h2_to_s2(M, ha, hb, hu, nu, eta, gamma, delta, p, q, r, s)
  return(generate_sumstats_bimr_varbeta(Na, Nb, M, vars$sa, vars$sb, vars$su, nu, eta, gamma, delta, p, q, r, s, af_ref_file = af_ref_file))
}

#' Helper function for turning the h2 parameterization into the variance one.
#'
#' @param M Total number of SNPs to simulate, including null SNPs.
#' @param ha Heritability of phenotype A.
#' @param hb Heritability of phenotype B.
#' @param hu Heritability of phenotype U.
#' @param nu Causal effect of U on A.
#' @param eta Causal effect of U on B.
#' @param gamma Causal effect of A on B.
#' @param delta Causal effect of B on A.
#' @param p Proportion of non-null SNPs.
#' @param q Proportion of non-null SNPs effecting U.
#' @param r Proportion of non-null SNPs effecting A.
#' @param s Proportion of non-null SNPs effecting B.
h2_to_s2 <- function(M, ha, hb, hu, nu, eta, gamma, delta, p, q, r, s){
  # Total effect shorthand.
  R2 = (1/(1-gamma*delta))**2
  Rab2 = gamma**2*R2
  Rba2 = delta**2*R2
  Rua2 = (nu + eta*delta)**2*R2
  Rub2 = (eta + nu*gamma)**2*R2

  # Convert heritability to effect size variance.
  Qa = ha - hu*Rua2
  Qb = hb - hu*Rub2
  ga = (R2*Qa - Rba2*Qb)/(R2**2 - Rab2*Rba2)
  gb = (R2*Qb - Rab2*Qa)/(R2**2 - Rab2*Rba2)
  su = ifelse(q == 0, 0, hu/(M*q*p))
  sa = ifelse(r == 0, 0, ga/(M*r*p))
  sb = ifelse(s == 0, 0, gb/(M*s*p))
  return(list(sa = sa, sb = sb, su = su))
}


#' A simple helper function to count the number of instuments for each pair.
#'
#' @param selected A list of lists, output of `select_snps`.
count_instruments <- function(selected) {
  return(do.call(cbind, purrr::map(selected, function(pheno) {
    return(purrr::map(pheno[-1], function(x){sum(x>0, na.rm=T)}))
  })))
}


#' Helper function to create a dataset and compare the results of many methods.
#'
#' @importFrom foreach %do%
#'
#'
#' @param name A name for this simulation setup.
#' @param iter The current iteration.
#' @param Na Number of samples in each study of phenotype A.
#' @param Nb Number of samples in each study of phenotype B.
#' @param M Total number of SNPs to simulate, including null SNPs.
#' @param ha Heritability of phenotype A.
#' @param hb Heritability of phenotype B.
#' @param hu Heritability of phenotype U.
#' @param nu Causal effect of U on A.
#' @param eta Causal effect of U on B.
#' @param gamma Causal effect of A on B.
#' @param delta Causal effect of B on A.
#' @param p Proportion of non-null SNPs.
#' @param q Proportion of non-null SNPs effecting U.
#' @param r Proportion of non-null SNPs effecting A.
#' @param s Proportion of non-null SNPs effecting B.
#' @param af_ref_file String giving the location of a gzipped tsv with a column
#'   called `MAF`, giving a list of reference allele frequencies to use in
#'   simulation. Usually just an HM3 LD-score file for a particular chromosome.
run_sim <- function(name, iter, Na, Nb, M, ha, hb, hu, nu, eta, gamma, delta, p, q, r, s, af_ref_file = NULL){
  dataset <- generate_sumstats_bimr(Na, Nb, M, ha, hb, hu, nu, eta, gamma, delta, p, q, r, s, af_ref_file)
  rg = stats::cor(dataset$beta)[1,2]

  mr_methods = c("egger", "egger_wf", "egger_wa", "egger_s", "ivw", "aps", "raps",
                 "mbe", "median", "mr_mix", "cause", "mr_presso")
  p_thresh = list("5e-6" = 5e-6, "5e-8" = 5e-8)

  sumstats_select <- dataset$sumstats_1
  sumstats_fit <- dataset$sumstats_2

  select_std <- purrr::map(p_thresh, ~ select_snps(sumstats_select, p_thresh = .x, exclusive = FALSE, weight = FALSE, filter = NULL))
  if("cause" %in% mr_methods){
    select_cause <- list("1e-3" = select_snps(sumstats_select, p_thresh = 1e-3, exclusive = FALSE, weight = FALSE, filter = NULL))
  }
  if("egger_wf" %in% mr_methods){
    select_wf <- purrr::map(p_thresh, ~ select_snps(sumstats_select, p_thresh = .x, exclusive = FALSE, weight = TRUE, filter = 1.6))
  }
  if("egger_wa" %in% mr_methods){
    select_wa <- purrr::map(p_thresh, ~ select_snps(sumstats_select, p_thresh = .x, exclusive = FALSE, weight = TRUE, filter = 0))
  }
  if("egger_s" %in% mr_methods){
    select_s <- purrr::map(p_thresh, ~ select_snps(sumstats_select, p_thresh = .x, exclusive = FALSE, weight = FALSE, filter = 1.6))
  }

  select_res <- foreach::foreach(mr_method = mr_methods, .combine = dplyr::bind_rows) %do% {
    selected <- switch(mr_method, cause = select_cause, egger_wf = select_wf, egger_wa = select_wa, egger_s = select_s, select_std)
    pv_res <- foreach::foreach(selected = selected, p_thresh = names(selected), .combine = dplyr::bind_rows) %do% {
      print(c(iter, mr_method, p_thresh))
      inst_count <- count_instruments(selected)
      start.time <- proc.time()[[3]]
      tce_res <- tryCatch(
        {
          fit_tce(sumstats_fit, selected, mr_method = switch(mr_method, egger_wf = "egger_w", egger_wa = "egger_w", egger_s = "egger", mr_method))
        }, error = function(cond){
          return(list(R_tce = matrix(nrow=2, ncol=2), SE_tce = matrix(nrow=2, ncol=2)))
        })
      runtime <- proc.time()[[3]] - start.time
      tibble::tibble(
        name = name, iter=iter, mr_method = mr_method, p_thresh = p_thresh, runtime = runtime,
        Rab = dataset$Rab, Rba = dataset$Rba,
        Rab_hat = tce_res$R_tce[1, 2], Rba_hat = tce_res$R_tce[2, 1],
        Sab_hat = tce_res$SE_tce[1, 2], Sba_hat = tce_res$SE_tce[2, 1],
        pab = tce_res$p_val[1,2], pba = tce_res$p_val[2,1],
        Na = Na, Nb = Nb, M = M, ha = dataset$ha, hb = dataset$hb, hu = dataset$hu,
        ga = dataset$ga, gb = dataset$gb, ea = dataset$ea, eb = dataset$eb,
        nu = nu, eta = eta, gamma = gamma, delta = delta, q = q, r = r, s = s,
        nab = inst_count[1,2][[1]], nba = inst_count[2,1][[1]], rg = rg)
    }
  }
  select_res <- tidyr::unite(
    select_res, "method", c("mr_method", "p_thresh"), sep = "_", remove = FALSE)
  return(select_res)
}


