#' @title
#' Analysis of the simulated data according to the cross-validate adaptive enrichment risk scores (CADEN) design
#'
#' @description The design has two stages. At the end of stage 1, an interim analysis is performed to test the efficacy of the treatment compared to the control in the overall trial population. If the test in the overall population is significanct, the trial proceeds into stage 2 by enrolling all patients (the "unselected" strategy) implemented here via simulating the data for stage 2 patients. If the test in the overall population is not significant then the trial proceeds by testing whether there is a subgroup of patients (sensitive group) that show a promising treatment effect. Depending on the results of the test in the sensitive group, the design proceeds into stage 2 according to one of the mutually exclusive strategies "stop" or "enrichment".
#'
#' @param param_txt A name of the text file with the parameters as described in simulate.R
#' @param simdata_stage1 A list with the simulated data according to the output of simulate.R
#'
#' @return An object of class \code{"caden"}.
#' @return datalist_stage2: a list of simulated data for stage 2
#' @return sens_status_predicted: A vector of the predicted sensitivity status
#' @return sensitivity: A vector of values for sensitivity for identifying the sensitive group
#' @return specificity: A vector of values for specificity for identifying the sensitive group
#' @return noneligible: Number of patients non-eligible for the trial (for the enrichment strategy).
#' @return pval_overall: P-value for the difference between the arm in the overall trial population (for "unselected" strategy only)
#' @return pval_fisher_group: P-value for the difference between the arms in the sensitive group (as assessed using Fisher exact test)
#'
#' @examples
#' data(simdata_stage1)
#' param_file <- "data/param.txt"
#' sim_res <- analyse_simdata(param_file, simdata_stage1)
#' @seealso
#' \code{\link{simulate_data}}, function.
#' @author Svetlana Cherlin, James Wason
#' @export analyse_simdata
#' @importFrom "poolr" "stouffer"

analyse_simdata <- function(param_file, datalist_stage1) {
  datalist <- NULL
  datalist_stage2 <- NULL
  sensitivity <- NA
  specificity <- NA
  sensitivity_stage2 = NA
  specificity_stage2 = NA
  noneligible <- NA
  pval_sens_group <- NA
  pval_sens_group_stage2 <- NA
  pval_overall <- NA
  ## get simulation parameters
  input <- read.table(param_file, header = F, fill = TRUE)
  input <- na.omit(input)
  cvrs <- NA
  cvrs_stage1 <- NA
  cvrs_stage2 <- NA

  size_stage1 <- as.numeric(input[input == "size_stage1", 2])
  size_stage2 <- as.numeric(input[input == "size_stage2", 2])
  num_all_var <- as.numeric(input[input == "num_all_var", 2])
  num_sens_var <- as.numeric(input[input == "num_sens_var", 2])
  mu1 <- as.numeric(input[input == "mu1", 2])
  mu2 <- as.numeric(input[input == "mu2", 2])
  mu0 <- as.numeric(input[input == "mu0", 2])
  sigma1 <- as.numeric(input[input == "sigma1", 2])
  sigma2 <- as.numeric(input[input == "sigma2", 2])
  sigma0 <- as.numeric(input[input == "sigma0", 2])
  rho1 <- as.numeric(input[input == "rho1", 2])
  rho2 <- as.numeric(input[input == "rho2", 2])
  rho0 <- as.numeric(input[input == "rho0", 2])
  perc_sens <- as.numeric(input[input == "perc_sens", 2])
  resp_rate_treat <- as.numeric(input[input == "resp_rate_treat", 2])
  resp_rate_con <- as.numeric(input[input == "resp_rate_con", 2])
  resp_rate_sens_treat <- as.numeric(input[input == "resp_rate_sens_treat", 2])
  seed = input[input == "seed", 2]
  if (length(seed) == 0) {seed = as.null()}
  if (!is.null (seed))
  {
     seed <- as.numeric(input[input == "seed", 2])
  }
  threshold_overall <- as.numeric(input[input == "threshold_overall", 2])
  threshold_group <- as.numeric(input[input == "threshold_group", 2])
  standardise_cvrs <- as.logical(input[input == "standardise_cvrs", 2])
  full_model <- as.logical(input[input == "full_model", 2])

  print("Input parameters:")
  print(paste("size_stage1: ", size_stage1, sep = ""))
  print(paste("size_stage2: ", size_stage2, sep = ""))
  print(paste("num_all_var: ", num_all_var, sep = ""))
  print(paste("num_sens_var: ", num_sens_var, sep = ""))
  print(paste("mu1: ", mu1, sep = ""))
  print(paste("mu2: ", mu2, sep = ""))
  print(paste("mu0: ", mu0, sep = ""))
  print(paste("sigma1: ", sigma1, sep = ""))
  print(paste("sigma2: ", sigma2, sep = ""))
  print(paste("sigma0: ", sigma0, sep = ""))
  print(paste("rho1: ", rho1, sep = ""))
  print(paste("rho2: ", rho2, sep = ""))
  print(paste("rho0: ", rho0, sep = ""))
  print(paste("perc_sens: ", perc_sens, sep = ""))
  print(paste("resp_rate_treat: ", resp_rate_treat, sep = ""))
  print(paste("resp_rate_con: ", resp_rate_con, sep = ""))
  print(paste("resp_rate_sens_treat: ", resp_rate_sens_treat, sep = ""))
  print(paste("seed: ", seed, sep = ""))
  print(paste("threshold_overall: ", threshold_overall, sep = ""))
  print(paste("threshold_group: ", threshold_group, sep = ""))
  print(paste("standardise_cvrs: ", standardise_cvrs, sep = ""))
  print(paste("full_model: ", full_model, sep = ""))

  ## sensitive covariates for sensitive patients
  m1 <- numeric(length = num_sens_var)
  m1[] <- mu1
  sigma1_mtx <- matrix(nrow = num_sens_var, ncol = num_sens_var, data = sigma1 * sigma1 * rho1)
  diag(sigma1_mtx) <- sigma1^2

  ## sensitive covariates for non-sensitive patients
  m2 <- numeric(length = num_sens_var)
  m2[] <- mu2
  sigma2_mtx <- matrix(nrow = num_sens_var, ncol = num_sens_var, data = sigma2 * sigma2 * rho2)
  diag(sigma2_mtx) <- sigma2^2

  ## non-sensitive covariates for  all patients
  m0 <- numeric(length = num_all_var - num_sens_var)
  m0[0] <- mu0
  sigma0_mtx <- matrix(nrow = num_all_var - num_sens_var, ncol = num_all_var - num_sens_var, data = sigma0 * sigma0 * rho0)
  diag(sigma0_mtx) <- sigma0^2

  mu <- log(resp_rate_con / (1 - resp_rate_con)) # intercept corresponding to response rate for controls
  lambda <- log(resp_rate_treat / (1 - resp_rate_treat)) - mu # main treatment effect
  interaction.scaling <- log(resp_rate_sens_treat / (1 - resp_rate_sens_treat)) - mu - lambda # interaction scaling
  gamma <- sapply(m1, function(x) {
    ifelse(x == 0, 0, interaction.scaling / (num_sens_var))
  }) ## covariate-treatment interation term

  params <- list(num_all_var = num_all_var, num_sens_var = num_sens_var, size_stage1 = size_stage1, size_stage2 = size_stage2, m1 = m1, m2 = m2, m0 = m0, sigma1_mtx = sigma1_mtx, sigma2_mtx = sigma2_mtx, sigma0_mtx = sigma0_mtx, perc_sens = perc_sens, mu = mu, lambda = lambda, gamma = gamma, seed = seed)
  
  ## Stage 1
  
  res1 <- analyse_stage1(datalist_stage1, threshold_overall, threshold_group, seed, standardise_cvrs, full_model)
  sens_status_predicted <- res1$sens_status_predicted
  sensitivity_stage1 = res1$sensitivity
  specificity_stage1 = res1$specificity
  cvrs_stage1 <- res1$cvrs
  
  ## Stage 2
  
  if (res1$decision != "stop") {
    res2 <- analyse_stage2(res1$decision, params, res1$model, standardise_cvrs)
    noneligible <- res2$noneligible
    datalist_stage2 <- res2$datalist_stage2
    cvrs_stage2 <- res2$cvrs
  }

  ## Final analysis
  
  # for unselected, compute the overall effect and the effect in the sensitive group
  if (res1$decision == "unselected") {
    datalist <- list(patients = rbind(datalist_stage1$patients, datalist_stage2$patients), 
                     covar = rbind(datalist_stage1$covar, datalist_stage2$covar), 
                     resp_rate = c(datalist_stage1$resp_rate, datalist_stage2$resp_rate), 
                     response = c(datalist_stage1$response, datalist_stage2$response))
    out <- get_subgroup(datalist, seed, standardise_cvrs, full_model)
    sens_status_predicted <- out$sens_status_predicted
    sensitivity <- out$sensitivity
    specificity <- out$specificity
    cvrs = out$cvrs
    standardise_cvrs = out$standardise_cvrs
    size_treat <- nrow(datalist$patients[datalist$patients$treat == 1, ])
    size_con <- nrow(datalist$patients[datalist$patients$treat == 0, ])
    pval_overall <- prop.test(
      x = c(sum(datalist$response[datalist$patients$treat == 0]), sum(datalist$response[datalist$patients$treat == 1])),
      n = c(size_con, size_treat), alternative = "two.sided"
    )$p.value
    pval_sens_group <- analyse_fisher(datalist, sens_status_predicted)
  }

  ## for enrichment, analyse 2 stages separately and combine the p-values
  if (res1$decision == "enrichment") {
    sensitivity_stage2 = res2$sensitivity
    specificity_stage2 = res2$specificity
    sens_status_predicted <- c(sens_status_predicted, res2$sens_status_predicted)
    pval_sens_group_stage2 <- analyse_fisher(datalist_stage2, res2$sens_status_predicted)
    pval_sens_group = stouffer(c(res1$pval_sens_group, pval_sens_group_stage2))[[1]]
  }
  
  output <- list(
    datalist_stage2 = datalist_stage2, decision = res1$decision,
    sens_status_predicted = sens_status_predicted,
    sensitivity_stage1 = sensitivity_stage1,
    specificity_stage1 = specificity_stage1,
    sensitivity_stage2 = sensitivity_stage2,
    specificity_stage2 = specificity_stage2,
    sensitivity = sensitivity, specificity = specificity,
    noneligible = noneligible, pval_overall = pval_overall,
    pval_overall_stage1 = res1$pval_overall,
    cvrs = cvrs, cvrs_stage1 = cvrs_stage1, cvrs_stage2 = cvrs_stage2,
    pval_sens_group_stage1 = res1$pval_sens_group,
    pval_sens_group_stage2 = pval_sens_group_stage2,
    pval_sens_group = pval_sens_group
  )

  class(output) <- "caden"

  return(output)
}

