#' @title
#' Simulate data for the cross-validated adaptive enrichment risc scores (CADEN) design.
#'
#' @description
#' The data is simulated assuming that the response to treatment is influenced by a subset of K unknown covariates (the sensitive covariates) through the following model:
#'
#' logit(p_i)= mu+lambda*t_i+gamma_1*t_i*x_i1+...+gamma_K*t_i*x_iK,
#'
#' where p_i is the probability of response to treatment for the i-th patient; mu is the intercept; lambda is the treatment main effect that all patients experience regardless of the values of the covariates; t_i is the treatment that the i-th patient receives (t_i = 0 for the control arm and t_i=1 for the treatment arm); x_i1,...,x_iK are the values for the K unknown sensitive covariates; gamma_1,...,gamma_K are treatment-covariate interaction effects for the K covariates.
#' The model assumes that there is a subset of patients (the sensitive group) with a higher probability of response when treated with the new treatment, compared with the control treatment.
#'
#'
#' @param param_file A name of the parameters' text file. The file sould have a row for each parameter with the name of the parameters followed by a space followed by a value of the parameter. The list of the parameters is as follows:
#'
#'   size_stage1 - sample size for stage 1
#'
#'   size_stage2 - sample size for stage 2
#'
#'   num_all_var - number of covariates
#'
#'   num_sens_var - number of sensitive covariates
#'
#'   mu1 - mean for sensitive covariates in the sensitive group
#'
#'   mu2 - mean for the sensitive covariates in the non-sensitive  group
#'
#'   mu0 - mean for the non-senstive covariates
#'
#'   sigma1 - sd for sensitive covariates in the sensitive group
#'
#'   sigma2 - sd for the sensitive covariates in the non-sensitive group
#'
#'   sigma0 - sd for the non-senstive covariates
#'
#'   rho1 - correlation for sensitive covariates in the sensitive group
#'
#'   rho2 - correlation for the sensitive covariates in the non-sensitive group
#'
#'   rho0 - correlation for the non-senstive covariates
#'
#'   perc_sens - prevalence of the sensitive group
#'
#'   resp_rate_treat  - response rate for everyone on treatment
#'
#'   resp_rate_con  - response rate for everyone on control
#'
#'   resp_rate_sens_treat  - response rate for the sensitive group on treatment
#'
#'   seed  - seed for random number generating
#'
#'   threshold_overall - p-value threshold for the test for the differences in the treatment effect in the overall trial population
#'
#'   threshold_group - p-value threshold for the test for the treatment effect in the sensitive group
#'   
#'   standardise_cvrs - A logical flag (0/1) for the standardisation of the risk scores during the analysis. The standardisation is performed with respect to the training data sets, per cross-validation fold
#'
#'   full_model - A logical flag for the full model (treatment effect, covariate effect and the interaction effect) for the analysis. When full_model is 0, only interaction effect is included in the model
#'
#'   The example of the file is as follows (see also the example of the file in the "data" directory).
#'
#'   size_stage1 500
#'
#'   size_stage2 500
#'
#'   num_all_var 100
#'
#'   num_sens_var 10
#'
#'   mu1 1
#'
#'   mu2 0
#'
#'   mu0 0
#'
#'   sigma1 0.5
#'
#'   sigma2 0.1
#'
#'   sigma0 0.5
#'
#'   rho1 0
#'
#'   rho2 0
#'
#'   rho0 0
#'
#'   perc_sens 0.1
#'
#'   resp_rate_treat 0.25
#'
#'   resp_rate_con 0.25
#'
#'   resp_rate_sens_treat 0.7
#'
#'   seed 123
#'
#'   threshold_overall 0.04
#'
#'   threshold_group 0.01
#' 
#'   standardise_cvrs 0
#'
#'   full_model 0

#'
#' @return A list with two data frames (patients, covar) and two vectors (resp.rate, response)
#' @return patients: a data frame with one row per patient and the following columns:
#'         FID (family ID), IID (individual ID), treat (1 for treatment and 0 for control),
#'         sens.status (true sensitivity status), stage (1)
#' @return covar: a data frame with covariate data for L covariates
#' @return resp.rate: a vector of response rates
#' @return response: a vector of simulated binary responses
#'
#'
#' @examples
#' param_file = "data/param.txt"
#' simdata_stage1 <- simulate_data(param_file)
#' @seealso
#' \code{\link[pkg]{analyse_simdata.R}} function.
#' @author Svetlana Cherlin, James Wason
#' @export simulate_data
#' @importFrom "MASS" "mvrnorm"



simulate_data <- function(param_file) {
  input <- read.table(param_file, header = F)

  size_stage1 <- input[input == "size_stage1", 2]
  num_all_var <- input[input == "num_all_var", 2]
  num_sens_var <- input[input == "num_sens_var", 2]
  mu1 <- input[input == "mu1", 2]
  mu2 <- input[input == "mu2", 2]
  mu0 <- input[input == "mu0", 2]
  sigma1 <- input[input == "sigma1", 2]
  sigma2 <- input[input == "sigma2", 2]
  sigma0 <- input[input == "sigma0", 2]
  rho1 <- input[input == "rho1", 2]
  rho2 <- input[input == "rho2", 2]
  rho0 <- input[input == "rho0", 2]
  perc_sens <- input[input == "perc_sens", 2]
  resp_rate_treat <- input[input == "resp_rate_treat", 2]
  resp_rate_con <- input[input == "resp_rate_con", 2]
  resp_rate_sens_treat <- input[input == "resp_rate_sens_treat", 2]
  seed <- input[input == "seed", 2]

  print("Input parameters:")
  print(paste("size_stage1: ", size_stage1, sep = ""))
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

  mu <- log(resp_rate_con / (1 - resp_rate_con)) # intercept corresponding to response rate for controls
  lambda <- log(resp_rate_treat / (1 - resp_rate_treat)) - mu # main treatment effect
  interaction.scaling <- log(resp_rate_sens_treat / (1 - resp_rate_sens_treat)) - mu - lambda # interaction scaling

  # sensitive covariates for the sensitive group
  m1 <- numeric(length = num_sens_var)
  m1[] <- mu1
  sigma1_mtx <- matrix(nrow = num_sens_var, ncol = num_sens_var, data = sigma1 * sigma1 * rho1)
  diag(sigma1_mtx) <- sigma1^2

  # sensitive covariates for non-sensitive patients
  m2 <- numeric(length = num_sens_var)
  m2[] <- mu2
  sigma2_mtx <- matrix(nrow = num_sens_var, ncol = num_sens_var, data = sigma2 * sigma2 * rho2)
  diag(sigma2_mtx) <- sigma2^2

  # non-sensitive covariates, all patients
  m0 <- numeric(length = num_all_var - num_sens_var)
  m0[0] <- mu0
  sigma0_mtx <- matrix(nrow = num_all_var - num_sens_var, ncol = num_all_var - num_sens_var, data = sigma0 * sigma0 * rho0)
  diag(sigma0_mtx) <- sigma0^2

  mu <- log(resp_rate_con / (1 - resp_rate_con)) # intercept corresponding to response rate for controls
  lambda <- log(resp_rate_treat / (1 - resp_rate_treat)) - mu # main treatment effect
  interaction.scaling <- log(resp_rate_sens_treat / (1 - resp_rate_sens_treat)) - mu - lambda # interaction scaling
  gamma <- sapply(m1, function(x) {
    ifelse(x == 0, 0, interaction.scaling / (num_sens_var * x))
  }) ## covariate-treatment interation term

  sens_status <- numeric(length = size_stage1)

  ## Simulate covariates for the sensitive group
  covar_sens <- cbind(
    mvrnorm(n = size_stage1 * perc_sens, m1, sigma1_mtx, tol = 1e-6), # sensitive covariates
    mvrnorm(n = size_stage1 * perc_sens, m0, sigma0_mtx, tol = 1e-6)
  ) # non-sensitive covariates

  covar_sens <- as.data.frame(covar_sens)
  sens_status[1:size_stage1 * perc_sens] <- 1 # subgroup membership status (sensitivity status)

  ## Simulate covariates for patients who do not belong to the sensitive subgroup
  covar_nosens <- cbind(
    mvrnorm(n = size_stage1 * (1 - perc_sens), m2, sigma2_mtx, tol = 1e-6), # sensitive covariates
    mvrnorm(n = size_stage1 * (1 - perc_sens), m0, sigma0_mtx, tol = 1e-6)
  ) # non-sensitive covariates
  covar_nosens <- as.data.frame(covar_nosens)
  sens_status[(size_stage1 * perc_sens) + 1:size_stage1] <- 0 # subgroup membership status

  ## Combine the covariates
  covar <- as.data.frame(rbind(covar_sens, covar_nosens))
  colnames(covar) <- paste("Covar", 1:num_all_var, sep = "")
  covar <- as.matrix(covar)

  ## Equal randomisation to control/treatment arm
  treat <- seq(size_stage1) %% 2 # alternates between 1 and 0, N times

  ## Compute response rate
  levels <- covar[, 1:num_sens_var] %*% as.matrix(gamma, nrow = num_sens_var) ## covariate-treatment interation term multiplied by the value of the covariate
  levels <- levels[, 1]
  linpred <- mu + treat * lambda + treat * levels
  resp_rate <- as.vector(exp(linpred) / (1 + exp(linpred)))

  ## Simulate response
  response <- rbinom(size_stage1, 1, resp_rate)

  ## Sort the data according to the treatment assignment (treatment first)
  covar <- as.matrix(covar[order(treat, decreasing = TRUE), ])
  rownames(covar) <- seq(nrow(covar))
  resp_rate <- resp_rate[order(treat, decreasing = TRUE)]
  response <- response[order(treat, decreasing = TRUE)]
  sens_status <- sens_status[order(treat, decreasing = TRUE)]
  treat <- treat[order(treat, decreasing = TRUE)]

  patients <- data.frame(FID = seq(1:size_stage1), IID = seq(1:size_stage1), treat = treat, sens_status = sens_status, stage = rep(1, size_stage1))

  return(list(patients = patients, covar = covar, resp_rate = resp_rate, response = response))
}
