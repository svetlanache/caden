
#############################################################
#' @title analyse_stage1
#'
#' @description Perform interim analysis
#'             
#'
#' @param  datalist_stage1
#'         threshold_overall
#'         threshold_group
#'         standardise_cvrs
#'         full_model
#'
#' @return  A list with the following elements
#'          decision
#'          sens_status_predicted
#'          model 
#'
#' @author Svetlana Cherlin, James Wason
#############################################################

analyse_stage1 <- function(datalist_stage1, threshold_overall, threshold_group, seed, standardise_cvrs, full_model) {
  patients <- datalist_stage1$patients
  response <- datalist_stage1$response
  model <- NULL
  sens_status_predicted <- NULL
  sensitivity <- NA
  specificity <- NA
  pval_sens_group = NA
  cvrs <- NA

  ## power for the overall arm comparison
  size_treat <- nrow(patients[patients$treat == 1, ])
  size_con <- nrow(patients[patients$treat == 0, ])
  pval_overall <- prop.test(
    x = c(sum(response[patients$treat == 0]), sum(response[patients$treat == 1])),
    n = c(size_con, size_treat), alternative = "two.sided"
  )$p.value

  if (pval_overall > threshold_overall) {
    res <- get_subgroup(datalist_stage1, seed, standardise_cvrs, full_model)
    sensitivity = res$sensitivity
    specificity = res$specificity
    sens_status_predicted <- res$sens_status_predicted
    pval_sens_group <- analyse_fisher(datalist_stage1, sens_status_predicted)
    if (pval_sens_group < threshold_group) {
      model <- get_model(datalist_stage1, standardise_cvrs, full_model)
      decision <- "enrichment"
    } else {
      decision <- "stop"
    }
  } else {
    decision <- "unselected"
  }
  output <- list(decision = decision, sens_status_predicted = sens_status_predicted, 
                 model = model, sensitivity = sensitivity, specificity = specificity,
                 pval_overall = pval_overall, pval_sens_group = pval_sens_group, cvrs = res$cvrs)
  return(output)
}


#############################################################
#' @title get_model
#'
#' @description Compute the coefficients for the risk scores.
#'             
#'
#' @param  datalist_stage1 
#'         threshold_overall
#'         threshold_group
#'         standardise_cvrs
#'         full_model
#'
#' @return  A list with the following elements
#'          decision
#'          sens_status_predicted
#'          model 
#'
#' @author Svetlana Cherlin, James Wason
#############################################################


get_model <- function(datalist, standardise_cvrs, full_model) {
  mean_cvrs <- NA
  sd_cvrs <- NA
  beta_hat <- apply(as.matrix(datalist$covar), 2, function(x) {
    mod <- tryCatch({
        if (full_model) {
           glm(datalist$response ~ datalist$patients$treat + x + datalist$patients$treat:x)
        } else {
           glm(datalist$response ~ datalist$patients$treat:x)
        }
      },
      error = function(e) e,
      warning = function(w) w
    )
    if (is(mod, "error") | is(mod, "warning")) {
      0
    } else {
      if (is.na(mod$coeff[names(mod$coeff) == "datalist$patients$treat:x"])) {
        0
      } else {
        mod$coeff[names(mod$coeff) == "datalist$patients$treat:x"]
      }
    }
  })
  beta_hat[is.na(beta_hat)] <- 0

  cvrs <- apply(as.matrix(datalist$covar), 1, function(x) {
    sum(x * beta_hat)
  })

  if (standardise_cvrs) {
    mean_cvrs <- mean(cvrs)
    sd_cvrs <- sd(cvrs)
    cvrs <- (cvrs - mean_cvrs) / sd_cvrs
  }

  km <- tryCatch({
      kmeans(cvrs, 2)
    },
    error = function(e) e,
    warning = function(w) w
  )
  if (is(km, "error") | is(km, "warning")) {
    message(km)
    centers <- NULL
 } else {
    centers <- sort(km$centers, decreasing = TRUE)
 }

  model <- list(beta_hat = beta_hat, centers = centers, mean_cvrs = mean_cvrs, sd_cvrs = sd_cvrs)
  return(model)
}

#############################################################
#' @title check_eligibility
#'
#' @description Compute the eligibility status for stage 2
#'             
#'
#' @param  covar
#'         model
#'         standardise_cvrs
#'
#' @return  A list with the following elements
#'          eligible
#'          sens_status_predicted
#'          cvrs
#'
#' @author Svetlana Cherlin, James Wason
#############################################################

check_eligibility <- function(covar, model, standardise_cvrs) {
  cvrs <- sum(covar * model$beta_hat) # risk score
  if (standardise_cvrs) {
    cvrs <- (cvrs - model$mean_cvrs) / model$sd_cvrs
  }
  dist <- abs(cvrs - model$centers) # distance from cluster centers

  if (dist[1] < dist[2]) {
    eligible <- 1
  } else {
    eligible <- 0
  } # If new risk score is close to the first cluster then eligible = TRUE
  sens_status_predicted <- eligible

  output <- list(eligible = eligible, sens_status_predicted = sens_status_predicted, cvrs = cvrs)
  return(output)
}

#############################################################
#' @title analyse_subgroup
#'
#' @description Compute the treatment effect in the sensitive group
#'
#' @param  datalist 
#'         sens_status_predicted 
#'
#' @return  pval_treat_group P-value for the treatment effect in the sensitive group
#'
#' @author Svetlana Cherlin, James Wason
#############################################################

# analyse_subgroup <- function(datalist, sens_status_predicted) {
#   mod <- tryCatch({
#       glm(datalist$response ~ datalist$patients$treat + sens_status_predicted + datalist$patients$treat:sens_status_predicted)
#     },
#     error = function(e) e,
#     warning = function(w) w
#   )
# 
#   if (is(mod, "error") | is(mod, "warning")) {
#     pval_treat_group <- 1
#   } else {
#     theta <- as.vector(mod$coeff)
#     temp <- c(0, 1, 0, 1)
#     treat_group <- as.numeric(temp %*% theta) ## treatment effect in the sensitive group
#     var_treat_group <- as.numeric(t(as.matrix(temp)) %*% vcov(mod) %*% as.matrix(temp)) ## sd of the treatment effect in the sensitive group
#     statistic_treat_group <- treat_group / (sqrt(var_treat_group))
#     pval_treat_group <- 2 * pnorm(-abs(statistic_treat_group))
#   }
#   
#   print(paste("pval_treat_group is ", pval_treat_group, sep = ""))
#   return(pval_treat_group)
# }

#############################################################
#' @title get_subgroup
#'
#' @description Compute the subgroup memebership, sensitivity, 
#'              specificity and risk scores.
#'
#' @param  datalist 
#'         seed 
#'         standardise_cvrs
#'         full_model
#'
#' @return  A list with the following elements
#'          cvrs Risk scores
#'          sens_status_predicted
#'          sensitivity
#'          specificity
#'
#' @author Svetlana Cherlin, James Wason
#############################################################

get_subgroup <- function(datalist, seed, standardise_cvrs, full_model) {
  sensitivity <- NA
  specificity <- NA
  patients <- datalist$patients
  covar <- datalist$covar
  response <- datalist$response
  if (!is.null(seed)) {
    set.seed(seed)
  }

  ## Divide to folds and keep prevalence of responders/non-responders within folds
  nfolds <- 10
  patients_nr <- patients[response == 0, ]
  patients_r <- patients[response == 1, ]

  covar_nr <- covar[response == 0, ]
  covar_r <- covar[response == 1, ]

  response_nr <- response[response == 0]
  response_r <- response[response == 1]

  foldid_nr <- sample(rep(seq(nfolds), length = length(response_nr))) # non-responders
  foldid_r <- sample(rep(seq(nfolds), length = length(response_r))) # responders
  res <- data.frame()

  ## CV loop
  for (i in seq(nfolds)) {
    which_nr <- foldid_nr == i
    which_r <- foldid_r == i
    patients_train <- rbind(patients_nr[!which_nr, ], patients_r[!which_r, ])
    patients_test <- rbind(patients_nr[which_nr, ], patients_r[which_r, ])
    covar_train <- rbind(covar_nr[!which_nr, ], covar_r[!which_r, ])
    covar_test <- rbind(covar_nr[which_nr, ], covar_r[which_r, ])
    response_train <- c(response_nr[!which_nr], response_r[!which_r])
    ## Fit a single-covariate regression model for the training data
    beta_hat <- apply(as.matrix(covar_train), 2, function(x) {
      mod <- tryCatch({
          if (full_model) {
             glm(response_train ~ patients_train$treat + x + patients_train$treat:x)
          } else {
             glm(response_train ~ patients_train$treat:x)
          }
        },
        error = function(e) e,
        warning = function(w) w
      )
      if (is(mod, "error") | is(mod, "warning")) {
        0
      } else {
        if (is.na(mod$coeff[names(mod$coeff) == "patients_train$treat:x"])) {
          0
        } else {
          mod$coeff[names(mod$coeff) == "patients_train$treat:x"]
        }
      }
    })

    beta_hat[is.na(beta_hat)] <- 0

    ## compute risk scores from training data
    cvrs_train = apply(as.matrix(covar_train), 1, function (x)
      {
         sum(x*beta_hat)
      })

    ## Compute the  risk scores in the testing data
    cvrs <- apply(as.matrix(covar_test), 1, function(x) {
      sum(x * beta_hat)
    })

    ## Standardize cvrs according to the training data set
    if (standardise_cvrs) {
      cvrs <- (cvrs - mean(cvrs_train)) / sd(cvrs_train)
    }

    temp <- data.frame(FID = patients_test$FID, IID = patients_test$IID, cvrs = cvrs, sens_status_predicted = 0)

    ## Divide testing data to 2 clusters
    km <- tryCatch({
        kmeans(temp$cvrs, 2)
      },
      error = function(e) e,
      warning = function(w) w
    )

    if (is(km, "error") | is(km, "warning")) {
      message(km)
    } else {
      if (km$centers[1] > km$centers[2]) {
        temp$sens_status_predicted[km$cluster == 1] <- 1
      } else {
        temp$sens_status_predicted[km$cluster == 2] <- 1
      }
      res <- rbind(res, temp)
    }
  } # End of CV loop

  m <- match(
    paste(as.character(patients$FID), as.character(patients$IID), sep = ":"),
    paste(as.character(res$FID), as.character(res$IID), sep = ":")
  )
  res <- res[m, ]
  ## Compute sensitivity and specificity of the sensitive group selection algorithm for simulated data
  if (with(patients, exists('sens_status'))) {
    conf <- matrix(
      nrow = 2, ncol = 2,
      data = c(
        sum(!res$sens_status_predicted & !patients$sens_status), # predicted non.sens, true non.sens [1,1]
        sum(!res$sens_status_predicted & patients$sens_status), # predicted non.sens, true sens [1,2]
        sum(res$sens_status_predicted & !patients$sens_status), # predicted sens, true non. sens [2,1]
        sum(res$sens_status_predicted & patients$sens_status) # predicted sens, true sens [2,2]
      ), 
      byrow = TRUE
    )
    sensitivity <- conf[2, 2] / (conf[2, 2] + conf[1, 2])
    specificity <- conf[1, 1] / (conf[1, 1] + conf[2, 1])
  }

  ret <- list(cvrs = res$cvrs, sens_status_predicted = res$sens_status_predicted, sensitivity = sensitivity, specificity = specificity)
  return(ret)
}

#############################################################
#' @title analyse_stage2
#'
#' @description Analyse the stage 2 of the trial
#'              
#'
#' @param  decision
#'         params
#'         model
#'         standardise_cvrs
#'
#' @return  A list with the following elements
#'          datalist_stage2
#'          noneligible
#'          sensitivity
#'          specificity
#'          sens_status_predicted
#'
#' @author Svetlana Cherlin, James Wason
#############################################################

analyse_stage2 <- function(decision, params, model, standardise_cvrs) {
  noneligible <- 0
  response <- numeric()
  covar <- matrix(ncol = params$num_all_var, nrow = 0)
  resp_rate <- numeric(length = params$size_stage2)
  treat <- numeric(length = params$size_stage2)
  sens_status_predicted <- NULL
  sens_status <- NULL
  noneligible <- 0
  conf <- NULL
  sensitivity <- NA
  specificity <- NA
  eligible <- NULL
  cvrs <- NA
  i <- 0
  while (i < params$size_stage2) {

    ## simulate sensitivity status and covariates
    sens_status <- c(sens_status, rbinom(1, 1, params$perc_sens)) ## simulate sensitivity status
    if (sens_status[length(sens_status)]) {
      covar_simulated <- c(mvrnorm(n = 1, params$m1, params$sigma1_mtx, tol = 1e-6), mvrnorm(n = 1, params$m0, params$sigma0_mtx, tol = 1e-6))
    } else {
      covar_simulated <- c(mvrnorm(n = 1, params$m2, params$sigma2_mtx, tol = 1e-6), mvrnorm(n = 1, params$m0, params$sigma0_mtx, tol = 1e-6))
    }

    ## check eligibility
    if (decision == "unselected") {
      eligible <- c(eligible, 1)
      sens_status_predicted <- c(sens_status_predicted, NA)
    } else if (decision == "enrichment") {
      res <- check_eligibility(covar_simulated, model, standardise_cvrs)
      eligible <- c(eligible, res$eligible)
      sens_status_predicted <- c(sens_status_predicted, res$sens_status_predicted)
      if (res$eligible){
        cvrs <- c(cvrs, res$cvrs)
      }
    }

    if (eligible[length(eligible)]) {
      i <- i + 1
      covar <- rbind(covar, covar_simulated)
      treat[i] <- rbinom(1, 1, 0.5) ## randomise
      levels <- covar[i, 1:params$num_sens_var] %*% as.matrix(params$gamma, nrow = params$num_sens_var) ## covariate-treatment interation term multiplied by the value of the covariate
      linpred <- params$mu + treat[i] * params$lambda + treat[i] * levels
      resp_rate[i] <- as.vector(exp(linpred) / (1 + exp(linpred))) ## compute response rate
      response[i] <- rbinom(1, 1, resp_rate[i]) ## simulate response
    } else {
      noneligible <- noneligible + 1
    }
  } ## end while loop
  names(covar) <- NULL
  ## For enrichment strategy, compute sensitivity and specificity of the sensitive group
  if (decision == "enrichment") {
      conf <- matrix(
      nrow = 2, ncol = 2,
      data = c(
        sum(!sens_status_predicted & !sens_status), # predicted non.sens, true non.sens [1,1]
        sum(!sens_status_predicted & sens_status), # predicted non.sens, true sens [1,2]
        sum(sens_status_predicted & !sens_status), # predicted sens, true non. sens [2,1]
        sum(sens_status_predicted & sens_status)
      ), # predicted sens, true sens [2,2]
      byrow = TRUE
    )
    sensitivity <- conf[2, 2] / (conf[2, 2] + conf[1, 2])
    specificity <- conf[1, 1] / (conf[1, 1] + conf[2, 1])
  }

  sens_status <- sens_status[eligible == 1]
  sens_status_predicted <- sens_status_predicted[eligible == 1]

  ## sort the data according to the treatment assignment (treatment first)
  rownames(covar) <- params$size_stage1 + seq(nrow(covar))
  names(covar) <- NULL
  colnames(covar) <- paste("Covar", seq(ncol(covar)), sep = "")

  patients <- data.frame(FID = params$size_stage1 + seq(params$size_stage2), IID = params$size_stage1 + seq(params$size_stage2), treat = treat, sens_status = sens_status, stage = rep(2, params$size_stage2))
  datalist_stage2 <- list(patients = patients, covar = covar, resp_rate = resp_rate, response = response)
  output <- list(datalist_stage2 = datalist_stage2, noneligible = noneligible, sensitivity = sensitivity, 
                 specificity = specificity, sens_status_predicted = sens_status_predicted, cvrs = cvrs)
  return(output)
}

#############################################################
#' @title analyse_fisher
#'
#' @description Perform Fisher's exact test for the sensitive group
#'              
#'
#' @param  datalist
#'         sens_status_predicted
#'
#' @return  pval P-value from Fisher's exact test
#'
#' @author Svetlana Cherlin, James Wason
#############################################################

analyse_fisher <- function(datalist, sens_status_predicted) {
  response <- datalist$response[sens_status_predicted == 1] # response for sensitive group
  treat <- datalist$patients$treat[sens_status_predicted == 1] # treatment allocation for sensitive group
  conf <- matrix(
    nrow = 2, ncol = 2,
    data = c(
      sum(response & !treat), # number of responders in control
      sum(!response & !treat), # number of non-responders in control
      sum(response & treat), # number of responders in treatment
      sum(!response & treat) # number of non-responders in treatment
    ), 
    byrow = TRUE
  )
  pval <- fisher.test(conf, alternative = "two.sided")$p.value
  return(pval)
}

