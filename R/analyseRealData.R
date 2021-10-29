#' @title
#' Analysis of the the real data (retrospective) according to the cross-validate adaptive enrichment risk scores (CADEN) design
#'
#' @description The design has two stages. At the end of stage 1, an interim analysis is performed to test the efficacy of the treatment compared to the control in the overall trial population. If the test in the overall population is significanct, the trial proceeds into stage 2 by analysing all patients (the "unselected" strategy). If the test in the overall population is not significant then the trial proceeds by testing whether there is a subgroup of patients (sensitive group) that show a promising treatment effect. Depending on the results of the test in the sensitive group, the design proceeds into stage 2 according to one of the mutually exclusive strategies "stop" or "enrichment".
#'
#' @param realdata_stage1 A list with two data frames (patients, covar) and two vectors (resp.rate, response)
#'         patients: a data frame with one row per patient and the following columns:
#'         FID (family ID), IID (individual ID), treat (1 for treatment and 0 for control),
#'         sens_status (true sensitivity status), stage (1)
#'         covar: a data frame with covariate data
#'         response: a vector of simulated binary responses
#' @param realdata_stage2 A list with two data frames with the same structure as realdata_stage1
#' @param threshold_overall P-value threshold for the test for the differences in the treatment effect in the overall trial population
#' @param threshold_group P-value threshold for the test for the treatment effect in the sensitive group
#' @param seed A seed for random number generating
#' @param standardise_cvrs A logical flag for the standardisation of the risk scores. Default is 'standardise_cvrs=TRUE'.
#'        The standardisation is performed with respect to the training data sets, per cross-validation fold.
#' @param full_model A logical flag for the full model (treatment effect, covariate effect and the interaction effect). Default is 1full_model = TRUE. When 'full_model = FALSE', only interaction effect is included in the model
#'
#'
#' @return An object of class \code{"caden"}.
#' @return decision: Enrichment, stop or unselected.
#' @return cvrs: A vector of the risk scores
#' @return sens_status_predicted: A vector of the predicted sensitivity status.
#' @return noneligible: Number of patients non-eligible for the trial (for the enrichment strategy).
#' @return pval_overall: P-value for the difference between the arm in the overall trial population.
#' @return pval_fisher_group: P-value for the difference between the arms in the sensitive group (as assessed using Fisher exact test).
#'
#' @examples
#' data(realdata_stage1)
#' data(realdata_stage2)
#' threshold_overall = 0.04
#' threshold_group = 0.1
#' seed = 123
#' standardise_cvrs = 0
#' full_model = 0
#' real_res <- analyse_realdata(realdata_stage1, realdata_stage2, threshold_overall, threshold_group, seed, standardise_cvrs, full_model)
#' @seealso
#' \code{\link{simulate_data}}, function.
#' @author Svetlana Cherlin, James Wason
#' @export analyse_realdata
#'

analyse_realdata <- function(datalist_stage1, datalist_stage2, threshold_overall, 
                             threshold_group, seed, standardise_cvrs = TRUE, full_model = TRUE) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  sens_status_predicted <- NA
  pval_overall <- NA
  pval_fisher_group <- NA

  ## Analyse stage 1
  res1 <- analyse_stage1(datalist_stage1, threshold_overall, threshold_group, seed, standardise_cvrs)
  cvrs <- res1$cvrs
  sens_status_predicted <- res1$sens_status_predicted

  ## Analyse stage 2

  if (res1$decision == "enrichment") {
    sens_status_predicted2 <- vector(length = nrow(datalist_stage2$patients))
    eligible <- vector(length = nrow(datalist_stage2$patients))
    cvrs2 <- vector(length = nrow(datalist_stage2$patients))

    for (i in 1:nrow(datalist_stage2$patients)) {
      res <- check_eligibility(datalist_stage2$covar[i, ], res1$model, standardise_cvrs)
      sens_status_predicted2[i] <- res$sens_status_predicted
      eligible[i] <- res$eligible
      cvrs2[i] <- res$cvrs
    }

    if (sum(eligible == 1) > 0) {
      patients2 <- datalist_stage2$patients[eligible == 1, ]
      covar2 <- datalist_stage2$covar[eligible == 1, ]
      response2 <- datalist_stage2$response[eligible == 1]
      sens_status_predicted <- c(sens_status_predicted, sens_status_predicted2[eligible == 1])
      cvrs <- c(cvrs, cvrs2[eligible == 1])
      ## combine data from the two stages
      datalist <- list(patients = rbind(datalist_stage1$patients, patients2), covar = rbind(datalist_stage1$covar, covar2), response = c(datalist_stage1$response, response2))
    }
  }

  if (res1$decision == "unselected") {
    datalist <- list(patients = rbind(datalist_stage1$patients, datalist_stage2$patients), covar = rbind(datalist_stage1$covar, datalist_stage2$covar), response = c(datalist_stage1$response, datalist_stage2$response))
    out <- get_subgroup(datalist, seed, standardise_cvrs, full_model)
    sens_status_predicted <- out$sens_status_predicted
    size_treat <- nrow(datalist$patients[datalist$patients$treat == 1, ])
    size_con <- nrow(datalist$patients[datalist$patients$treat == 0, ])
    cvrs <- out$cvrs
    pval_overall <- prop.test(
      x = c(sum(datalist$response[datalist$patients$treat == 0]), sum(datalist$response[datalist$patients$treat == 1])),
      n = c(size_con, size_treat), alternative = "two.sided"
    )$p.value
  }

  if (res1$decision != "stop") {
    ## Fisher's exact test
    pval_fisher_group <- analyse_fisher(datalist, sens_status_predicted)
  }

  output <- list(decision = res1$decision, sens_status_predicted = sens_status_predicted, pval_overall = pval_overall, pval_fisher_group = pval_fisher_group, cvrs = cvrs, standardise_cvrs = standardise_cvrs, full_model = full_model)

  class(output) <- "caden"
  return(output)
}


