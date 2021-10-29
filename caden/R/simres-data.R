#' Results
#'
#' @usage data(simres)
#'
#' @format An object of class \code{"caden"} with the following elements:
#'         datalist.stage2: a list of simulated data for stage 2.
#'         sens.status.predicted: A vector of the predicted sensitivity status.
#'         sensitivity: A vector of values for sensitivity for identifying the sensitive group.
#'         specificity: A vector of values for specificity for identifying the sensitive group.
#'         noneligible: Number of patients non-eligible for the trial (for the enrichment strategy).
#'         pval_overall: P-value for the difference between the arm in the overall trial population.
#'         pval_fisher_group: P-value for the difference between the arms in the sensitive group (as assessed using Fisher exact test).
#'         resp_rate: a vector of response rates
#'         response: a vector of simulated binary responses
#'
#' @keywords datasets
#'
#' @examples
#'
#' data(simdata_stage1)
#' \donttest{
#' simres <- analyse_simdata("param.file", simdata_stage1)
#' }
#'
"simres"
