#' Simulated data
#'
#' Simulated data for the cross-validate adaptive enrichment risk scores (CADEN) design
#' @docType data
#'
#' @usage data(simdata_stage1)
#'
#' @format A list with two data frames (patients, covar) and two vectors (resp.rate, response)
#'         A list with two data frames (patients, covar) and two vectors (resp.rate, response)
#'         patients: a data frame with one row per patient and the following columns:
#'         FID (family ID), IID (individual ID), treat (1 for treatment and 0 for control),
#'         sens.status (true sensitivity status), stage (1)
#'         covar: a data frame with covariate data
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
"simdata_stage1"
