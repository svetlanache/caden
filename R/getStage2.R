#' @title
#' Simulate a data set for stage 2 for the the cross-validate adaptive enrichment risk scores (CADEN) design
#'
#' @description Create a data set for stage 2 by sampling with replacement from the stage 1 data set
#'
#' @param datalist1 A list with two data frames (patients, covar) and avector of binary responses
#'         patients: a data frame with one row per patient and the following columns:
#'         FID (family ID), IID (individual ID), treat (1 for treatment and 0 for control)
#'         covar: a data frame with covariate data
#'         response: a vector of simulated binary responses
#' @param seed
#'
#'
#' @return datalist2 A list with two data frames (patients, covar) and avector of binary responses
#'         patients: a data frame with one row per patient and the following columns:
#'         FID (family ID), IID (individual ID), treat (1 for treatment and 0 for control)
#'         covar: a data frame with covariate data
#'         response: a vector of simulated binary responses
#'
#' @examples
#' data(realdata_stage1)
#' data(realdata_stage2)
#' real <- analyseReal(realdata_stage1, realdata_stage2, 0.04, 0.01, 123)
#' @seealso
#' \code{\link{analyse_realdata}}, function.
#' @author Svetlana Cherlin, James Wason
#' @export get_stage2
#'

get_stage2 <- function(datalist1, seed) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  ## sample with replacement

  nc <- nrow(datalist1$patients[datalist1$patients$treat == 0, ])
  nt <- nrow(datalist1$patients[datalist1$patients$treat == 1, ])

  wc <- sample(seq(nc), nc, replace = TRUE)
  wt <- sample(seq(nt), nt, replace = TRUE)

  idc <- as.character(datalist1$patients$IID[datalist1$patients$treat == 0])
  idt <- as.character(datalist1$patients$IID[datalist1$patients$treat == 1])

  id2c <- idc[wc]
  id2t <- idt[wt]
  id2 <- c(id2c, id2t)

  m <- match(id2, as.character(datalist1$patients$IID))

  patients2 <- datalist1$patients[m, ]
  covar2 <- datalist1$covar[m, ]
  response2 <- datalist1$response[m]
  datalist2 <- list(patients = patients2, covar = covar2, response = response2)

  return(datalist2)
}
