% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/realdata_stage1-data.R
\docType{data}
\name{realdata_stage1}
\alias{realdata_stage1}
\title{An example of the real data}
\format{
A list with two data frames (patients, covar) and two vectors (resp.rate, response)
        A list with two data frames (patients, covar) and two vectors (resp.rate, response)
        patients: a data frame with one row per patient and the following columns:
        FID (family ID), IID (individual ID), treat (1 for treatment and 0 for control),
        covar: a data frame with covariate data for L covariates
        response: a vector of simulated binary responses
}
\usage{
data(realdata_stage1)
}
\description{
An example of the real data for the cross-validate adaptive enrichment risk scores (CADEN) design
}
\examples{

data(realdata_stage1)
\donttest{
simres <- analyse_realdata(realdata_stage1, realdata_stage2, 0.04, 0.01, 123)
}

}
\keyword{datasets}
