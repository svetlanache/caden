% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{get_model}
\alias{get_model}
\title{get_model}
\usage{
get_model(datalist, standardise_cvrs, full_model)
}
\arguments{
\item{datalist_stage1}{threshold_overall
threshold_group
standardise_cvrs
full_model}
}
\value{
A list with the following elements
         decision
         sens_status_predicted
         model
}
\description{
Compute the coefficients for the risk scores.
}
\author{
Svetlana Cherlin, James Wason
}
