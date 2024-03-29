% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{get_subgroup}
\alias{get_subgroup}
\title{analyse_subgroup}
\usage{
get_subgroup(datalist, seed, standardise_cvrs, full_model)
}
\arguments{
\item{datalist}{seed 
standardise_cvrs
full_model}
}
\value{
pval_treat_group P-value for the treatment effect in the sensitive group

A list with the following elements
         cvrs Risk scores
         sens_status_predicted
         sensitivity
         specificity
}
\description{
Compute the treatment effect in the sensitive group

Compute the subgroup memebership, sensitivity, 
             specificity and risk scores.
}
\author{
Svetlana Cherlin, James Wason
}
