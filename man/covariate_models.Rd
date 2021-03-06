% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_models.R
\name{covariate_models}
\alias{covariate_models}
\title{Create a covariate model list}
\usage{
covariate_models(model = "pevGARCH", arg_checks = TRUE)
}
\arguments{
\item{model}{\code{character} name for covariate models. Currently only 
\code{"pevGARCH"} is supported.}

\item{arg_checks}{\code{logical} value of if the arguments should be
checked using standard protocols via \code{\link{check_args}}. The 
default (\code{arg_checks = TRUE}) ensures that all inputs are 
formatted correctly and provides directed error messages if not.}
}
\value{
\code{list} of covariate model structures.
}
\description{
Convenience function for creating covariate model 
 \code{list}s.
}
\examples{
 covariate_models()

}
