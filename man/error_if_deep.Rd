% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{error_if_deep}
\alias{error_if_deep}
\title{Error if a function's request is deeper than can be handled}
\usage{
error_if_deep(lev)
}
\arguments{
\item{lev}{The number of frames back in the stack where the request needs
to be able to be evaluated.}
}
\value{
Throws an error if the function is called in a place where it 
 cannot operate and returns \code{NULL} otherwise.
}
\description{
Produces an informative error message when a function 
 that should only be called inside of other functions is called outside
 of a function (hence the request to the function is too deep for
 what it can handle).
}
\examples{
 \dontrun{
 # will error:
 # error_if_deep(-10)
 }
 error_if_deep(0)

}
