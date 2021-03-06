% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_rodents.R
\name{prefab_data_sets}
\alias{prefab_data_sets}
\title{Provide the names of the prefab data sets}
\usage{
prefab_data_sets(interpolate = NULL)
}
\arguments{
\item{interpolate}{\code{logical} value indicating if the interpolated data
only should be listed (\code{interpolate = TRUE}), if the non-interpolated
data only should be listed (\code{interpolate = FALSE}), or if both
should be listed (\code{interpolate = NULL}, default).}
}
\value{
\code{character} vector of data set names.
}
\description{
Create a \code{character} vector of the names of the rodents
 data sets to be included including the pre-fabricated (prefab) data sets,
 extracting them from the \code{prefab_rodents_controls} internal 
 function. \cr \cr
 The currently prefabricated options include (\code{"all"}, 
 \code{"controls"}, \code{"all_interp"}, and \code{"controls_interp"}).
}
\examples{
 prefab_data_sets()

}
