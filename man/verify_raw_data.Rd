% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download.R
\name{verify_raw_data}
\alias{verify_raw_data}
\title{Verify that the raw data folder exists}
\usage{
verify_raw_data(raw_path_data = "PortalData", main = ".")
}
\arguments{
\item{raw_path_data}{\code{character} value indicating the folder path
to the data within the \code{raw} subdirectory but above the files where
the raw data \emph{should} exist. For example, a standard portalcasting 
directory downloads the raw data files into \code{"raw\PortalData"}, 
so \code{raw_location_data = "PortalData"} (as \code{"raw/"} is implied).}

\item{main}{\code{character} value of the name of the main component of
the directory tree.}
}
\value{
\code{logical} indicator of whether or not the raw data folder
 exists.
}
\description{
Check that the raw data folder exists.
}
\examples{
 verify_raw_data()

}