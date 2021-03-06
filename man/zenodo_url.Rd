% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download.R
\name{zenodo_url}
\alias{zenodo_url}
\alias{zenodo_versions}
\title{Obtain the URL for a Zenodo record to be downloaded}
\usage{
zenodo_url(
  concept_rec_id = NULL,
  rec_version = "latest",
  rec_id = NULL,
  arg_checks = TRUE
)

zenodo_versions(concept_rec_id, arg_checks = TRUE)
}
\arguments{
\item{concept_rec_id}{Concept record identifier, a \code{character} value
corresponding to the Zenodo concept.}

\item{rec_version}{\code{character} value of the version number or 
\code{"latest"} (default) for the data to be download.}

\item{rec_id}{Optional input record identifier, a \code{character} value
corresponding to the Zenodo record.}

\item{arg_checks}{\code{logical} value of if the arguments should be
checked using standard protocols via \code{\link{check_args}}. The 
default (\code{arg_checks = TRUE}) ensures that all inputs are 
formatted correctly and provides directed error messages if not.}
}
\value{
\code{zenodo_url}: \code{character} value of the URL for the zip 
 to be downloaded. \cr \cr
 \code{zenodo_versions}: a \code{data.frame} of version number and record 
 identifier for each version available.
}
\description{
\code{zenodo_url} obtains the URL for a given Zenodo record, 
 identified either by the concept record identifier (\code{concept_rec_id}) 
 and version (\code{rec_version}) or record identifier (\code{rec_id}). 
 (\strong{Note}: if \code{rec_id} is used, it overrides 
 \code{concept_rec_id}). \cr \cr
 \code{zenodo_versions}: determines the available version numbers and the 
 corresponding record identifier for each version available for a given 
 Zenodo concept (group of records).
}
\examples{
 \donttest{
   zenodo_versions("1215988")
   zenodo_url("1215988", "latest")
   zenodo_url("1215988", "1.71.0")
   zenodo_url(rec_id = "1217163")
 }

}
