% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_models.R
\name{update_models}
\alias{update_models}
\title{Update models based on user input controls}
\usage{
update_models(
  main = ".",
  models = NULL,
  controls_model = NULL,
  update_prefab_models = FALSE,
  control_files = files_control(),
  bline = FALSE,
  quiet = FALSE,
  verbose = FALSE,
  arg_checks = TRUE
)
}
\arguments{
\item{main}{\code{character} value of the name of the main component of
the directory tree.}

\item{models}{\code{character} vector of name(s) of model(s), used to 
restrict the models to update from \code{controls_model}. The default
(\code{NULL}) translates to all models in the \code{controls_model} list,
which is recommended for general usage and is enforced within the main
\code{\link{portalcast}} pipeline.}

\item{controls_model}{Controls for models not in the prefab set or for 
overriding those in the prefab set. \cr 
A \code{list} of a single model's script-writing controls or a
\code{list} of \code{list}s, each of which is a single model's 
script-writing controls. \cr 
Presently, each model's script writing controls should include four 
elements: 
\itemize{
 \item \code{name}: a \code{character} value of the model name.
 \item \code{data_sets}: a \code{character} vector of the data set names
  that the model is applied to. 
 \item \code{covariatesTF}: a \code{logical} indicator of if the 
  model needs covariates.
 \item \code{lag}: an \code{integer}-conformable value of the lag to use 
  with the covariates or \code{NA} if \code{covariatesTF = FALSE}.
} 
If only a single model is added, the name of 
the model from the element \code{name} will be used to name the model's
\code{list} in the larger \code{list}. If multiple models are added, each
element \code{list} must be named according to the model and the
\code{name} element. \cr}

\item{update_prefab_models}{\code{logical} indicator if all of the models'
scripts should be updated, even if they do not have an explicit change
to their model options via \code{controls_model}. Default is
\code{FALSE}, which leads to only the models in \code{controls_model}
having their scripts re-written. Switching to \code{TRUE} results in the
models listed in \code{models} having their scripts re-written. \cr \cr
This is particularly helpful when one is changing the global (with respect
to the models) options \code{main}, \code{quiet}, \code{verbose}, or
\code{control_files}.}

\item{control_files}{\code{list} of names of the folders and files within
the sub directories and saving strategies (save, overwrite, append, etc.).
Generally shouldn't need to be edited. See \code{\link{files_control}}.}

\item{bline}{\code{logical} indicator if horizontal break lines should be
included in messaging.}

\item{quiet}{\code{logical} indicator if progress messages should be
quieted.}

\item{verbose}{\code{logical} indicator of whether or not to print out
all of the information or not (and thus just the tidy messages).}

\item{arg_checks}{\code{logical} value of if the arguments should be
checked using standard protocols via \code{\link{check_args}}. The 
default (\code{arg_checks = TRUE}) ensures that all inputs are 
formatted correctly and provides directed error messages if not.}
}
\value{
\code{NULL}.
}
\description{
Update model scripts based on the user-defined model control
 inputs. This allows users to define their own models or re-define prefab
 models within the \code{\link{portalcast}} pipeline.
}
\examples{
 \donttest{
  setup_dir()
  cm <- model_control(name = "AutoArima", data_sets = "all")
  update_models(controls_model = cm)
 }

}
