#' @title Combine (ensemble) casts
#'
#' @description Combine multiple casts' output into a single ensemble. 
#'  Presently, only a general average ensemble is available.
#'
#' @details A pre-loaded table of casts can be input, but if not (default), 
#'  the table will be efficiently (as defined by the inputs) loaded and 
#'  trimmed. \cr 
#'  The casts can be trimmed specifically using the \code{cast_ids} input,
#'  otherwise, all relevant casts from the stated \code{cast_groups}
#'  will be included. 
#'
#' @param main \code{character} value of the name of the main component of
#'  the directory tree.
#'
#' @param method \code{character} value of the name of the ensemble method
#'  to use. Presently, only \code{"unwtavg"} (unweighted average) is allowed.
#'
#' @param cast_groups \code{integer} (or integer \code{numeric}) value
#'  of the cast group to combine with an ensemble. If \code{NULL} (default),
#'  the most recent cast group is ensembled. 
#'
#' @param cast_ids \code{integer} (or integer \code{numeric}) values 
#'  representing the casts of interest for restricting ensembling, as indexed
#'  within the directory in the \code{casts} sub folder. 
#'  See the casts metadata file (\code{casts_metadata.csv}) for summary
#'  information.
#'
#' @param end_moon \code{integer} (or integer \code{numeric}) 
#'  newmoon number of the forecast origin. Default value is 
#'  \code{NULL}, which equates to no selection with respect to 
#'  \code{end_moon}.
#'
#' @param cast_tab Optional \code{data.frame} of cast table outputs. If not
#'  input, will be loaded.
#'
#' @param models \code{character} value(s) of the name of the model to 
#'  include. Default value is \code{NULL}, which equates to no selection with 
#'  respect to \code{model}. \code{NULL} translates to all \code{models}
#'  in the table.
#'
#' @param data_set \code{character} value of the rodent data set to include
#'  Default value is \code{NULL}, which equates to the first data set 
#'  encountered.
#'
#' @param include_interp \code{logical} indicator of if the basic data set
#'  names should also be inclusive of the associated interpolated data sets.
#'
#' @param arg_checks \code{logical} value of if the arguments should be
#'  checked using standard protocols via \code{\link{check_args}}. The 
#'  default (\code{arg_checks = TRUE}) ensures that all inputs are 
#'  formatted correctly and provides directed error messages if not. 
#'
#' @param species \code{character} vector of the species code(s) 
#'  or \code{"total"} for the total across species) to be plotted 
#'  \code{NULL} translates to the species defined by  
#'  \code{evalplot_species}.
#'
#' @return \code{NULL}. Plot is generated.
#' 
#' @examples
#'  \donttest{
#'   setup_dir()
#'   portalcast(models = c("AutoArima", "ESSS"), end_moon = 515)
#'   ensemble_casts()
#' }
#'
#' @export
#'
ensemble_casts <- function(main = ".", method = "unwtavg", 
                           cast_groups = NULL, cast_ids = NULL, 
                           cast_tab = NULL, end_moon = NULL, 
                           models = NULL, data_set = NULL, 
                           include_interp = TRUE,
                           species = NULL, arg_checks = TRUE){

  check_args(arg_checks = arg_checks)
  if(is.null(cast_tab)){
    cast_choices <- select_casts(main = main, cast_ids = cast_ids, 
                                 models = models, end_moons = end_moon, 
                                 data_sets = data_set, 
                                 include_interp = include_interp,
                                 arg_checks = arg_checks)
    if(NROW(cast_choices) == 0){
      stop("no casts available for request", call. = FALSE)
    }else{
      cast_tab <- read_cast_tabs(main = main, cast_ids = cast_choices$cast_id,
                                 arg_checks = arg_checks)
      cast_tab <- add_obs_to_cast_tab(main = main, cast_tab = cast_tab,
                                      arg_checks = arg_checks)
      cast_tab <- add_err_to_cast_tab(main = main, cast_tab = cast_tab,
                                      arg_checks = arg_checks)
      cast_tab <- add_lead_to_cast_tab(main = main, cast_tab = cast_tab,
                                       arg_checks = arg_checks)
      cast_tab <- add_covered_to_cast_tab(main = main, cast_tab = cast_tab,
                                          arg_checks = arg_checks)
    }
  }
  cast_tab$data_set <- gsub("_interp", "", cast_tab$data_set)
  cast_ids <- ifnull(cast_ids, unique(cast_tab$cast_id))
  models <- ifnull(models, unique(cast_tab$model))
  data_set <- ifnull(data_set, "controls")
  species <- ifnull(species, 
                    base_species(total = TRUE, arg_checks = arg_checks)) 
  end_moon <- ifnull(end_moon, max(unique(na.omit(cast_tab)$end_moon))) 
  cast_groups <- ifnull(cast_groups, unique(cast_tab$cast_group))
  cast_id_in <- cast_tab$cast_id %in% cast_ids
  model_in <- cast_tab$model %in% models
  data_set_in <- cast_tab$data_set == data_set
  species_in <- cast_tab$species %in% species
  end_moon_in <- cast_tab$end_moon %in% end_moon
  cast_group_in <- cast_tab$cast_group %in% cast_groups
  all_in <- cast_id_in & model_in & data_set_in & species_in & end_moon_in &
            cast_group_in
  if(sum(all_in) == 0){
    stop("no casts available for requested", call. = FALSE)
  }
  cast_tab <- cast_tab[all_in, ]
  nspecies <- length(species)
  moons <- unique(cast_tab$moon)
  nmoons <- length(moons)
  nmodels <- length(models)
  weight <- 1/nmodels
  estimate <- rep(NA, nspecies * nmoons)
  mvar <- rep(NA, nspecies * nmoons)
  l_pi <- rep(NA, nspecies * nmoons)
  u_pi <- rep(NA, nspecies * nmoons)
  obs <- rep(NA, nspecies * nmoons)
  error <- rep(NA, nspecies * nmoons)
  end_moon_id <- rep(NA, nspecies * nmoons)
  ecast_id <- rep(NA, nspecies * nmoons)
  species_id <- rep(NA, nspecies * nmoons)
  moon_id <- rep(NA, nspecies * nmoons)
  covered <- rep(NA, nspecies * nmoons)
  nmodels <- length(models)
  counter <- 1
  for(i in 1:nspecies){
    for(j in 1:nmoons){
      species_in <- cast_tab$species %in% species[i]
      moon_in <- cast_tab$moon %in% moons[j]
      all_in <- species_in & moon_in
      pcast_tab <- cast_tab[all_in, ]
  
      estimates <- na.omit(pcast_tab$estimate)
      if(length(estimates) > 0){
        mvars <- ((na.omit(pcast_tab$upper_pi) - na.omit(pcast_tab$estimate))/
                  (na.omit(pcast_tab$confidence_level))) ^2
        estimate[counter] <- sum(estimates * weight)
        wt_ss <- sum(weight * (estimates - estimate[counter]) ^ 2)
        if(nmodels > 1){
          mvar[counter] <- sum((mvars * weight) + (wt_ss / (nmodels - 1)))
        }else{
          mvar[counter] <- mvars
        }
        CL <- unique(na.omit(pcast_tab$confidence_level))
        u_pi[counter] <- estimate[counter] + sqrt(mvar[counter] * CL)
        l_pi[counter] <- estimate[counter] - sqrt(mvar[counter] * CL)
        obs[counter] <- unique(pcast_tab$obs)
        error[counter] <- estimate[counter] - obs[counter]
        end_moon_id[counter] <- unique(pcast_tab$end_moon)
        ecast_id[counter] <- as.numeric(paste0(9999, min(pcast_tab$cast_id)))
        moon_id[counter] <- moons[j]
        species_id[counter] <- species[i]
        covered[counter] <- estimate[counter] >= l_pi[counter] &
                            estimate[counter] <= u_pi[counter]
      }
      counter <- counter + 1
      
    }  
  }

  ensemble_name <- paste0("ensemble_", method)

  lead <- moon_id - end_moon_id
  data.frame(model = ensemble_name, moon = moon_id, species = species_id,
             estimate = estimate, var = mvar, lower_pi = l_pi, 
             upper_pi = u_pi, obs = obs, error = error, 
             end_moon = end_moon_id, lead = lead, data_set = data_set,
             cast_id = ecast_id, covered = covered, stringsAsFactors = FALSE) 
}