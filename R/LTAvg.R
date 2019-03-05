#' @title Long Term Average model for Portal Predictions
#'
#' @description Fit an LTAvg model in the portalcasting pipeline. 
#'
#' @details Model "LTAvg" is a simple long-term average.
#'
#'  Because other models used currently cannot handle missing data, 
#'  we interpolate missing data here for comparison and ensemble building.
#'
#' @param abundances Class-\code{rodents} \code{data.frame} table of rodent 
#'   abundances and time measures.
#'
#' @param metadata Class-\code{metadata} model metadata \code{list}.
#'
#' @param level \code{character} value name of the type of plots included 
#'   (\code{"All"} or \code{"Controls"}).
#'
#' @param quiet \code{logical} value indicating if the function should be 
#'   quiet.
#'
#' @return \code{list} of [1] \code{"forecast"} (the forecasted abundances)
#'   and [2] \code{"all_model_aic"} (the model AIC values).
#'
#' @export
#'
LTAvg <- function(abundances, metadata, level = "All", quiet = FALSE){
  if (!("rodents" %in% class(abundances))){
    stop("`abundances` is not of class rodents")
  }
  if (!("logical" %in% class(quiet))){
    stop("`quiet` is not of class logical")
  }
  if (length(level) > 1){
    stop("`level` can only be of length = 1")
  }
  if (!is.character(level)){
    stop("`level` is not a character")
  }
  if (!any(c("All", "Controls") %in% level)){
    stop("`level` is not valid option")
  } 
  if (!("metadata" %in% class(metadata))){
    stop("`metadata` is not a metadata list")
  } 
  nfcnm <- length(metadata$rodent_forecast_newmoons)
  CL <- metadata$confidence_level
  abundances <- interpolate_abundance(abundances)
  species <- colnames(abundances)[-which(colnames(abundances) == "moons")]
  fcast <- data.frame()
  aic <- data.frame()

  for (s in species){

    ss <- gsub("NA.", "NA", s)
    if (!quiet){
      message(paste0("Fitting LTAvg model for ", ss))
    }
    abund_s <- extract2(abundances, s)

    if (sum(abund_s) == 0){
      model_fcast <- fcast0(nfcnm, "mean")
      model_aic <- 1e6
      estimate <- as.numeric(model_fcast$mean)
      LowerPI <- rep(0, nfcnm)
      UpperPI <- rep(0, nfcnm)

    } else{
      model <- lm(abund_s ~ 1)
      model_fcast <- suppressWarnings( 
                        predict(model, interval = "prediction", level = CL)
                     )

      model_aic <- AIC(model)
      estimate <- as.numeric(model_fcast[1:nfcnm, "fit"])


      LowerPI <- as.numeric(model_fcast[1:nfcnm, "lwr"])
      UpperPI <- as.numeric(model_fcast[1:nfcnm, "upr"])
    }    

    fcast_s <- data.frame(date = metadata$forecast_date, 
                 forecastmonth = metadata$rodent_forecast_months,
                 forecastyear = metadata$rodent_forecast_years, 
                 newmoonnumber = metadata$rodent_forecast_newmoons,
                 currency = "abundance", model = "LTAvg", level = level, 
                 species = ss, estimate = estimate, LowerPI = LowerPI, 
                 UpperPI = UpperPI, fit_start_newmoon = min(abundances$moons),
                 fit_end_newmoon = max(abundances$moons),
                 initial_newmoon = max(abundances$moons),
                 stringsAsFactors = FALSE)

    aic_s <- data.frame(date = metadata$forecast_date, currency = "abundance", 
               model = "LTAvg", level = level, species = ss, 
               aic = model_aic, fit_start_newmoon = min(abundances$moons),
               fit_end_newmoon = max(abundances$moons),
               initial_newmoon = max(abundances$moons),
               stringsAsFactors = FALSE)

    fcast <- rbind(fcast, fcast_s)
    aic <- rbind(aic, aic_s)
  }
  list("forecast" = fcast, "aic" = aic)
}
