#' Predict hazard, cumulative hazard or survival probability
#'
#' @param object An object of class \code{pam_xgb}
#' @param newdata A data set containing the same covariates as used for model
#' fitting. If of class \code{data.frame}, the function will try to transform
#' to the \code{xgb.DMatrix} format.
#' @param times A vector of times for which predictions should be generated
#' @param type The type of prediction desired. Either hazard (\code{type = "hazard"}),
#' cumulative hazard (\code{type = "cumu_hazard"}) or survival probability
#' (\code{type = "surv_prob"}).
#' @importFrom stats predict
#' @importFrom purrr map_dfr
#' @return A matrix of predictions containing one row per
#' observation (row in newdata) and 1 column per  specified time in the
#' \code{times} argument.
#' @export
predict.pamm <- function(
  object,
  newdata,
  type = c("hazard", "cumu_hazard", "surv_prob"),
  ...) {

  type <- match.arg(type)

  if (!is.ped(newdata)) {
    newdata <- as_ped(object, newdata)
  }

  newdata[["pred"]] <- predict(unpam(object), newdata, type = "response")
  if (type == "cumu_hazard") {
    newdata <- newdata %>%
      group_by(.data$id) %>%
      mutate(pred = cumsum(.data$pred * exp(.data$offset)))#TODO: is it correct to use offset here?
  }
  if (type == "surv_prob") {
     newdata <- newdata %>%
      group_by(.data[["id"]]) %>%
      mutate(pred = exp(-cumsum(.data$pred * exp(.data$offset))))
  }

  newdata %>%
    group_by(.data[["id"]]) %>%
    filter(row_number() == n()) %>%
    pull(.data[["pred"]]) # TODO: is the hazard/surv prob in the last available interval a useful return?

}


#' S3 method for pamm objects for compatibility with package pec
#'
#' @importFrom pec predictSurvProb
#' @importFrom purrr map
#' @export
predictSurvProb.pamm <- function(
  object,
  newdata,
  times) {

  if (!is.ped(newdata)) {

    trafo_args <- object[["trafo_args"]]
    id_var    <- trafo_args[["id"]]
    brks      <- trafo_args[["cut"]]
    if ( max(times) > max(brks) ) {
      stop("Can not predict beyond the last time point used during model estimation.
        Check the 'times' argument.")
    }
    ped_times <- sort(unique(union(c(0, brks), times)))
    # extract relevant intervals only, keeps data small
    ped_times <- ped_times[ped_times <= max(times)]
    # obtain interval information
    ped_info <- get_intervals(brks, ped_times[-1])
    # add adjusted offset such that cumulative hazard and survival probability
    # can be calculated correctly
    ped_info[["intlen"]] <- c(ped_info[["times"]][1], diff(ped_info[["times"]]))
    # create data set with interval/time + covariate info
    newdata[[id_var]] <- seq_len(nrow(newdata))
    newdata <- combine_df(ped_info, newdata)

  }

  newdata[["pred"]] <- predict(
    unpam(object),
    newdata = newdata,
    type    = "response")
  newdata <- newdata %>%
    arrange(.data$id, .data$times) %>%
    group_by(.data$id) %>%
    mutate(pred = exp(-cumsum(.data$pred * .data$intlen))) %>%
    ungroup() %>%
    filter(.data[["times"]] %in% .env[["times"]])

  id <- unique(newdata[[id_var]])
  pred_list <- map(
    id,
    ~ newdata[newdata[[id_var]] == .x, "pred"] %>% pull("pred"))

  do.call(rbind, pred_list)

}
