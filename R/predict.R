#' S3 method for pamm objects for compatibility with package pec
#'
#' @inheritParams pec::predictSurvProb
#' @importFrom pec predictSurvProb
#' @importFrom purrr map
#'
#' @export
predictSurvProb.pamm <- function(
  object,
  newdata,
  times,
  ...) {

  if (!is.ped(newdata)) {

    trafo_args <- object[["trafo_args"]]
    id_var     <- trafo_args[["id"]]
    brks       <- trafo_args[["cut"]]
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

  env_times <- times
  newdata[["pred"]] <- unname(predict(
    unpam(object),
    newdata = newdata,
    type    = "response"))
  newdata <- newdata %>%
    arrange(.data$id, .data$times) %>%
    group_by(.data$id) %>%
    mutate(pred = exp(-cumsum(.data$pred * .data$intlen))) %>%
    ungroup() %>%
    filter(.data[["times"]] %in% env_times)

  id <- unique(newdata[[id_var]])
  pred_list <- map(
      id,
      ~ newdata[newdata[[id_var]] == .x, "pred"] %>%
    pull("pred"))

  do.call(rbind, pred_list)

}
