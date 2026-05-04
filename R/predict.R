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

  trafo_args <- object[["trafo_args"]]
  id_var <- if (!is.null(trafo_args)) trafo_args[["id"]] else NULL
  if (is.null(id_var)) id_var <- object[["attr_ped"]][["id_var"]]
  if (is.null(id_var)) id_var <- attr(newdata, "id_var")
  if (is.null(id_var)) id_var <- "id"

  if (!is.ped(newdata)) {

    brks <- if (!is.null(trafo_args)) trafo_args[["cut"]] else NULL
    if (is.null(brks)) brks <- object[["attr_ped"]][["breaks"]]
    if (is.null(brks)) {
      stop("Can not derive interval cut points from the fitted model.
        Refit with transformation metadata or provide PED newdata.")
    }
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
    arrange(.data[[id_var]], .data$times) %>%
    group_by(.data[[id_var]]) %>%
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
