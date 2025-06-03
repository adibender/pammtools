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

#' S3 method for pamm objects for compatibility with package pec
#' 
#' This function is needed to use \code{pec::pec} in the competing risks setting.
#'
#' @inheritParams pec::predictEventProb
#' @importFrom pec predictEventProb
#' @importFrom purrr map
#' @importFrom pammtools get_intervals
#' @export
#' @rdname predictEventProb
predictEventProb.pamm <- function(
    object,
    newdata,
    times,
    cause,
    ...
) {
  # Make a copy of original newdata
  newdata_cs <- newdata
  # If is not PED, transform newdata
  if (!is.ped(newdata_cs)) {
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
    newdata_cs[[id_var]] <- seq_len(nrow(newdata_cs))
    newdata_cs <- combine_df(ped_info, newdata_cs)
  }
  # Recovers all causes
  all_causes <- object[["attr_ped"]]$risks
  # Predict values for each event
  for (cs in all_causes) {
    newdata_cs[["cause"]] <- cs
    newdata_cs[[paste0("csh", cs)]] <- predict(object, 
                                               newdata = newdata_cs, type = "response")
  }  
  # Vector to rename cause-specific hazards
  cause_vars <- paste0("csh", all_causes)
  # Calculate CIF
  newdata_cs <- newdata_cs %>%
    arrange(id, times) %>%
    group_by(id) %>%
    mutate(
      sp_all_cause = exp(-Reduce(
        `+`,
        lapply(cause_vars, function(var) cumsum(get(var) * intlen))
      ))
    ) %>%
    mutate(across(
      all_of(cause_vars),
      ~ cumsum(.x * (sp_all_cause - 1e-5) * intlen),
      .names = "cif{.col}"
    )) %>%
    rename_with(~ sub("cifcsh", "cif", .x), starts_with("cifcsh")) %>% 
    ungroup() %>%
    filter(.data[["times"]] %in% .env[["times"]])
  # Generate IDs
  id <- unique(newdata_cs[[id_var]])
  # Filter CIF for cause of interest
  pred_list <- map(id, ~newdata_cs[newdata_cs[[id_var]] == .x,
                                paste0("cif", cause)] %>%
                     pull(paste0("cif", cause)))
  # Alloc all individuals by all times in one dataframe
  do.call(rbind, pred_list)
}

#' S3 method for compatibility with package pec
#'
#' This function is needed to use \code{pec::pec} in the competing risks setting.
#'
#' @importFrom pec predictEventProb
#' @importFrom purrr map
#' @importFrom pammtools get_intervals
#' @export
#' @rdname predictEventProb
predictEventProb.list <- function(
    object=pam_csh,
    newdata=test_fourD,
    times=times_eval,
    cause,
    ...
) {
  # Make a copy of original newdata
  newdata_cs <- newdata
  # If is not PED, transform newdata
  if (!is.ped(newdata_cs)) {
    trafo_args <- object[[1]][["trafo_args"]]
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
    newdata_cs[[id_var]] <- seq_len(nrow(newdata_cs))
    newdata_cs <- combine_df(ped_info, newdata_cs)
    
  }
  # Recovers all causes
  #all_causes <- c(1:length(object))
  all_causes <- names(object)
  # Predict values for each event
  for (cs in all_causes) {
    name_cause <- sub("cause = ", "", cs)
    newdata_cs[[paste0("csh", name_cause)]] <- predict(object[[cs]], 
                                                       newdata = newdata_cs, 
                                                       type = "response")
  }
  # Vector to rename cause-specific hazards
  cause_vars <- paste0("csh", sub("cause = ", "", all_causes))
  # Calculate CIF
  newdata_cs <- newdata_cs %>%
    arrange(id, times) %>%
    group_by(id) %>%
    mutate(
      sp_all_cause = exp(-Reduce(
        `+`,
        lapply(cause_vars, function(var) cumsum(get(var) * intlen))
      ))
    ) %>%
    mutate(across(
      all_of(cause_vars),
      ~ cumsum(.x * (sp_all_cause - 1e-5) * intlen),
      .names = "cif{.col}"
    )) %>%
    rename_with(~ sub("cifcsh", "cif", .x), starts_with("cifcsh")) %>% 
    ungroup() %>%
    filter(.data[["times"]] %in% .env[["times"]])
  # Generate IDs
  id <- unique(newdata_cs[[id_var]])
  # Filter CIF for cause of interest
  pred_list <- map(id, ~newdata_cs[newdata_cs[[id_var]] == .x,
                                paste0("cif", cause)] %>%
                     pull(paste0("cif", cause)))
  # Alloc all individuals by all times in one dataframe
  do.call(rbind, pred_list)
  
}