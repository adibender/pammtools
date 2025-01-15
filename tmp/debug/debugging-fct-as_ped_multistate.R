debug_as_ped_multistate <- function(
    data,
    formula,
    cut          = NULL,
    max_time     = NULL,
    tdc_specials = c("concurrent", "cumulative"),
    censor_code  = 0L,
    transition  = character(),
    timescale    = c("gap", "calendar"),
    min_events   = 1L,
    ...
) {
  
  assert_character(transition, min.chars = 1L, min.len = 1L, any.missing = FALSE,
                   len = 1L)
  assert_integer(min_events, lower = 1L, len = 1L)
  
  status_error(data, formula, censor_code)
  assert_subset(tdc_specials, c("concurrent", "cumulative"))
  
  rhs_vars <- get_rhs_vars(formula)
  if (!(transition %in% rhs_vars)) {
    formula <- add_to_rhs(formula, transition)
  }
  
  dots            <- list(...)
  dots$data       <- data
  dots$formula    <- get_ped_form(formula, data = data, tdc_specials = tdc_specials)
  dots$cut        <- sort(unique(cut))
  dots$max_time   <- max_time
  dots$transition <- transition
  dots$min_events <- min_events
  dots$timescale  <- timescale
  
  # ped <- do.call(split_data_multistate, dots)
  #attr(ped, "time_var")   <- get_lhs_vars(dots$formula)[1]
  
  return(dots)
  
}


debug_split_data_multistate <- function(
    formula,
    data,
    transition    = character(),
    cut        = NULL,
    max_time   = NULL,
    event      = 1L,
    min_events = 1L,
    timescale = c("gap", "calendar"),
    ...) {
  
  assert_character(transition, min.chars = 1L, min.len = 1L, any.missing = FALSE,
                   len = 1L)
  assert_integer(min_events, lower = 1L, len = 1L)
  assert_character(timescale)
  timescale <- match.arg(timescale)
  
  dots_in <- list(...)
  dots_in$formula <- formula
  
  ## assert that inputs have correct formats
  assert_class(formula, "formula")
  assert_data_frame(data, min.rows = 1, min.cols = 2)
  assert_numeric(cut, lower = 0, finite = TRUE, any.missing = FALSE,
                 min.len = 1, null.ok = TRUE)
  assert_number(max_time, lower = 0, finite = TRUE, null.ok = TRUE)
  
  
  ## extract names for event time and status variables
  surv_vars <- all.vars(update(formula, .~0))
  vars <- if ("." %in% all.vars(formula)) {
    names(data)
  } else {
    all.vars(formula)
  }
  uvars <- union(surv_vars, vars)
  if (!all(uvars %in% vars)) {
    stop(paste("Variables provided in formula not in data set:",
               paste0(setdiff(uvars, vars), collapse = ", ")))
  }
  
  ## obtain interval breaks points for each spell
  if(timescale == "gap") {
    data <- mutate(data, !!!list(.time = quo(!!as.name(surv_vars[2]) - !!as.name(surv_vars[1]))))
    formula <- update(formula, Surv(.time, status)~.)
    formula <- update_formula(formula, proposed_names = c(".time", surv_vars[3]))
  }
  # split data for each spell
  data_list <- split(data, data[[transition]])
  # rm(data)
  # only keep spells with minimum number of events
  data_list <- data_list[map_dbl(data_list, ~sum(.x[[surv_vars[3]]])) >= min_events]
  cuts <- get_cut.list(data_list, formula, cut = cut, max_time = max_time,
                  event = event, timescale = timescale)

  
  ## create argument list to be passed to split_data
  dots <- list(...)
  dots$multiple_id <- TRUE # possible in case of multi-state models with back transitions

  # if id allready in the data set, remove id variable from dots but keep
  # id variable for later rearrangment
  if (!is.null(dots$id)) {
    id_var <- dots$id
  } else {
    id_var  <- "id"
    dots$id <- id_var
  }
  
  return(cuts)
  
  # 
  # split_df_list <- map2(
  #   .x = data_list,
  #   .y = ifelse(is.list(cuts), cuts, list(cuts)),
  #   .f = ~ {
  #     dots$data    <- .x
  #     dots$formula <- formula
  #     dots$cut     <- .y
  #     split_df     <- do.call(split_data, dots)
  #   }
  # )
  # 
  # split_df <- bind_rows(split_df_list)
  # split_df <- split_df %>%
  #   arrange(.data[[transition]], .data[[dots$id]], .data[["tstart"]])
  # 
  # # remove all obs beyond last observed event time
  # if (is.null(max_time)) {
  #   max_time <- max(split_df[["tend"]][split_df[["ped_status"]] == 1])
  #   split_df <- split_df %>% filter(.data[["tend"]] <= max_time)
  # }
  # 
  # # if (timescale == "calendar") {
  # #   split_check <- split_df %>%
  # #     group_by(.data[[dots$id]]) %>%
  # #     summarize(dups = sum(duplicated(.data[["tstart"]])))
  # 
  # #   if (any(split_check[["dups"]]) != 0) {
  # #     stop("Something went wrong during data transformation. \n Please post an issue at 'https://github.com/adibender/pammtools/issues' with your code and data")
  # #   }
  # # }
  # 
  # ## set class and and attributes
  # class(split_df) <- c("ped", class(split_df))
  # attr(split_df, "breaks") <- cuts
  # attr(split_df, "id_var") <- dots_in$id <- id_var
  # attr(split_df, "intvars") <- c(id_var, "tstart", "tend", "interval", "offset",
  #                                "ped_status")
  # dots_in$transition           <- transition
  # dots_in$timescale            <- timescale
  # dots_in$cut                  <- sort(unique(cuts))
  # dots_in$max_time             <- max_time
  # dots_in$event                <- event
  # dots_in$min_events           <- min_events
  # attr(split_df, "trafo_args") <- dots_in
  # class(split_df)              <- unique(class(split_df))
  # 
  # split_df
  
}

debug_get_cut.list <- function (
    data,
    formula,
    cut       = NULL,
    max_time  = NULL,
    event     = 1L,
    timescale = "gap",
    ...) {
  
  lhs_vars <- get_lhs_vars(formula)
  if (length(lhs_vars) == 3 & timescale == "gap") {
    rhs_vars <- get_rhs_vars(formula)
    formula_cuts <- as.formula(
      paste0("Surv(", lhs_vars[2], ",", lhs_vars[3], ") ~ ",
             paste(rhs_vars, collapse = "+")))
  } else {
    formula_cuts <- formula
  }
  cuts <- map(
    .x = data,
    .f = ~get_cut.default(
      data     = .x,
      formula  = formula_cuts,
      cut      = cut,
      max_time = max_time,
      event    = event,
      ...)
  )
  
  cuts <- Reduce(union, cuts)
  sort(unique(cuts))
  
}

helper_as_ped_multistate <- function(
    data,
    formula,
    cut          = NULL,
    max_time     = NULL,
    tdc_specials = c("concurrent", "cumulative"),
    censor_code  = 0L,
    transition  = character(),
    timescale    = c("gap", "calendar"),
    min_events   = 1L,
    ...
) {
  
  assert_character(transition, min.chars = 1L, min.len = 1L, any.missing = FALSE,
                   len = 1L)
  assert_integer(min_events, lower = 1L, len = 1L)
  
  status_error(data, formula, censor_code)
  assert_subset(tdc_specials, c("concurrent", "cumulative"))
  
  rhs_vars <- get_rhs_vars(formula)
  if (!(transition %in% rhs_vars)) {
    formula <- add_to_rhs(formula, transition)
  }
  if (!(transition %in% rhs_vars)) {
    formula <- add_to_rhs(formula, transition)
  }
  
  # check if transformation contains state or trafo
  if (is.character(transition)) {
    status <- as.name(get_lhs_vars(formula)[3])
    msg_if <- "If-clause works"
  }
  
  dots            <- list(...)
  dots$data       <- data
  dots$formula    <- get_ped_form(formula, data = data, tdc_specials = tdc_specials)
  dots$cut        <- sort(unique(cut))
  dots$max_time   <- max_time
  dots$transition <- transition
  dots$min_events <- min_events
  dots$timescale  <- timescale
  
  # ped <- do.call(split_data_multistate, dots)
  #attr(ped, "time_var")   <- get_lhs_vars(dots$formula)[1]
  
  return(msg_if)
  
}
test_df <- data.frame(
  id     = c(1,1,1, 2,2,2),
  tstart = c(0, .5, 1, 0, 0, 2),
  tstop  = c(.5, 1, 1.5, 1.5, 2.5, 2.5),
  status = c(1, 1, 1, 1, 0, 1),
  state   = c(1, 2, 3, 1, 2, 3),
  age    = c(50,50,50, 24,24,24))

print(test_df)

ped_test_df <- helper_as_ped_multistate(
  data = test_df,
  formula = Surv(tstart, tstop, status) ~ state + age, 
  transition = "state",
  timescale = "calendar",
  id = "id"
)


test_df2 <- test_df %>%
  filter(status == 1) %>%
  group_by(id, status) %>%
  mutate(
    state_start = coalesce(lag(state), state),
    transition = paste0(state_start, "->", state)) %>%
  select(id, tstart, tstop, state_start, state, transition, status, age)
print(test_df2)

ped_test_df2 <- as_ped_multistate(
  data = test_df2,
  formula = Surv(tstart, tstop, status) ~ state + age, 
  transition = "transition",
  timescale = "calendar",
  id = "id"
)

print(ped_test_df %>%filter(id==2))
