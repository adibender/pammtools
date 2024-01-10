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

test_df <- data.frame(
  id     = c(1,1, 2,2,2),
  tstart = c(0, .5, 0, .8, 1.2),
  tstop  = c(.5, 3, .8, 1.2, 3),
  status = c(1, 0, 1, 1, 1),
  enum   = c(1, 2, 1, 2, 3),
  age    = c(50, 50, 24, 24, 24))

print(test_df)

# debug multistate
data       <- test_df
formula    <- Surv(tstart, tstop, status)~ enum + age
transition <- "enum"
id         <- "id"
timescale  <- "gap"
min_events <- 1L
censor_code<- 0L
max_time <- 3

# add transition column to formula if not already specified
rhs_vars <- get_rhs_vars(formula)
if (!(transition %in% rhs_vars)) {
  formula <- add_to_rhs(formula, transition)
}

# get dots object
dots <- debug_as_ped_multistate(
  data       = test_df,
  formula    = Surv(tstart, tstop, status)~ enum + age,
  transition = "enum",
  id         = "id",
  timescale  = "gap",
  max_time   = 3)

# # get ped object
# ped <- debug_as_ped_multistate(
#   data       = test_df,
#   formula    = Surv(tstart, tstop, status)~ enum + age,
#   transition = "enum",
#   id         = "id",
#   timescale  = "gap",
#   max_time   = 3)

# call split function
ped <- do.call(split_data_multistate, dots)
# error: splits not at correct time
# -> debug split_data_multistate

# ----------------------------------------------------------------------------#
# copied code from split-data.R
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------------------#
# line 185 onwards
# ----------------------------------------------------------------------------#
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
# works correctly

# ----------------------------------------------------------------------------#
# line 197 onwards
# ----------------------------------------------------------------------------#
## obtain interval breaks points for each spell
#comment out if-clause as timescale == gap
#if(timescale == "gap") {
  data <- mutate(data, !!!list(.time = quo(!!as.name(surv_vars[2]) - !!as.name(surv_vars[1]))))
  formula <- update(formula, Surv(.time, status)~.)
  formula <- update_formula(formula, proposed_names = c(".time", surv_vars[3]))
  # JP 09.01.2024: What does the second update_formula do? -> try  to run without it
#}
# split data for each spell
data_list <- split(data, data[[transition]])
# rm(data)
# only keep spells with minimum number of events
data_list <- data_list[map_dbl(data_list, ~sum(.x[[surv_vars[3]]])) >= min_events]
cuts <-  do.call(debug_split_data_multistate, dots)
# error message: no applicable method for 'get_cut' applied to an object of class "list"
# JP 10.01.2024: changed get_cut -> get_cut.list

lhs_vars <- get_lhs_vars(formula)
if (length(lhs_vars) == 3 & timescale == "gap") {
  rhs_vars <- get_rhs_vars(formula)
  formula_cuts <- as.formula(
    paste0("Surv(", lhs_vars[2], ",", lhs_vars[3], ") ~ ",
           paste(rhs_vars, collapse = "+")))
} else {
  formula_cuts <- formula
}

cuts <- debug_split_data_multistate(
          formula = dots$formula,
          data = dots$data,
          transition = dots$transition,
          cut = dots$cut,
          max_time = dots$max_time,
          event = dots$event,
          min_events = dots$min_events,
          timescale = dots$timescale
)


split_df_list <- map2(
  .x = data_list,
  .y = ifelse(is.list(cuts), cuts, list(cuts)),
  .f = ~ {
    dots$data    <- .x
    dots$formula <- formula
    dots$cut     <- .y
    split_df     <- do.call(split_data, dots)
  }
)

split_df <- bind_rows(split_df_list)
split_df <- split_df %>%
  arrange(.data[[transition]], .data[[dots$id]], .data[["tstart"]])

# remove all obs beyond last observed event time
if (is.null(max_time)) {
  max_time <- max(split_df[["tend"]][split_df[["ped_status"]] == 1])
  split_df <- split_df %>% filter(.data[["tend"]] <= max_time)
}

# if (timescale == "calendar") {
#   split_check <- split_df %>%
#     group_by(.data[[dots$id]]) %>%
#     summarize(dups = sum(duplicated(.data[["tstart"]])))

#   if (any(split_check[["dups"]]) != 0) {
#     stop("Something went wrong during data transformation. \n Please post an issue at 'https://github.com/adibender/pammtools/issues' with your code and data")
#   }
# }

## set class and and attributes
class(split_df) <- c("ped", class(split_df))
attr(split_df, "breaks") <- cuts
attr(split_df, "id_var") <- dots_in$id <- id_var
attr(split_df, "intvars") <- c(id_var, "tstart", "tend", "interval", "offset",
                               "ped_status")
dots_in$transition           <- transition
dots_in$timescale            <- timescale
dots_in$cut                  <- sort(unique(cuts))
dots_in$max_time             <- max_time
dots_in$event                <- event
dots_in$min_events           <- min_events
attr(split_df, "trafo_args") <- dots_in
class(split_df)              <- unique(class(split_df))

split_df


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

