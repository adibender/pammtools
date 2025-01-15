setwd("C:/Users/ra63liw/Documents/98_git/pammtools-multi-state/pammtools")


# create test file
test_df <- data.frame(
  id     = c(1,1, 2,2,2),
  tstart = c(0, .5, 0, .8, 1.2),
  tstop  = c(.5, 3, .8, 1.2, 3),
  status = c(1, 0, 1, 1, 2),
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


# ----------------------------------------------------------------------------#
# copied code from as-ped.R
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------------------#
# line 362 onwards
# ----------------------------------------------------------------------------#
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

# # call split function
# ped <- do.call(split_data_multistate, dots)
# # error: splits not at correct time
# # -> debug split_data_multistate

# ----------------------------------------------------------------------------#
# copied code from split-data.R
# ----------------------------------------------------------------------------#

# initialize variables
id <- dots$id
data <- dots$data
formula <- dots$formula
max_time <- dots$max_time
transition <- dots$transition
min_events <- dots$min_events
timescale <- dots$timescale

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

# control:
print(data_list)

# rm(data)
# only keep spells with minimum number of events
data_list <- data_list[map_dbl(data_list, ~sum(.x[[surv_vars[3]]])) >= min_events]
# control:
print(data_list)

cuts <-  do.call(debug_split_data_multistate, dots)
# control:
print(cuts)
# error message: no applicable method for 'get_cut' applied to an object of class "list"
# JP 10.01.2024: changed get_cut -> get_cut.list
# cut point is 1.8 instead of 1.2, hence error lies within split_data_multistate
# within split_data_multistate the error lies in the function get_cut.list


# ----------------------------------------------------------------------------#
# copied code from get-cut-points.R
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------------------#
# debugging get_cut.list
# line 44 onwards
# ----------------------------------------------------------------------------#
lhs_vars <- get_lhs_vars(formula)
# control:
print(lhs_vars)

if (length(lhs_vars) == 3 & timescale == "gap") {
  rhs_vars <- get_rhs_vars(formula)
  formula_cuts <- as.formula(
    paste0("Surv(", lhs_vars[2], ",", lhs_vars[3], ") ~ ",
           paste(rhs_vars, collapse = "+")))
} else {
  formula_cuts <- formula
}
# control:
print(formula_cuts)

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
# error: Caused by error in `data[[outcome_vars[1]]]`:
# ! subscript out of bounds
# run manually

# ----------------------------------------------------------------------------#
# debugging get_cut.default
# line 23 onwards
# ----------------------------------------------------------------------------#

#initialize variables
data <- data
formula <- formula_cuts
cut <- dots$cut
max_time <- dots$max_time
event <- 1L

# comment out if-clause because cut == NULL
# if (is.null(cut)) {
  outcome_vars <- get_lhs_vars(formula)
  if (length(outcome_vars) == 2) {
    cut <- unique(data[[outcome_vars[1]]][1L * (data[[outcome_vars[2]]]) == event])
  } else {
    cut_start <- unique(data[[outcome_vars[1]]])
    cut_end   <- unique(data[[outcome_vars[2]]])
    cut       <- union(cut_start, cut_end)
  }
  if (!is.null(max_time)) {
    cut <- cut[cut < max_time]
    cut <- c(cut, max_time)
  }
# }












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
# control:
print(split_df_list)

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