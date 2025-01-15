add_trans_prob <- function(
    newdata
    , object
    , overwrite       = FALSE 
    , alpha           = 0.05
    , n_sim           = 500L
    , time_var        = NULL
    , interval_length = "intlen",
    ...
  ) {

  
  
  if (!overwrite) {
    if ("trans_prob" %in% names(newdata)) {
      stop("Data set already contains 'trans_prob' column.
        Set `overwrite=TRUE` to overwrite")
    }
  } else {
    rm.vars <- intersect(
      c("trans_prob"
        # , "surv_lower" # not yet implemented
        # , "surv_upper" # not yet implemented
      ),
      names(newdata))
    newdata <- newdata %>% select(-one_of(rm.vars))
  }
  
  # extract to simulate ci
  X             <- predict.gam(object, newdata = newdata, type = "lpmatrix")
  coefs         <- coef(object)
  V             <- object$Vp
  sim_coef_mat  <- mvtnorm::rmvnorm(n_sim, mean = coefs, sigma = V)
  
  # add cumu_hazard to newdata
  newdata <- newdata |> add_cumu_hazard(object, ci=F)
  
  old_groups <- dplyr::groups(newdata)
  res_data <- newdata %>% ungroup(transition)
  out_data <- group_split(res_data) |> 
    map(res_data, .f = ~ group_by(.x, transition))|> 
    map(res_data, .f = ~ get_trans_prob(.x)) |>
    map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
    bind_rows()
  
  return(out_data)
  
}



# ---------------------------------------------------------------------------- #
# build function to calculate the transition probabilities with pammtools, 
# using hazards
# ---------------------------------------------------------------------------- #
get_trans_prob <- function(
    newdata,
    # object,
    time_var   = NULL,
    interval_length = "intlen",
    transition = "transition",
    tend = "tend",
    cumu_hazard = "cumu_hazard",
    ...) {
  
  # interval_length
  assert_character(interval_length)
  assert_subset(interval_length, colnames(newdata))
  # transition
  assert_character(transition)
  assert_subset(transition, colnames(newdata))
  # time
  assert_character(tend)
  assert_subset(tend, colnames(newdata))
  # cumu_hazard
  assert_character(cumu_hazard)
  assert_subset(cumu_hazard, colnames(newdata))
  assert_data_frame(newdata, all.missing = FALSE)
  
  # assert_class(object, classes = "glm")
  
  interval_length <- sym(interval_length)
  transition <- sym(transition)
  
  # include from and to, to obtain transition probability in multidim array
  newdata <- newdata %>% 
    mutate(from = as.numeric(gsub("->.*", "", !!transition))
           , to = as.numeric(gsub(".*->", "", !!transition)))
  
  # get unique transitions to build transition matrix
  unique_transition <- data.frame(unique(newdata %>% select(!!transition, from, to)))
  # get unique time points
  unique_tend <- data.frame(unique(newdata %>% 
                                     ungroup(!!transition) %>% 
                                     select(!!tend)))
  
  # transition matrix
  m <- sapply(unique_transition[,c(2,3)], max) + 1 #transition starts at 0, integer of matrix at 1
  M <- array(0, dim=c(max(m), max(m), nrow(unique_transition))) 
  
  
  # create transition matrices to be used at every time point,
  # multiply matrices with "scalar" alpha_ij_k which is the delta cumu hazard at time t_k for transition i->j
  
  for (iter in 1:nrow(unique_transition)){
    M[unique_transition$from[iter] + 1, unique_transition$to[iter] + 1,iter] <- 1
    M[unique_transition$from[iter] + 1, unique_transition$from[iter] + 1,iter] <- -1
  }
  
  # add cumu hazards to dataset
  newdata <- newdata %>% 
    # group_by(!!transition) %>%
    mutate(delta_cumu_hazard = cumu_hazard - ifelse(is.na(lag(cumu_hazard)), 0, lag(cumu_hazard)))
  
  # create dA array, to calculate transition probabilities
  alpha <- array(rep(0, nrow(unique_tend)*nrow(unique_transition)), dim=c(nrow(unique_tend), nrow(unique_transition)))
  I <- array(rep(diag(max(m)), nrow(unique_tend))
             , dim=c( max(m), max(m), nrow(unique_tend)))
  A <- array(0, dim=c(max(m), max(m), nrow(unique_tend)))
  cum_A <- array(0, dim=c(max(m), max(m), nrow(unique_tend)))
  
  # calculate differences in hazards
  alpha <- sapply(1:nrow(unique_transition), function(iter){
    val <- newdata %>% ungroup() %>% filter(transition == unique_transition[iter,1]) %>% arrange(tend)
    val$delta_cumu_hazard
  })
  for (t in 1:nrow(unique_tend)) {
    for (trans in 1:nrow(unique_transition)){
      A[,,t] <- A[,,t] + M[,,trans] * alpha[t, trans]
    }
  }
  
  # #slow code, can be optimized but lack of time for now 
  # for (trans in 1:nrow(unique_transition)){
  #   for (iter in 1:nrow(unique_tend)) {
  #     M[,,trans,iter] <- M[,,trans,iter] * alpha[iter, trans]
  #   }
  # }
  
  # prepare transition probabilities
  A <- I + A
  
  for (iter in 1:nrow(unique_tend)) {
    if (iter == 1) {
      cum_A[,,iter] = A[,,iter]
    } else {
      cum_A[,,iter] = round(cum_A[,,iter-1] %*% A[,,iter],10) #use matrix multiplikation
    }
  }
  
  # transform array so that transition probability can be joined via tend and transition
  tmp <- cbind(unique_tend
               , sapply(1:nrow(unique_transition), function(row) cum_A[unique_transition$from[row] + 1, unique_transition$to[row] + 1, ]))
  colnames(tmp) <- c("tend", as.character(unique_transition$transition))
  trans_prob_df <- tmp %>%
    pivot_longer(cols = c(as.character(unique_transition$transition)), names_to = "transition", values_to = "trans_prob")
  
  # join probabilities and return matrix
  newdata <- newdata %>%
    left_join(trans_prob_df, by=c("tend", "transition")) %>%
    select(-delta_cumu_hazard, -from, -to)
  
  return(newdata)
  
}

# add_trans_prob <- function(
#     newdata
#     #    , object
#     , overwrite       = FALSE 
#     , time_var        = NULL
#     , interval_length = "intlen",
#     ...) {
#   
#   interval_length <- quo_name(enquo(interval_length))
#   
#   if (!overwrite) {
#     if ("trans_prob" %in% names(newdata)) {
#       stop("Data set already contains 'trans_prob' column.
#         Set `overwrite=TRUE` to overwrite")
#     }
#   } else {
#     rm.vars <- intersect(
#       c("trans_prob"
#         # , "surv_lower" # not yet implemented
#         # , "surv_upper" # not yet implemented
#       ),
#       names(newdata))
#     newdata <- newdata %>% select(-one_of(rm.vars))
#   }
#   
#   get_trans_prob(newdata
#                  # , object
#                  , time_var = time_var, interval_length = interval_length, ...)
#   
# }
