#' 
#' #' Simulate data for multi-state models
#' #'
#' #' @rdname sim_pexp
#' #' @inherit sim_pexp
#' #' @param tmat Transition matrix containing the log-hazard definition for each transition.
#' #' Transitions that can not occurr should contain a \code{NA}.
#' #' @param proportion_censoring The proportion of censored observations
#' #' @param keep_transitions_at_risk Logical. If \code{TRUE} (default), all transitions
#' #' for which a subject was at risk will be included in the data set.
#' #' @export
#' #' @keywords internal
#' #' @examples
#' #' t_mat <- matrix(data = NA, nrow = 3, ncol = 3)
#' #' t_mat[1,2] <- "log(0.7) + x1"
#' #' t_mat[1,3] <- "log(0.5) - 0.25 * x2"
#' #' t_mat[2,3] <- "log(0.9)"
#' #'
#' #' n = 100
#' #' data <- cbind.data.frame(
#' #'  id = 1:n,
#' #'  x1 = runif(n, -3, 3),
#' #'  x2 = runif(n, 0, 6),
#' #'  from = 1,
#' #'  t = 0)
#' #'
#' #' cut =  seq(0, 3, by = 0.01)
#' #'
#' #' msm_df <- sim_pexp_msm(
#' #'  t_mat = t_mat,
#' #'  data = data,
#' #'  cut = cut,
#' #'  proportoin_censoring = 0.3,
#' #'  keep_transitions_at_risk = TRUE)
#' #' head(msm_df)
#' #'
#' sim_pexp_msm <- function(
#'     t_mat,
#'     data,
#'     cut,
#'     proportion_censoring = 0.3,
#'     keep_transitions_at_risk = FALSE ) {
#'   
#'   # Aus der Transition Matrix können die Transition-States abgelesen werden
#'   trans_states <- which(apply(t_mat,1, function(x) !all(is.na(x))))
#'   
#'   # cens_time <- c(
#'   #   rep(max(cut), times = floor( (1 - proportoin_censoring)*nrow(data))),
#'   #   runif((nrow(data) - floor((1 - proportoin_censoring)*nrow(data))), 0, max(cut))
#'   # )
#'   # data$cens_time <- sample(cens_time)
#'   
#'   # Initialize object where results are stored
#'   results <- NULL
#'   sim_df <- NULL
#'   trans_not_made <- NULL
#'   
#'   
#'   while (nrow(data) > 0) {
#'     
#'     # Für jeden Transition state
#'     for (i in seq_along(trans_states)) {
#'       # Welche Beobachtungen befinden sich in diesem Status
#'       if (nrow(data[data$from == i,]) == 0) next
#'       # Welche Hazards führen von diesem Status weg
#'       form <- as.Formula(paste("~ ", paste(t_mat[i,which(!is.na(t_mat[i,]))], collapse = "|")))
#'       # Führe für jede der Beobachtungen in dem state ein CR Experiment durch
#'       # Speichere die Ergebnisse in einem Dataframe
#'       results <- rbind(results,
#'                        sim_pexp_cr(form, data[data$from == i,], cut = cut) %>%
#'                          mutate(
#'                            time = time + t,
#'                            to =  which(!is.na(t_mat[i,]))[type]
#'                          ) %>%
#'                          select(-one_of(paste0("hazard",seq_len(length(attr(form,"rhs"))))))
#'       )
#'       
#'       
#'       if (keep_transitions_at_risk) {
#'         # Alle Transitions, die nicht gemacht wurden, aber theoretisch möglich gewesen wären, sollen als zensiert im Dataframe auftauchen (Vgl. mstate)
#'         trans_not_made <- rbind(trans_not_made,
#'                                 results %>%
#'                                   slice(rep(row_number(),length(attr(form, "rhs")) )) %>%
#'                                   mutate(
#'                                     pos_types = seq_len(length(attr(form, "rhs"))),
#'                                     status = if_else(type != pos_types, 0, 1),
#'                                     to =  which(!is.na(t_mat[i,]))[pos_types]
#'                                   ) %>%
#'                                   filter(
#'                                     type != pos_types
#'                                   ) %>%
#'                                   mutate(
#'                                     pos_types = NULL
#'                                   )
#'                                 
#'         )
#'         
#'       }
#'       
#'     }
#'     
#'     sim_df <- rbind(sim_df, results, trans_not_made)
#'     
#'     data <- data %>%
#'       left_join(
#'         results[,c("id", "to", "time","status")], by = "id") %>%
#'       filter(
#'         status != 0 & to %in% trans_states) %>%
#'       mutate(
#'         from   = to,
#'         to     = NULL,
#'         t      = t + time,
#'         time   = NULL,
#'         status = NULL
#'       )
#'     # results und trans_not_made wieder leeren
#'     results <- NULL
#'     trans_not_made <- NULL
#'     
#'     # Ergebnisframe sortieren, ids untereinander und nach Transitionzeiten sortieren
#'     sim_df <- sim_df %>%
#'       group_by(id) %>%
#'       arrange(time, .by_group = TRUE)
#'     # Solange, wie sich Individuen in data befinden gibt es id's in transient States, und das Vorgehen wird wiederholt
#'   }
#'   
#'   # sim_df aufräumen
#'   sim_df <- sim_df %>%
#'     mutate(
#'       status = status * (time <= max(cut)),
#'       time   = pmin(time, max(cut)),
#'       # to     = ifelse(status, to, from),
#'       tstart = t,
#'       tstop  = time,
#'       gap    = tstop - tstart,
#'       time   = NULL,
#'       t      = NULL,
#'       type   = NULL,
#'       transition = as.factor(paste0(from, "->", to))
#'     ) %>%
#'     relocate(to, .after = from)
#'   
#'   return(sim_df)
#'   
#' }
#' 
#' 
#' 
#' #' A formula special used to handle cumulative effect specifications
#' #'
#' #' Can be used in the second part of the formula specification provided
#' #' to \code{\link[pammtools]{sim_pexp}} and should only be used in this
#' #' context.
#' #'
#' #' @importFrom purrr map
#' #' @export
#' #' @keywords internal
#' fcumu <- function(..., by = NULL, f_xyz, ll_fun) {
#'   
#'   vars   <- as.list(substitute(list(...)))[-1] %>%
#'     map(~as.character(.x)) %>%
#'     unlist()
#'   vars <- vars[vars != "t"]
#'   
#'   list(
#'     vars   = vars,
#'     f_xyz  = f_xyz,
#'     ll_fun = ll_fun)
#'   
#' }
#' 
#' #' @import dplyr
#' #' @importFrom tidyr unnest
#' #' @importFrom rlang sym :=
#' #' @keywords internal
#' eta_cumu <- function(data, fcumu, cut, ...) {
#'   
#'   vars   <- fcumu$vars
#'   f_xyz  <- fcumu$f_xyz
#'   ll_fun <- fcumu$ll_fun
#'   eta_name <- paste0("eta_", vars[2])
#'   comb_df <- combine_df(
#'     data.frame(t = cut),
#'     select(data, one_of("id", vars)))
#'   comb_df <- comb_df %>% unnest(cols = -one_of("id"))
#'   comb_df %>%
#'     group_by(.data$id, .data$t) %>%
#'     mutate(
#'       LL = ll_fun(t, !!sym(vars[1])) * 1,
#'       delta = c(mean(abs(diff(!!sym(vars[1])))), abs(diff(!!sym(vars[1]))))) %>%
#'     ungroup() %>%
#'     filter(.data$LL != 0) %>%
#'     group_by(.data$id, .data$t) %>%
#'     summarize(!!eta_name :=
#'                 sum(.data$delta * f_xyz(.data$t, .data[[vars[1]]], .data[[vars[2]]])))
#'   
#' }
#' 
#' #' Simulate data for competing risks scenario
#' #'
#' #'
#' #' @keywords internal
#' sim_pexp_cr <- function(formula, data, cut) {
#'   
#'   # Formula extends the base class formula by allowing for multiple responses and multiple parts of regressors
#'   Form    <- Formula(formula)
#'   # Extract the right handside of the Formula
#'   F_rhs   <- attr(Form, "rhs")
#'   l_rhs   <- length(F_rhs)
#'   seq_rhs <- seq_len(l_rhs)
#'   
#'   if (!("id" %in% names(data))) {
#'     data$id <- 1:(nrow(data))
#'   }
#'   
#'   if (!("t" %in% names(data))) {
#'     data$t <- 0
#'   }
#'   
#'   data <- data %>%
#'     mutate(
#'       time   = max(cut),
#'       status = 1
#'     )
#'   
#'   # construct eta for time-constant part
#'   # offset (the log of the duration during which the subject was under risk in that interval)
#'   
#'   ped  <- split_data(
#'     formula = Surv(time, status)~.,
#'     data    = select_if(data, is_atomic),
#'     cut     = cut,
#'     id      = "id") %>%
#'     mutate(
#'       t = t + tstart
#'     )
#'   
#'   # calculate cause specific hazards
#'   
#'   for (i in seq_rhs) {
#'     ped[[paste0("hazard", i)]] <-  exp(eval(F_rhs[[i]], ped))
#'   }
#'   ped[["rate"]] <- reduce(ped[paste0("hazard", seq_rhs)], `+`)
#'   
#'   # simulate survival times
#'   
#'   sim_df <- ped %>%
#'     group_by(id) %>%
#'     mutate(
#'       time   = rpexp(rate = .data$rate, t = .data$tstart),
#'       status = 1L * (.data$time <= max(cut)),
#'       time   = pmin(.data$time, max(cut)),
#'       # t wieder ins "Original" zurückrechnen, muss später auf die Waitingtime drauf gerechnet werden
#'       t = .data$t - .data$tstart
#'     ) %>%
#'     filter(.data$tstart < .data$time & .data$time <= .data$tend)
#'   
#'   
#'   
#'   # Ziehe aus den möglichen hazards eins mit den entsprechenden Wahrscheinlichkeiten
#'   sim_df$type <- apply(sim_df[paste0("hazard", seq_rhs)], 1,
#'                        function(probs)
#'                          sample(seq_rhs, 1, prob = probs))
#'   
#'   sim_df %>%
#'     select(-one_of(c("tstart", "tend", "interval", "offset", "ped_status", "rate")))
#'   
#' }
#' 
#' 
#' # ------------------------------------------------------------------------------
#' #
#' # EXAMPLE FROM VIGNETTE / R HELPER
#' #
#' # ------------------------------------------------------------------------------
#' library(survival)
#' library(dplyr)
#' # library(pammtools)
#' 
#' # set number of observations/subjects
#' n <- 250
#' # create data set with variables which will affect the hazard rate.
#' df <- cbind.data.frame(x1 = runif (n, -3, 3), x2 = runif (n, 0, 6)) %>%
#'   as_tibble()
#' # the formula which specifies how covariates affect the hazard rate
#' f0 <- function(t) {
#'   dgamma(t, 8, 2) *6
#' }
#' form <- ~ -3.5 + f0(t) -0.5*x1 + sqrt(x2)
#' set.seed(24032018)
#' sim_df <- sim_pexp(form, df, 1:10)
#' head(sim_df)
#' plot(survfit(Surv(time, status)~1, data = sim_df ))
#' 
#' # for control, estimate with Cox PH
#' mod <- coxph(Surv(time, status) ~ x1 + pspline(x2), data=sim_df)
#' coef(mod)[1]
#' layout(matrix(1:2, nrow=1))
#' termplot(mod, se = TRUE)
#' 
#' # and using PAMs
#' layout(1)
#' ped <- sim_df %>% as_ped(Surv(time, status)~., max_time=10)
#' library(mgcv)
#' pam <- gam(ped_status ~ s(tend) + x1 + s(x2), data=ped, family=poisson, offset=offset)
#' coef(pam)[2]
#' plot(pam, page=1)


# ------------------------------------------------------------------------------
#
# DEBUG VERSION OF SIM_PEXP_MSM FOR TESTING
#
# ------------------------------------------------------------------------------

sim_pexp_msm_debug <- function(
    t_mat,
    data,
    cut,
    proportion_censoring = 0.3,
    keep_transitions_at_risk = FALSE ) {
  
  # Aus der Transition Matrix können die Transition-States abgelesen werden
  trans_states <- which(apply(t_mat,1, function(x) !all(is.na(x))))
  
  # cens_time <- c(
  #   rep(max(cut), times = floor( (1 - proportoin_censoring)*nrow(data))),
  #   runif((nrow(data) - floor((1 - proportoin_censoring)*nrow(data))), 0, max(cut))
  # )
  # data$cens_time <- sample(cens_time)
  
  # Initialize object where results are stored
  results <- NULL
  sim_df <- NULL
  trans_not_made <- NULL
  
  
  while (nrow(data) > 0) {
    
    # Für jeden Transition state
    for (i in seq_along(trans_states)) {
      # Welche Beobachtungen befinden sich in diesem Status
      if (nrow(data[data$from == i,]) == 0) next
      # Welche Hazards führen von diesem Status weg
      form <- as.Formula(paste("~ ", paste(t_mat[i,which(!is.na(t_mat[i,]))], collapse = "|")))
      # Führe für jede der Beobachtungen in dem state ein CR Experiment durch
      # Speichere die Ergebnisse in einem Dataframe
      results <- rbind(results,
                       sim_pexp_cr(form, data[data$from == i,], cut = cut) %>%
                         mutate(
                           time = time + t,
                           to =  which(!is.na(t_mat[i,]))[type]
                         ) %>%
                         select(-one_of(paste0("hazard",seq_len(length(attr(form,"rhs"))))))
      )
      
      
      if (keep_transitions_at_risk) {
        # Alle Transitions, die nicht gemacht wurden, aber theoretisch möglich gewesen wären, sollen als zensiert im Dataframe auftauchen (Vgl. mstate)
        trans_not_made <- rbind(trans_not_made,
                                results %>%
                                  slice(rep(row_number(),length(attr(form, "rhs")) )) %>%
                                  mutate(
                                    pos_types = seq_len(length(attr(form, "rhs"))),
                                    status = if_else(type != pos_types, 0, 1),
                                    to =  which(!is.na(t_mat[i,]))[pos_types]
                                  ) %>%
                                  filter(
                                    type != pos_types
                                  ) %>%
                                  mutate(
                                    pos_types = NULL
                                  )
                                
        )
        
      }
      
    }
    
    sim_df <- rbind(sim_df, results, trans_not_made)

    data <- data %>%
      left_join(
        results[,c("id", "to", "time","status")], by = "id") %>%
      filter(
        status != 0 & to %in% trans_states) %>%
      mutate(
        from   = to,
        to     = NULL,
        t      = t + time,
        time   = NULL,
        status = NULL
      )
    # results und trans_not_made wieder leeren
    results <- NULL
    trans_not_made <- NULL
    
    # Ergebnisframe sortieren, ids untereinander und nach Transitionzeiten sortieren
    sim_df <- sim_df %>%
      group_by(id) %>%
      arrange(time, .by_group = TRUE)
    # Solange, wie sich Individuen in data befinden gibt es id's in transient States, und das Vorgehen wird wiederholt
  }
  
  # sim_df aufräumen
  sim_df <- sim_df %>%
    mutate(
      status = status * (time <= max(cut)),
      time   = pmin(time, max(cut)),
      # to     = ifelse(status, to, from),
      tstart = t,
      tstop  = time,
      gap    = tstop - tstart,
      time   = NULL,
      t      = NULL,
      type   = NULL,
      transition = as.factor(paste0(from, "->", to))
    ) %>%
    relocate(to, .after = from)
  
  return(sim_df)
  
}


