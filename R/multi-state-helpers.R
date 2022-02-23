#' @examples
#' t_mat <- matrix(data = NA, nrow = 4, ncol = 4)
#' t_mat[1,2] <- "-1 - sin(x2)"
#' t_mat[1,3] <- "-1.3 + 0.1*x2*x3"
#' t_mat[2,3] <- "-1.5 + 0.8*x1 + 0.3*x2 - 0.1*x3"
#' get_possible_transitions(t_mat)
#' @keywords internal
from_to_pairs <- function(t_mat) {


  n         <- nrow(t_mat)
  from      <- rep(seq_len(n), times = n)
  to        <- rep(seq_len(n), each = n)
  ind_avail <- which(!is.na(t_mat))
  from      <- from[ind_avail]
  to        <- to[ind_avail]

  list(from = from[order(from)], to = to[order(from)])

}

#' @examples
#' ftp <- from_to_pairs2(t_mat)
from_to_pairs2 <- function(t_mat) {

  res <- apply(t_mat, 1, function(x) which(!is.na(x)))
  names(res) <- rownames(t_mat)

  res

}


get_transitions <- function(from_to_list) {

  map_chr(from_to_list, ~paste0(.x, " -> ", .y))

}

#' Add counterfactual observations for possible transitions
#'
#' If data only contains one row per transition that took place, this function
#' adds additional rows for each transition that was possible at that time
#' (for each subject in the data).
#' @param data Data set that only contains rows for transitions that took place.
#' @param from_to_pairs A list with one element for each possible initial state.
#' The values of each list element indicate possible transitions from that state.
#' Will be calculated from the data if unspecified.
#' @param from_col Name of the column that stores initial state.
#' @param to_col Name of the column that stores end state.
#' @param transition_col Name of the column that contains the transition identifier (factor variable).
#' @examples
#'
#' res <- add_counterfactual_transitions(data, list(c(2,3), 3))
#' res <- add_counterfactual_transitions(data)
#' res %>% select(id, from, to, status, tstart, tstop, transition, trans)
#' @export
add_counterfactual_transitions <- function(
  data,
  from_to_pairs  = list(),
  from_col       = "from",
  to_col         = "to",
  transition_col = "transition") {

  if(length(from_to_pairs) == 0L) {
    from_to_pairs <- map(
      .x = unique(data[[from_col]]),
      .f = ~{
        data %>% filter(.data[[from_col]] == .x) %>%
          pull(to_col) %>%
          unique()
      })
  }
  from_to_pairs <- from_to_pairs %>% discard(~length(.x) == 0)

  l_from <- split(data, data[[from_col]])

  orig_status <- data %>%
    select(id, from, to, tstart, status) %>%
    rename(
      orig_to = to,
      orig_status = status)

  res <- map2_dfr(
    from_to_pairs,
    l_from, ~{
      n_to <- length(.x)
      .y %>%
        group_by(id) %>%
        mutate(initial_id = seq_len(n())) %>%
        ungroup() %>%
        slice(rep(row_number(), n_to)) %>%
        arrange(id, from, to) %>%
        group_by(id, initial_id) -> temp
        temp %>%
        mutate(
          to = .x,
          transition = paste0(from, "->", to)
        ) %>%
        ungroup() %>%
        mutate(
          initial_id = NULL,
          transition = as.factor(transition))
    }
  )

  res %>%
    ungroup() %>%
    left_join(orig_status) %>%
    mutate(status = status * (to == orig_to)) %>%
    select(-orig_status, -orig_to) %>%
    arrange(id, tstart, tstop, from, to, status)


}

# ped_msm <- res %>% as_ped(
#   formula = Surv(tstart, tstop, status)~.,
#   transition = "transition",
#   timescale = "calendar")
# data %>%
#   select(id, tstart, tstop, status, transition)
# res %>%
#   select(id, tstart, tstop, status, from, to, transition)
# ped_msm %>%
#   group_by(id, transition) %>%
#   slice(1,n()) %>%
#   select(id, tstart, tend, ped_status, transition) %>%
#   filter(id %in% c(1, 4, 5))

# split_df  %>%
#   mutate(row_number = row_number()) %>%
#   group_by(id, transition) %>%
#   slice(1,n()) %>%
#   select(row_number, id, tstart, tend, ped_status, transition) %>%
#   filter(id %in% c(1, 4, 5))
