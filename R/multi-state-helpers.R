#' Extract transition information from different objects
#'
#' @rdname from_to_pairs
#' @param t_mat an object that contains information about possible transitions.
#' @keywords internal
from_to_pairs <- function(t_mat, ...) {

  UseMethod("from_to_pairs", t_mat)

}

#' @rdname from_to_pairs
#' @keywords internal
from_to_pairs2 <- function(t_mat, ...) {

  res <- apply(t_mat, 1, function(x) which(x) - 1)
  names(res) <- seq_len(nrow(t_mat)) - 1
  res <- res[vapply(res, length, 0) != 0]

  res

}

#' @rdname from_to_pairs
#' @param from_col The name of the column in the data frame that contains "from" states.
#' @param to_col The name of the column in the data frame that contains "to" states.
#' @keywords internal
#' @examples
#' \dontrun{
#' df = data.frame(id = c(1,1, 2,2), from = c(1, 1, 2, 2), to = c(2, 3, 2, 2))
#' from_to_pairs(df)
#' }
from_to_pairs.data.frame <- function(t_mat, from_col = "from", to_col = "to", ...) {

  map(
    .x = sort(unique(t_mat[[from_col]])),
    .f = ~{
      t_mat %>%
        filter(.data[[from_col]] == .x) %>%
        pull(to_col) %>%
        unique()
    })

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
#' @export
add_counterfactual_transitions <- function(
  data,
  from_to_pairs  = list(),
  from_col       = "from",
  to_col         = "to",
  transition_col = "transition") {


  if(length(from_to_pairs) == 0) {
    from_to_pairs <- from_to_pairs.data.frame(data, from_col, to_col) %>%
      discard(~length(.x) == 0)
  }

  l_from <- split(data, data[[from_col]])

  orig_status <- data %>%
    select(one_of(c("id", "from", "to", "tstart", "status"))) %>%
    rename(
      "orig_to" = "to",
      "orig_status" = "status")

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
    mutate(status = .data$status * (.data$to == .data$orig_to)) %>%
    select(-one_of("orig_status", "orig_to")) %>%
    arrange(.data$id, .data$tstart, .data$tstop, .data$from, .data$to, .data$status)


}
