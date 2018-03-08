#' Create nested data frame from event and time-dependent data
#'
#' @inheritParams tidyr::nest
#' @inherit split_tdc
#' @import checkmate dplyr
#' @importFrom tidyr nest
#' @importFrom tibble as_tibble
#' @importFrom purrr invoke_map map map_int reduce
#' @export
nest_tdc <- function(
  event_df,
  tdc_df,
  te_var,
  id_var     = "id",
  time_var   = "time",
  status_var = "status",
  entry_time = 0,
  cens_value = 0,
  .key       = "tdc") {


  assert_data_frame(event_df)
  assert_data_frame(tdc_df)
  assert_subset(c(id_var, time_var, status_var), colnames(event_df))
  assert_subset(c(id_var, te_var), colnames(tdc_df))


  common_id <- warn_partial_overlap(event_df[[id_var]], tdc_df[[id_var]])
  event_df  <- filter(event_df, !!sym(id_var) %in% common_id)
  tdc_df    <- filter(tdc_df, !!sym(id_var) %in% common_id)

  tdc_vars <- setdiff(
    get_tdc(tdc_df, id_var),
    c(id_var, status_var, time_var, te_var))

  nested_tdc <- tdc_df %>%
    as_tibble() %>%
    select(!!id_var, !!te_var, one_of(tdc_vars)) %>%
    nest(one_of(te_var, tdc_vars), .key = !!.key)

  tdc_in_event_df <- intersect(tdc_vars, names(event_df))
  if (!(length(tdc_in_event_df)==0)) {
    event_df <- event_df %>%  select(-one_of(tdc_vars))
  }

  # nested_tdc <- left_join(event_df, nested_tdc)

  nested_tdc <- event_df %>% left_join(nested_tdc, by = id_var)

  # important to calculate unique event times after subsetting
  utime <- unique(nested_tdc[[time_var]]*nested_tdc[[status_var]]) %>% sort()
  ute_time <- tdc_df %>% pull(te_var) %>% unique() %>% sort()

  # utime <- union(utime, ute_time) %>% unique() %>% sort()

  attr(nested_tdc, "id_var")     <- id_var
  attr(nested_tdc, "time_var")   <- time_var
  attr(nested_tdc, "status_var") <- status_var
  attr(nested_tdc, "te_var")     <- te_var
  attr(nested_tdc, "tdc_col")    <- .key
  attr(nested_tdc, "cens_value") <- cens_value
  attr(nested_tdc, "breaks")     <- utime
  attr(nested_tdc, "te")         <- ute_time
  attr(nested_tdc, "id_n") <- nested_tdc %>% pull(time_var) %>%
    pmin(max(utime)) %>%
    map_int(findInterval, vec=utime, left.open=TRUE, rightmost.closed=TRUE)
  # attr(nested_tdc, "id_te_n") <- map_int(nested_tdc[["tdc"]], ~nrow(.))
  attr(nested_tdc, "id_tseq") <- attr(nested_tdc, "id_n") %>%
    map(seq_len) %>% unlist()
  attr(nested_tdc, "id_teseq") <- rep(seq_along(common_id),
    times=attr(nested_tdc, "id_n"))

  class(nested_tdc) <- c("ntdc_df", class(nested_tdc))

  nested_tdc

}


# #' Transform a nested data frame to data frame with matrix columns
# #'
# #' @param nested_df A nested data frame with list columns that contain
# #' a matrix per row.
# #' @importFrom dplyr select_if slice
# #' @importFrom purrr compose
# #' @keywords internal
# as_matdf <- function(nested_df) {

# 	ndf  <- select_if(nested_df, compose("!", is.list))
# 	ldf  <- select_if(nested_df, is.list)
# 	nrep <- sapply(ldf[[1]], function(z) nrow(z))
# 	rind <- lapply(seq_len(nrow(ndf)), function(z) rep(z, nrep[z])) %>%
# 		unlist()

# 	ldf <- lapply(ldf, function(z) do.call(rbind, z))
# 	ndf <- ndf %>% slice(rind) %>% as.data.frame()
# 	for(i in seq_along(ldf)) {
# 		ndf[[names(ldf)[i]]] <- ldf[[i]]
# 	}

# 	return(ndf)

# }
