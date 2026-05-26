gg_state_occupation <- function(
    newdata, 
    init_state, 
    group_var = NULL, 
    group_labels = NULL,
    time_var = "tend",
    ncol = NULL,
    nrow = NULL
    ) {
  # newdata: your data object with attributes(newdata)$matrix
  # group_var: column name to facet by, e.g., "treat"
  # init_state: numeric vector for initial state, e.g., c(0,1,0,0)

  # Extract attribute matrix
  mat_df <- attributes(newdata)$matrix
  time_var <- sym(time_var)
  
  labels <- attr(newdata[[group_var]], "labels")
  ncol <- length(unique(newdata[[group_var]]))
  if (!is.null(labels)) {
    newdata[[group_var]] <- factor(
      newdata[[group_var]],
      levels = names(labels),
      labels = labels
    )
  }
  
  if (!is.null(group_var)) {
  # Iterate over groups
  df_all <- mat_df %>%
    mutate(
      df_long = map(trans_prob_matrix, function(x) {
        # x is a 4 × 4 × T array for ONE group
        
        res <- apply(x, 3, function(mat) init_state %*% mat)
        
        df <- as.data.frame(t(res))
        colnames(df) <- paste0("state_", seq_len(ncol(df)))
        df$time <- unique(sort(newdata |> dplyr::pull(time_var)))
        
        df |>
          pivot_longer(
            cols = starts_with("state_"),
            names_to = "state",
            values_to = "prob"
          )
      })
    ) %>%
    select(-trans_prob_matrix) %>%
    unnest(df_long)
  
    # plot
    p <- ggplot(df_all, aes(x = time, y = prob, fill = state)) +
      geom_area(color = "black", alpha = 0.8) +
      facet_wrap(vars(.data[[group_var]]), ncol = ncol) +
      labs(
        x = "Time",
        y = "State occupation probability",
        fill = "State"
      ) +
      theme_minimal()
  } else {
    # process without grouping
    x <- mat_df$trans_prob_matrix[[1]]  # assume single matrix
    res <- apply(x, 3, function(mat) init_state %*% mat)
    df_all <- as.data.frame(t(res))
    colnames(df_all) <- paste0("state_", seq_len(ncol(df_all)))
    df_all$time <- unique(sort(newdata |> dplyr::pull(time_var)))
    
    df_all <- df_all |> pivot_longer(
      cols = starts_with("state_"),
      names_to = "state",
      values_to = "prob"
    )
    
    p <- ggplot(df_all, aes(x = time, y = prob, fill = state)) +
      geom_area(color = "black", alpha = 0.8) +
      labs(
        x = "Time",
        y = "State occupation probability",
        fill = "State"
      ) +
      theme_minimal()
  }
  
  return(p)
  
}

# Example usage:
# gg_state_occupation(newdata, group_var = "treat", init_state = c(0,1,0,0))