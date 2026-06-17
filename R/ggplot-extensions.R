# Stolen from the \code{RmcdrPlugin.KMggplot2} (slightly modified)

#' Step ribbon plots.
#'
#' \code{geom_stepribbon} is an extension of the \code{geom_ribbon}, and
#' is optimized for Kaplan-Meier plots with pointwise confidence intervals
#' or a confidence band. The default \code{direction}-argument \code{"hv"} is 
#' appropriate for right-continuous step functions like the hazard rates etc 
#' returned by \code{pammtools}.
#'
#' @rdname geom_stepribbon
#' @importFrom ggplot2 layer GeomRibbon
#' @inheritParams ggplot2::geom_ribbon
#' @inheritParams ggplot2::geom_step
#' @seealso
#'   \code{\link[ggplot2]{geom_ribbon}} \code{geom_stepribbon}
#' @examples
#' library(ggplot2)
#' huron <- data.frame(year = 1875:1972, level = as.vector(LakeHuron))
#' h <- ggplot(huron, aes(year))
#' h + geom_stepribbon(aes(ymin = level - 1, ymax = level + 1), fill = "grey70") +
#'     geom_step(aes(y = level))
#' h + geom_ribbon(aes(ymin = level - 1, ymax = level + 1), fill = "grey70") +
#'     geom_line(aes(y = level))

#' @export
geom_stepribbon <- function(
    mapping     = NULL,
    data        = NULL,
    stat        = "identity",
    position    = "identity",
    direction   = "hv",
    na.rm       = FALSE,
    show.legend = NA,
    inherit.aes = TRUE, ...) {

  layer(
    data        = data,
    mapping     = mapping,
    stat        = stat,
    geom        = GeomStepribbon,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(na.rm = na.rm, direction = direction, ... )
  )

}

#' @rdname geom_stepribbon
#' @importFrom ggplot2 ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomStepribbon <- ggproto(
  "GeomStepribbon", GeomRibbon,

  extra_params = c("na.rm"),

  draw_group = function(data, panel_scales, coord, na.rm = FALSE, direction = "hv") {

    if (na.rm) data <- data[complete.cases(data[c("x", "ymin", "ymax")]), ]
    data   <- rbind(data, data)
    data   <- data[order(data$x), ]
    data   <- ggplot2_stairstep(data[complete.cases(data["x"]), ], 
                                direction = direction)
    GeomRibbon$draw_group(data, panel_scales, coord, na.rm = na.rm)
  }

)


# code adapted from
# https://github.com/tidyverse/ggplot2/blob/9741da5050f81b7b5c012c56d02f45fc93d68f89/R/geom-path.r#L320
ggplot2_stairstep <- function(data, direction =  c("hv", "vh", "mid")) {
  direction <- match.arg(direction)
  data <- as.data.frame(data)[order(data$x), ]
  n <- nrow(data)
  if (n <= 1) {
    return(data[0, , drop = FALSE])
  }
  if (direction == "vh") {
    xs <- rep(1:n, each = 2)[-2 * n]
    ys <- c(1, rep(2:n, each = 2))
  }
  if (direction == "hv") {
    xs <- c(1, rep(2:n, each = 2))
    ys <- rep(1:n, each = 2)[-2 * n]
  }
  if (direction == "mid") {
    xs <- rep(1:(n - 1), each = 2)
    ys <- rep(1:n, each = 2)
  }

  ymin <- c(data$ymin[ys])
  ymax <- c(data$ymax[ys])
  if (direction == "mid") {
    gaps <- data$x[-1] - data$x[-n]
    mid_x <- data$x[-n] + gaps/2
    x <- c(data$x[1], mid_x[xs], data$x[n])
    data_attr <- data[c(1, xs, n),
                      setdiff(names(data), c("x", "ymin", "ymax"))]
  } else {
    x <- data$x[xs]
    ymin <- data$ymin[ys]
    ymax <- data$ymax[ys]
    data_attr <- data[xs, setdiff(names(data), c("x", "ymin", "ymax"))]
  }
  cbind(data.frame(x = x, ymin = ymin, ymax = ymax), data_attr)
}

#' Plot State Occupation Probabilities
#'
#' Creates a stacked area plot of state occupation probabilities over time,
#' computed from transition probability matrices stored as an attribute of
#' the input data. Optionally facets by a grouping variable.
#'
#' @param newdata A data frame with an attribute \code{matrix} containing
#'   a data frame with a column \code{trans_prob_matrix}. Each element of
#'   \code{trans_prob_matrix} should be a 3-dimensional array of dimensions
#'   \code{n_states x n_states x n_timepoints}.
#' @param init_state A numeric vector specifying the initial state distribution.
#'   Should sum to 1 and have length equal to the number of states. For example,
#'   \code{c(0, 1, 0, 0)} places all subjects in state 2 at baseline.
#' @param group_var A character string giving the name of the column in
#'   \code{newdata} to facet by (e.g., \code{"treat"}). If \code{NULL}
#'   (default), no faceting is applied.
#' @param time_var A character string giving the name of the time variable in
#'   \code{newdata}. Defaults to \code{"tend"}.
#' @param ncol An integer specifying the number of columns in the facet wrap.
#'   If \code{NULL} (default), defaults to the number of unique groups.
#'
#' @return A \code{ggplot} object showing stacked-area state occupation
#'   probabilities over time, optionally faceted by \code{group_var}.
#' @export
gg_state_occupation <- function(
    newdata,
    init_state,
    group_var = NULL,
    time_var  = "tend",
    ncol      = NULL
) {
  # Extract attribute matrix
  mat_df <- attributes(newdata)$matrix
  time_var <- sym(time_var)

  # error handler: init state number must align with states in trans prob
  n_states <- dim(mat_df$trans_prob_matrix[[1]])[1]
  if (length(init_state) != n_states) {
    stop("`init_state` must have length equal to the number of states (", n_states, ").")
  }

  if (!is.null(group_var)) {

    labels <- attr(newdata[[group_var]], "labels")
    if (is.null(ncol)) ncol <- length(unique(newdata[[group_var]]))

    if (!is.null(labels)) {
      newdata[[group_var]] <- factor(
        newdata[[group_var]],
        levels = names(labels),
        labels = labels
      )
    }
    
    # Iterate over groups
    df_all <- mat_df %>%
      mutate(
        df_long = map(.data[["trans_prob_matrix"]], function(x) {
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
      select(-dplyr::all_of("trans_prob_matrix")) %>%
      unnest(cols = "df_long")
    
    # plot
    p <- ggplot(df_all, aes(x = .data[["time"]], y = .data[["prob"]], fill = .data[["state"]])) +
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
    
    p <- ggplot(df_all, aes(x = .data[["time"]], y = .data[["prob"]], fill = .data[["state"]])) +
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

