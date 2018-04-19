#' Visualize effect estimates for specific covariate combinations
#'
#' Depending on the plot type and input, creates either a bivariate effect
#' surface or 1-dimensional slices.
#'
#' @import ggplot2
#' @importFrom rlang quos
#' @export
#' @keywords internaldo
gg_partial <- function(data, mod, term, ..., reference = NULL) {

  expressions <- quos(...)
  vars <- names(expressions)
  n_vars <- length(expressions)

  ndf <- make_newdata(data, ...) %>%
    add_term2(mod_dlnm, term, reference=reference)


  gg_base <- ggplot(ndf, aes_string(x = vars[1])) +
    xlab(vars[1])
  if(n_vars == 1) {
    gg_out <- gg_base +
      geom_ribbon(aes_string(ymin = "low", ymax = "high"), alpha = 0.3) +
      geom_line(aes_string(y = "fit"))
  } else {
    if(n_vars == 2) {
      gg_out <- gg_base + aes_string(y = vars[2], z = "fit") +
        geom_tile(aes_string(fill = "fit")) +
        geom_contour(col = "grey30") +
        scale_y_reverse(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0)) +
        scale_fill_gradient2(high = "firebrick2", low = "steelblue")
    }
  }

  gg_out + theme(legend.position = "bottom")
}

#' @rdname gg_partial
#' @importFrom purrr map_int
#' @importFrom rlang quos
#' @export
#' @keywords internal
gg_slice <- function(data, mod, term, ..., reference=NULL) {

  expressions <- quos(...)
  vars        <- names(expressions)
  ndf         <- make_newdata(data, ...) %>%
    add_term2(mod, term, reference = reference)

  n_unique <- map_int(vars, ~length(unique(ndf[[.x]])))
  vars<- vars[rev(order(n_unique))]
  vars <- vars[n_unique[rev(order(n_unique))] > 1]
  ndf <- ndf %>% mutate_at(vars[-1], ~as.factor(.x))
  n_vars <- length(vars)

  gg_out <- ggplot(ndf, aes_string(x = vars[1], y = "fit")) +
    geom_ribbon(aes_string(ymin = "low", ymax = "high"), alpha = 0.3) +
    geom_line()
  if(n_vars > 1) {
    gg_out <- gg_out + aes_string(group=vars[2], fill=vars[2]) +
      geom_line(aes_string(col = vars[2]))
    if(n_vars > 2) {
      form   <- as.formula(paste0("~", vars[-1:-2], collapse="+"))
      gg_out <- gg_out + facet_wrap(form, labeller=label_both)
    }
  }

  gg_out + theme(legend.position = "bottom")

}
