#'  Transform PED data back to classic format
#'
#' @param ped A data set in PED format
#' @param id_var Name of the id column
#' @export
#' @keywords internal
unsplit <- function(ped, id_var = "id") {
  last_ind <- cumsum(rle(ped[[id_var]])$lengths)

  uns_df <- ped[last_ind, ]
  uns_df$time <- uns_df$tstart + exp(uns_df$offset)
  uns_df$status <- uns_df$ped_status
  uns_df$interval <- uns_df$start <- uns_df$offset <- uns_df$tend <-
    uns_df$ped_status <- uns_df$intmid <- NULL

  return(uns_df)

}
