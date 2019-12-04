#' Add New Scale
#' Source: eliocamp/new_aes.R
#' @export
new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}