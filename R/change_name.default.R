#' Add New Scale
#' Source: eliocamp/new_aes.R
#' @export
change_name.default <- function(list, old, new) {
  nam <- names(list)
  nam[nam %in% old] <- new
  names(list) <- nam
  list
}