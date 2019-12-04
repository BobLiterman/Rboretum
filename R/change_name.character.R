#' Add New Scale
#' Source: eliocamp/new_aes.R
#' @export
change_name.character <- function(list, old, new) {
  list[list %in% old] <- new
  list
}