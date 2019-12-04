#' Add New Scale
#' Source: eliocamp/new_aes.R
#' @export
remove_new <- function(aes) {
  stringi::stri_replace_all(aes, "", regex = "(_new)*")
}