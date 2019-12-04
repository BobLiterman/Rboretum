#' Add New Scale
#' Source: eliocamp/new_aes.R
#' @export
bump_aes.Layer <- function(layer, new_aes) {
  original_aes <- new_aes
  
  old_aes <- names(layer$mapping)[remove_new(names(layer$mapping)) %in% new_aes]
  new_aes <- paste0(old_aes, "_new")
  
  old_geom <- layer$geom
  
  old_setup <- old_geom$handle_na
  new_setup <- function(self, data, params) {
    colnames(data)[colnames(data) %in% new_aes] <- original_aes
    old_setup(data, params)
  }
  
  new_geom <- ggplot2::ggproto(paste0("New", class(old_geom)[1]), old_geom,
                               handle_na = new_setup)
  
  new_geom$default_aes <- change_name(new_geom$default_aes, old_aes, new_aes)
  new_geom$non_missing_aes <- change_name(new_geom$non_missing_aes, old_aes, new_aes)
  new_geom$required_aes <- change_name(new_geom$required_aes, old_aes, new_aes)
  new_geom$optional_aes <- change_name(new_geom$optional_aes, old_aes, new_aes)
  
  layer$geom <- new_geom
  
  old_stat <- layer$stat
  
  old_setup2 <- old_stat$handle_na
  new_setup <- function(self, data, params) {
    colnames(data)[colnames(data) %in% new_aes] <- original_aes
    old_setup2(data, params)
  }
  
  new_stat <- ggplot2::ggproto(paste0("New", class(old_stat)[1]), old_stat,
                               handle_na = new_setup)
  
  new_stat$default_aes <- change_name(new_stat$default_aes, old_aes, new_aes)
  new_stat$non_missing_aes <- change_name(new_stat$non_missing_aes, old_aes, new_aes)
  new_stat$required_aes <- change_name(new_stat$required_aes, old_aes, new_aes)
  new_stat$optional_aes <- change_name(new_stat$optional_aes, old_aes, new_aes)
  
  layer$stat <- new_stat
  
  layer$mapping <- change_name(layer$mapping, old_aes, new_aes)
  layer
}