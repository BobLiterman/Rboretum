#' Get All Clades
#'
#' This function takes a rooted tree and returns a vector containing each monophyletic group, including both root groups
#' @param rooted_tree Rooted phylo object
#' @param exclude_root OPTIONAL: If TRUE, exclude root clades from returned list. Otherwise, return all groups (Default: FALSE)
#' @return Character vector of monophyletic clades
#' @export
#' @examples
#' getAllClades(rooted_tree)
#'

getAllClades <- function(rooted_tree,exclude_root){

  if(missing(exclude_root)){
    exclude_root <- FALSE
  }
  if(exclude_root != TRUE){
    exclude_root <- FALSE
  }

  if(!(ape::is.rooted.phylo(rooted_tree))){
    stop("Function getAllClades() requires rooted tree.")
  }

  tree_spp <- sort(rooted_tree$tip.label)
  if(length(tree_spp) <= 3){
    stop("Tree must contain at least 3 species")
  }

  splits <- Rboretum::getAllSplits(rooted_tree)
  temp_nonroot <- splits %>% filter(!is.na(Split_Node)) %>% pull(Clade) %>% unlist() %>% as.character()

  if(!exclude_root){
    temp_root_1 <- splits %>% filter(is.na(Split_Node)) %>% pull(Clade) %>% unlist() %>% as.character()
    temp_root_2 <- splits %>% filter(is.na(Split_Node)) %>% pull(Mirror_Clade) %>% unlist() %>% as.character()
    all_clades <- c(temp_nonroot,temp_root_1,temp_root_2)
    return(all_clades)
  }

  else{
    return(temp_nonroot)
  }
}
