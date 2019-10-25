#' Get All Splits
#'
#' This function takes a rooted tree and returns information about each monophyletic group/split
#' @param rooted_tree Rooted phylo object
#' @param exclude_root OPTIONAL: If TRUE, exclude root clades. Otherwise, return all groups (Default: FALSE)
#' @return Three-column dataframe; 1: Monophyletic clade; 2: 'Mirror Clade'; 3: Node ID
#' @export
#' @examples
#' getAllSplits(rooted_tree)
#'

getAllSplits <- function(rooted_tree,exclude_root){

  if(missing(exclude_root)){
    exclude_root <- FALSE
  }
  if(exclude_root != TRUE){
    exclude_root <- FALSE
  }

  if(!(ape::is.rooted.phylo(rooted_tree))){
    stop("Function getAllSplits() requires rooted tree.")
  }

  # Get species
  tree_species <- rooted_tree$tip.label
  species_count <- length(tree_species)

  non_root_clades <- c()
  non_root_mirror <- c()
  node_list <- c()

  # Get all subtrees
  for(j in ape::subtrees(rooted_tree)){
    temp_clade <- c((j$tip.label))
    mirror_clade <- dplyr::setdiff(tree_species, temp_clade)

    temp_length <- length(temp_clade)
    mirror_length <- length(mirror_clade)

    # Remove subtree that is whole tree
    if(temp_length != species_count & mirror_length != species_count){

      # Find monophyletic group
      mono_A <- ape::is.monophyletic(rooted_tree,temp_clade)
      mono_B <- ape::is.monophyletic(rooted_tree,mirror_clade)

      # Note actual monophyletic clades and add bootrap values if appropriate
      if(mono_A & !(mono_B)){
        non_root_clades <- c(non_root_clades,sort(temp_clade) %>% paste(collapse = ";"))
        non_root_mirror <- c(non_root_mirror,sort(mirror_clade) %>% paste(collapse = ";"))
        node_list <- c(node_list,ape::getMRCA(rooted_tree,temp_clade))
      }
      if(mono_B & !(mono_A)){
        non_root_clades <- c(non_root_clades,sort(mirror_clade) %>% paste(collapse = ";"))
        non_root_mirror <- c(non_root_mirror,sort(temp_clade) %>% paste(collapse = ";"))
        node_list <- c(node_list,ape::getMRCA(rooted_tree,mirror_clade))
      }
      # If root...
      if(mono_A & mono_B){
        non_root_clades <- c(non_root_clades,sort(temp_clade) %>% paste(collapse = ";"))
        non_root_mirror <- c(non_root_mirror,sort(mirror_clade) %>% paste(collapse = ";"))
        node_list <- c(node_list,NA)
      }
    }
  }

  clade_node_df <- data.frame("Clade"=as.character(non_root_clades),"Mirror_Clade"=as.character(non_root_mirror),"Split_Node"=as.integer(node_list)) %>%
    filter((!duplicated(Split_Node))) # Two entries for root reduced to one

  if(exclude_root){
    clade_node_df <- clade_node_df %>% filter(!is.na(Split_Node))
  }

  return(clade_node_df)
}
