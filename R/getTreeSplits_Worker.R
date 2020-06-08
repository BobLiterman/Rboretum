#' Rboretum Tree Splitter
#'
#' This function breaks down a rooted phylo or multiPhylo object into its respective set of splits
#' @param tree Rooted phylo object
#' @return Dataframe: Column 1: Semi-colon separated monophyletic clade; Column 2: 'Mirror Clade'; Column 3: Phylo Node ID (NA for root split)
#' @export

getTreeSplits_Worker <- function(tree){
  
  if(!Rboretum::isPhylo(tree,check_rooted = TRUE)){
    stop("'getTreeSplits_Worker' requires a rooted tree.")
  }
  
  # Get species
  tree_species <- sort(tree$tip.label)
  species_count <- length(tree_species)
  
  mono_clades <- c()
  mirror_clades <- c()
  node_list <- c()
  
  tree_subtrees <- ape::subtrees(tree)
  
  tip_list <- list()
  
  # Get all subtrees and create named list for monophyly check
  for(subtree in 1:length(tree_subtrees)){
    temp_subtree <- tree_subtrees[[subtree]]
    temp_tips <- temp_subtree[[2]]
    temp_name <- Rboretum::vectorSemi(sort(temp_tips))
    temp_length <- length(temp_tips)
    tip_list[[temp_name]] <- temp_tips
  }
  
  # Process subtrees
  for(j in tree_subtrees){
    temp_clade <- c((j$tip.label))
    mirror_clade <- tip_list[[Rboretum::vectorSemi(dplyr::setdiff(tree_species, temp_clade))]]
    
    temp_length <- length(temp_clade)
    mirror_length <- length(mirror_clade)
    
    # Remove subtree that is whole tree
    if(temp_length != species_count & mirror_length != species_count){
      
      # Find monophyletic group
      mono_A <- ape::is.monophyletic(tree,temp_clade)
      mono_B <- ape::is.monophyletic(tree,mirror_clade)
      
      # Note actual monophyletic clades and add bootrap values if appropriate
      if(mono_A & !(mono_B)){
        mono_clades <- c(mono_clades,naturalsort(temp_clade) %>% paste(collapse = ";"))
        mirror_clades <- c(mirror_clades,naturalsort(mirror_clade) %>% paste(collapse = ";"))
        node_list <- c(node_list,ape::getMRCA(tree,temp_clade))
      } else if(mono_B & !(mono_A)){
        mono_clades <- c(mono_clades,naturalsort(mirror_clade) %>% paste(collapse = ";"))
        mirror_clades <- c(mirror_clades,naturalsort(temp_clade) %>% paste(collapse = ";"))
        node_list <- c(node_list,ape::getMRCA(tree,mirror_clade))
      } else if(mono_A & mono_B){ # Root clade
        mono_clades <- c(mono_clades,naturalsort(temp_clade) %>% paste(collapse = ";"))
        mirror_clades <- c(mirror_clades,naturalsort(mirror_clade) %>% paste(collapse = ";"))
        node_list <- c(node_list,NA)
      }
    }
  }
  
  split_df <- data.frame("Clade"=as.character(mono_clades),"Mirror_Clade"=as.character(mirror_clades),"Split_Node"=as.integer(node_list)) %>%
    filter((!duplicated(Split_Node))) # Two entries for root reduced to one
  
  return(split_df)
}