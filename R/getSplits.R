#' Rboretum Tree Splitter
#'
#' This function takes a phylo object and returns information about each monophyletic group/split
#' @param tree Phylo object
#' @return Four-column dataframe; 1: Semi-colon separated monophyletic clade; 2: 'Mirror Clade'; 3: Phylo Node ID; 4: Node Boostrap (if present)
#' @export
#' @examples
#' getSplits(tree)
#'

getSplits <- function(tree){

  # Get species
  tree_species <- sort(tree$tip.label)
  species_count <- length(tree_species)

  non_root_clades <- c()
  non_root_mirror <- c()
  node_list <- c()
  bootstrap_list <- c()

  # Get all subtrees
  for(j in ape::subtrees(tree)){
    temp_clade <- c((j$tip.label))
    mirror_clade <- dplyr::setdiff(tree_species, temp_clade)

    temp_length <- length(temp_clade)
    mirror_length <- length(mirror_clade)

    # Remove subtree that is whole tree
    if(temp_length != species_count & mirror_length != species_count){

      # Find monophyletic group
      mono_A <- ape::is.monophyletic(tree,temp_clade)
      mono_B <- ape::is.monophyletic(tree,mirror_clade)
      
      # Note actual monophyletic clades and add bootrap values if appropriate
      if(mono_A & !(mono_B)){
        non_root_clades <- c(non_root_clades,sort(temp_clade) %>% paste(collapse = ";"))
        non_root_mirror <- c(non_root_mirror,sort(mirror_clade) %>% paste(collapse = ";"))
        node_list <- c(node_list,ape::getMRCA(tree,temp_clade))
        
        bs_tree <- Rboretum::getTrimmedTree(tree,temp_clade)
        node_bs <- bs_tree$node.label[1]
        if(!is.na(as.numeric(node_bs))){
          bootstrap_list <- c(bootstrap_list,as.numeric(node_bs))
        } else{ bootstrap_list <- c(bootstrap_list,NA) }
      }
      if(mono_B & !(mono_A)){
        non_root_clades <- c(non_root_clades,sort(mirror_clade) %>% paste(collapse = ";"))
        non_root_mirror <- c(non_root_mirror,sort(temp_clade) %>% paste(collapse = ";"))
        node_list <- c(node_list,ape::getMRCA(tree,mirror_clade))
        
        bs_tree <- Rboretum::getTrimmedTree(tree,temp_clade)
        node_bs <- bs_tree$node.label[1]
        if(!is.na(as.numeric(node_bs))){
          bootstrap_list <- c(bootstrap_list,as.numeric(node_bs))
        } else{ bootstrap_list <- c(bootstrap_list,NA) }
      }
      # If root...
      if(mono_A & mono_B){
        non_root_clades <- c(non_root_clades,sort(temp_clade) %>% paste(collapse = ";"))
        non_root_mirror <- c(non_root_mirror,sort(mirror_clade) %>% paste(collapse = ";"))
        node_list <- c(node_list,NA)
        
        root_1 <- sort(temp_clade)
        root_tree_1 <- Rboretum::getTrimmedTree(tree,root_1)
        root_1_BS <- root_tree_1$node.label[1]
        
        root_2 <- sort(mirror_clade)
        root_tree_2 <- Rboretum::getTrimmedTree(tree,root_2)
        root_2_BS <- root_tree_2$node.label[1]
        
        if(!is.na(as.numeric(root_1_BS))){
          bootstrap_list <- c(bootstrap_list,as.numeric(root_1_BS))
        } else if(!is.na(as.numeric(root_2_BS))){
          bootstrap_list <- c(bootstrap_list,as.numeric(root_2_BS))
        } else{ bootstrap_list <- c(bootstrap_list,NA) }
      }
    }
  }
  

  clade_node_df <- data.frame("Clade"=as.character(non_root_clades),"Mirror_Clade"=as.character(non_root_mirror),"Split_Node"=as.integer(node_list),"Split_Bootstrap"=as.numeric(bootstrap_list)) %>%
    filter((!duplicated(Split_Node))) # Two entries for root reduced to one

  return(clade_node_df)
}
