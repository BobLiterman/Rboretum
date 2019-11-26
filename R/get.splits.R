#' Rboretum Tree Splitter
#'
#' This function takes a rooted phylo object and returns information about each monophyletic group/split
#' @param tree Rooted phylo object
#' @return Four-column dataframe; 1: Semi-colon separated monophyletic clade; 2: 'Mirror Clade'; 3: Phylo Node ID; 4: Node Boostrap (if present)
#' @export
#' @examples
#' get.splits(tree)
#'

get.splits <- function(tree){
  
  if(has_error(ape::is.rooted(tree))){
    stop("Error in ape::is.rooted. Is 'tree' a phylo object?")
  } else if(!ape::is.rooted(tree)){
    stop("Tree must be rooted for get.splits")}
  
  # Get species
  tree_species <- sort(tree$tip.label)
  species_count <- length(tree_species)
  
  mono_clades <- c()
  mirror_clades <- c()
  node_list <- c()
  bootstrap_list <- c()
  
  if(is.null(tree$node.label)){
    hasBS <- FALSE
  } else { hasBS <- TRUE }
  
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
        mono_clades <- c(mono_clades,sort(temp_clade) %>% paste(collapse = ";"))
        mirror_clades <- c(mirror_clades,sort(mirror_clade) %>% paste(collapse = ";"))
        node_list <- c(node_list,ape::getMRCA(tree,temp_clade))
        
        if(hasBS){
          bs_tree <- Rboretum::trim.tree(tree,temp_clade)
          node_bs <- bs_tree$node.label[1]
          if(!is.na(as.numeric(node_bs))){
            bootstrap_list <- c(bootstrap_list,as.numeric(node_bs))
          } else{ bootstrap_list <- c(bootstrap_list,NA) }
        }
      } else if(mono_B & !(mono_A)){
        mono_clades <- c(mono_clades,sort(mirror_clade) %>% paste(collapse = ";"))
        mirror_clades <- c(mirror_clades,sort(temp_clade) %>% paste(collapse = ";"))
        node_list <- c(node_list,ape::getMRCA(tree,mirror_clade))
        
        if(hasBS){
          bs_tree <- Rboretum::trim.tree(tree,temp_clade)
          node_bs <- bs_tree$node.label[1]
          if(!is.na(as.numeric(node_bs))){
            bootstrap_list <- c(bootstrap_list,as.numeric(node_bs))
          } else{ bootstrap_list <- c(bootstrap_list,NA) }
        }
      } else if(mono_A & mono_B){
        
        mono_clades <- c(mono_clades,sort(temp_clade) %>% paste(collapse = ";"))
        mirror_clades <- c(mirror_clades,sort(mirror_clade) %>% paste(collapse = ";"))
        node_list <- c(node_list,NA)
        
        if(temp_length > 1 & mirror_length > 1){
        
            if(hasBS){
            root_1 <- sort(temp_clade)
            root_tree_1 <- Rboretum::trim.tree(tree,root_1)
            root_1_BS <- root_tree_1$node.label[1]
            
            root_2 <- sort(mirror_clade)
            root_tree_2 <- Rboretum::trim.tree(tree,root_2)
            root_2_BS <- root_tree_2$node.label[1]
            
            if(!is.na(as.numeric(root_1_BS))){
              bootstrap_list <- c(bootstrap_list,as.numeric(root_1_BS))
            } else if(!is.na(as.numeric(root_2_BS))){
              bootstrap_list <- c(bootstrap_list,as.numeric(root_2_BS))
            } else{ bootstrap_list <- c(bootstrap_list,NA) }
          }
        } else{
          
          if(hasBS){
            bootstrap_list <- c(bootstrap_list,NA)
          }
        }
      }
    }
  }
  
  if(hasBS){
    clade_node_df <- data.frame("Clade"=as.character(mono_clades),"Mirror_Clade"=as.character(mirror_clades),"Split_Node"=as.integer(node_list),"Split_Bootstrap"=as.numeric(bootstrap_list)) %>%
      filter((!duplicated(Split_Node))) # Two entries for root reduced to one
  } else{
    clade_node_df <- data.frame("Clade"=as.character(mono_clades),"Mirror_Clade"=as.character(mirror_clades),"Split_Node"=as.integer(node_list)) %>%
      mutate("Split_Bootstrap" = NA) %>%
      filter((!duplicated(Split_Node))) # Two entries for root reduced to one
  }
  
  return(clade_node_df)
}