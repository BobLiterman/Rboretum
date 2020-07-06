#' Rboretum Rooted Tree Reader
#'
#' This function simultaenously reads in and roots one or more trees at a common root clade
#' @param to_root Where to find tree files. Options include:
#' \itemize{
#'   \item A character vector of one or more tree file paths
#'   \item A path to a single directory containing all tree files
#' }
#' @param root_taxa Character vector containing outgroup species IDs (Must be in tree(s) and monophyletic)
#' @return A phylo object, rooted at specified taxa
#' @export

readRooted_Worker <- function(to_root_worker,root_taxa){
  
  # Ensure that a path and root taxa are provided as character vectors
  if(missing(to_root_worker)){
    return(NA) # No tree file or directories indicated with 'to_root_worker'
  } else if(!is.character(to_root_worker)){
    return(NA) # 'to_root_worker' should be a character path to a tree file.
  } else if(!file.exists(to_root_worker)){ # 'to_root_worker' must exist
    return(NA)
  } else if(dir.exists(to_root_worker)){ # 'to_root_worker' should not be a directory
    return(NA)
  } else if(missing(root_taxa)){
    return(NA) # No root taxa provided
  } else if(!is.character(root_taxa)){
    return(NA) # 'root_taxa' should be a character vector of tip labels
  }

  # If file exists and root taxa are provided, assess tree type and root status
  if(!has_error(silent=TRUE,expr=treeio::read.raxml(to_root_worker))){ # Check if tree can be read with treeio::read.raxml. Trees without RAxML-like branch labels will have an error
    tree_type <- 'raxml'
    if(ape::is.rooted(attributes(treeio::read.raxml(to_root_worker))$phylo)){
      root_status <- 'rooted'
    } else{
      root_status <- 'unrooted'
    }
  } else if(!has_error(silent=TRUE,expr=ape::read.nexus(to_root_worker))){ # Check if tree can be read with ape::read.nexus. Non-nexus trees will have an error
    tree_type <- 'non_raxml'
    tree_subtype <- 'nexus'
    if(ape::is.rooted(ape::read.nexus(to_root_worker))){
      root_status <- 'rooted'
    } else{
      root_status <- 'unrooted'
    }
  } else if(!has_error(silent=TRUE,expr=ape::read.tree(to_root_worker))){ # Check if tree can be read with ape::read.tree
    tree_type <- 'non_raxml'
    tree_subtype <- 'typical'
    if(ape::is.rooted(ape::read.tree(to_root_worker))){
      root_status <- 'rooted'
    } else{
      root_status <- 'unrooted'
    }
  } else{ 
    return(NA) # 'to_root_worker' does not point to file that can be read in by treeio::read.raxml, ape::read.tree() or ape::read.nexus()
  }
  
  # Process non-raxml trees
  if(tree_type == 'non_raxml'){
    
    # Read in tree
    if(tree_subtype == 'nexus'){
      raw_tree <- ape::read.nexus(to_root_worker)
    } else{
      raw_tree <- ape::read.tree(to_root_worker)
    }
    
    # Get tree species
    tree_species <- raw_tree$tip.label
    
    # Ensure at least one root taxon is present in the tree
    if(!any(root_taxa %in% tree_species)){
      return(NA) # 'root_taxa' completely absent from tree
    }
    
    # Root at all avaiable taxa
    root_taxa <- root_taxa[root_taxa %in% tree_species]
    
    # Mirror taxa is all taxa not in 'root_taxa'
    mirror_taxa <- tree_species[!tree_species %in% root_taxa]
    
    # Ensure root_taxa are monophyletic
    if(!Rboretum::checkTips(ape::unroot.phylo(raw_tree),root_taxa,check_mono=TRUE)){
      return(NA)
    }
    
    # If tree is already rooted at root_taxa, return unchanged
    if(root_status == 'rooted'){
      if(Rboretum::checkTips(raw_tree,root_taxa,check_root=TRUE)){
        raw_tree$node.label[[1]] <- "Root"
        return(raw_tree)
      }
    }
    
    # Ensure tree can be rooted at root taxa before proceding
    if(has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE)) & has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE))){
      return(NA) # Ape cannot root tree on these taxa
    }
    
    # Unroot rooted trees
    if(root_status == 'rooted'){
      raw_tree <- ape::unroot.phylo(raw_tree)
    }
    
    if(!has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))){
      rooted_tree <- ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE)
      rooted_tree$node.label[[1]] <- "Root"
      return(rooted_tree)
    } else if(!has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE))){
      rooted_tree <- ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE)
      rooted_tree$node.label[[1]] <- "Root"
      return(rooted_tree)
    } else{ 
      return(NA) # Ape cannot root tree on these taxa, shouldn't be possible
    }
    
  } else{ # Process RAxML Trees
    
    rax_tree <- treeio::read.raxml(to_root_worker)
    
    # Grab embedded Phylo
    raw_tree <- attributes(rax_tree)$phylo
    
    # Get tree species
    tree_species <- raw_tree$tip.label
    
    # Ensure at least one root taxon is present in the tree
    if(!any(root_taxa %in% tree_species)){
      return(NA) # 'root_taxa' completely absent from tree
    }
    
    # Root at all avaiable taxa
    root_taxa <- root_taxa[root_taxa %in% tree_species]
    mirror_taxa <- tree_species[!tree_species %in% root_taxa]
    
    root_clade <- semiSorter(root_taxa)
    mirror_clade <- semiSorter(mirror_taxa)
    
    # Ensure root_taxa are monophyletic
    if(!Rboretum::checkTips(ape::unroot.phylo(raw_tree),root_taxa,check_mono=TRUE)){
      return(NA)
    }
    
    # Get subtree order from raw_tree
    raw_subtree <- subtrees(raw_tree)
    raw_subtree_length <- length(raw_subtree)
    
    raw_clades <- purrr::map(.x=raw_subtree,.f=function(x){semiSorter(x$tip.label)}) %>% unlist()
    
    # Add taxa information to bootstrap_tibble
    bootstrap_tibble <- attributes(rax_tree)$data
    node_taxa_tibble <- tibble(node=integer(),taxa=character())
    
    for(i in 1:length(subtrees(raw_tree))){
      subtree <- subtrees(raw_tree)[[i]]
      node_taxa_tibble <- node_taxa_tibble %>% add_row(node=as.integer(subtree$node.label[1]),taxa=vectorSemi(subtree$tip.label) %>% semiSorter())
    }
    
    bootstrap_tibble <- left_join(bootstrap_tibble,node_taxa_tibble,by='node')
    bootstrap_names <- pull(bootstrap_tibble,taxa)
    named_bootstraps <- pull(bootstrap_tibble,bootstrap) %>% `names<-`(bootstrap_names)
    
    # Add node labels to phylo object
    raw_tree$node.label <- purrr::map(.x=raw_clades,.f=function(x){named_bootstraps[x]}) %>% unlist()
    
    # If tree is already rooted at root_taxa, return unchanged
    if(root_status == 'rooted'){
      if(Rboretum::checkTips(raw_tree,root_taxa,check_root=TRUE)){
        raw_tree$node.label[[1]] <- "Root"
        return(raw_tree)
      }
    }
    
    # Unroot rooted trees
    if(root_status == 'rooted'){
      raw_tree <- ape::unroot.phylo(raw_tree)
    }
    
    # Ensure tree can be rooted at root taxa before proceding
    if(has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE)) & has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE))){
      return(NA) # Ape cannot root tree on these taxa
    }
      
    if(!has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))){
      rooted_tree <- ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE)
      rooted_tree$node.label[[1]] <- "Root"
      return(rooted_tree)
    } else if(!has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE))){
      rooted_tree <- ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE)
      rooted_tree$node.label[[1]] <- "Root"
      return(rooted_tree)
    } else{ 
      return(NA) # Ape cannot root tree on these taxa, shouldn't be possible
    }
  }
}