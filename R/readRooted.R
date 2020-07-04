#' Rboretum Rooted Tree Reader
#'
#' This function simultaenously reads in and roots one or more trees at a common root clade
#' @param to_root Where to find tree files. Options include:
#' \itemize{
#'   \item A character vector of one or more tree file paths
#'   \item A path to a single directory containing all tree files
#' }
#' @param root_taxa Character vector containing outgroup species IDs (Must be in tree(s) and monophyletic)
#' @param tree_names OPTIONAL: If multiple tree paths are provided, a character vector of names to assign to trees. Length must equal the number of trees. [Default: Trees will be autonamed based on the filename]
#' @param dummy_names OPTIONAL: If TRUE, and multiple tree paths are provdied, trees will be named with placeholder names (e.g. Tree_1, Tree_2, etc.) [Default: Trees will be autonamed based on the filename]
#' @param prefix OPTIONAL: If 'to_root' is a directory, provide a character vector of file prefixes (e.g. all trees start with "RAxML")
#' @param suffix OPTIONAL: If 'to_root' is a directory, provide a character vector of file suffixes (e.g. all trees end with ".nwk")
#' @param disable_bs OPTIONAL: If TRUE, don't add a mirrored node label to unrooted trees [Default: FALSE; if trees are unrooted, mirror the missing node label]
#' @return A phylo object, rooted at specified taxa, or a named, rooted multiPhlyo
#' @examples 
#' # Read in one tree
#' root_taxa = c('Species_1','Species_2')
#' myTree <- readRooted('/path/to/tree.nwk',root_taxa)
#' 
#' # Read in multiple trees
#' tree_paths <- c('/path/to/tree1.nwk','/path/to/tree2.nwk')
#' tree_names <- c('Tree1','Tree2')
#' myTrees <- readRooted(tree_paths,root_taxa,tree_names=tree_names)
#' 
#' # Read all trees from a directory
#' myTrees <- readRooted('/path/to/tree/dir/',root_taxa) # Trees will be named based off their filenames
#' 
#' # Read all .nwk files from a directory
#' myTrees <- readRooted('/path/to/tree/dir/',root_taxa,suffix=".nwk") # Trees will be named based off their filenames
#' 
#' @export

readRooted_Worker <- function(to_root_worker,root_taxa,disable_bs){
  
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
  
  # Set disable_bs
  if(missing(disable_bs)){
    disable_bs <- FALSE
  } else if(!is.logical(disable_bs)){
    disable_bs <- FALSE
  } else if(length(disable_bs)!=1){
    disable_bs <- FALSE
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
        return(raw_tree)
      }
    }
    
    # Ensure tree can be rooted at root taxa before proceding
    if(has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE)) & has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE))){
      return(NA) # Ape cannot root tree on these taxa
    }
    
    if(root_status == 'rooted'){
      
      # Unroot tree if rooted at a different node
      raw_tree <- ape::unroot.phylo(raw_tree)
      
      if(!has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))){
        return(ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))
      } else if(!has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE))){
        return(ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE))
      } else{ 
        return(NA) # Ape cannot root tree on these taxa, shouldn't be possible
      }
    } else{
      
      # For trees that are unrooted, root tree and if node labels exist, add a node label for the mirror root split
      if(!has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))){
        rooted_tree <- ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE)
      } else if(!has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE))){
        rooted_tree <- ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE)
      } else{ 
        return(NA) # Ape cannot root tree on these taxa, shouldn't be possible
      }
      
      # If no node labels, root and return
      if(!"node.label" %in% attributes(rooted_tree)$names){
        return(rooted_tree)
      }
      
      # If rooted on a single taxon, root and return
      if(length(root_taxa)==1){
        return(rooted_tree)
      }
      
      # If disable_bs is set to TRUE
      if(disable_bs){
        return(rooted_tree)
      }
      
      # Otherwise, root and add a complementary node label for the newly created node
      root_clade <- semiSorter(root_taxa)
      mirror_clade <- semiSorter(mirror_taxa)
      
      subtree_length <- length(ape::subtrees(rooted_tree))
      clade_subtrees <- ape::subtrees(rooted_tree)[2:subtree_length]
      
      clades <- purrr::map(.x=clade_subtrees,.f=function(x){semiSorter(x$tip.label)}) %>% unlist()
      labels <- purrr::map(.x=clade_subtrees,.f=function(x){semiSorter(x$node.label[1])}) %>% unlist()
      
      names(labels) <- clades
      
      # Add mirror root label
      if(labels[[root_clade]] == ""){
        labels[[root_clade]] <- labels[[mirror_clade]]
      } else if(labels[[mirror_clade]] == ""){
        labels[[mirror_clade]] <- labels[[root_clade]]
      } else{
        return(rooted_tree) # Not sure why this would happen...
      }
      
      # Set root label to "Root", and add a node label to the mirror side of the root split
      rooted_tree$node.label <- c("Root",labels)
      return(rooted_tree)
    }
  } else{ # Process RAxML Trees
    
    rax_tree <- treeio::read.raxml(to_root_worker)
    
    # Grab embedded Phylo
    raw_tree <- attributes(rax_tree)$phylo
    
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
        return(raw_tree)
      }
    }
    
    # Ensure tree can be rooted at root taxa before proceding
    if(has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE)) & has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE))){
      return(NA) # Ape cannot root tree on these taxa
    }
    
    if(root_status == 'rooted'){
      
      # Unroot tree if rooted at a different node
      raw_tree <- ape::unroot.phylo(raw_tree)
      
      if(!has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))){
        return(ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))
      } else if(!has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE))){
        return(ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE))
      } else{ 
        return(NA) # Ape cannot root tree on these taxa, shouldn't be possible
      }
    } else{
      
      # For trees that are unrooted, root tree and  add a node label for the mirror root split
      
      if(!has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))){
        rooted_tree <- ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE)
      } else if(!has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE))){
        rooted_tree <- ape::root.phylo(raw_tree,outgroup = mirror_taxa,edgelabel = TRUE,resolve.root = TRUE)
      } else{ 
        return(NA) # Ape cannot root tree on these taxa, shouldn't be possible
      }
      
      # If rooted on a single taxon, root and return
      if(length(root_taxa)==1){
        return(rooted_tree)
      }
      
      # If disable_bs is set to TRUE
      if(disable_bs){
        return(rooted_tree)
      }
      
      # Otherwise, root and add a complementary node label for the newly created node
      root_clade <- semiSorter(root_taxa)
      mirror_clade <- semiSorter(mirror_taxa)
      
      subtree_length <- length(ape::subtrees(rooted_tree))
      clade_subtrees <- ape::subtrees(rooted_tree)[2:subtree_length]
      
      clades <- purrr::map(.x=clade_subtrees,.f=function(x){semiSorter(x$tip.label)}) %>% unlist()
      labels <- purrr::map(.x=clade_subtrees,.f=function(x){semiSorter(x$node.label[1])}) %>% unlist()
      
      names(labels) <- clades
      
      # Add mirror root label
      if(labels[[root_clade]] == ""){
        labels[[root_clade]] <- labels[[mirror_clade]]
      } else if(labels[[mirror_clade]] == ""){
        labels[[mirror_clade]] <- labels[[root_clade]]
      } else{
        return(rooted_tree) # Not sure why this would happen...
      }
      
      # Set root label to "Root", and add a node label to the mirror side of the root split
      rooted_tree$node.label <- c("Root",labels)
      return(rooted_tree)
    }
  }
}

readRooted <- function(to_root,root_taxa,tree_names,dummy_names,prefix,suffix,disable_bs){
  
  # Ensure that a path and root taxa are provided as character vectors
  if(missing(to_root)){
    stop("No tree file or directories indicated with 'to_root'")
  } else if(!is.character(to_root)){
    stop("'to_root' should be a character vector of file paths or the path to a tree directory.")
  } else if(missing(root_taxa)){
    stop("No root taxa provided")
  } else if(!is.character(root_taxa)){
    stop("'root_taxa' should be a character vector of tip labels")
  }
  
  # Create regex search pattern in case a directory is given
  if(missing(prefix)){
    prefix <- c()
  } else if(!is.character(prefix)){
    stop("'prefix' must be a character vector")
  } else{
    prefix <- unlist(purrr::map(.x=prefix,.f=function(x){paste(c("^",x),collapse = '')}))
    prefix <- paste(c("(",paste(prefix,collapse = "|"),")"),collapse = '')
  }
  
  if(missing(suffix)){
    suffix <- c()
  } else if(!is.character(suffix)){
    stop("'suffix' must be a character vector")
  } else{
    suffix <- unlist(purrr::map(.x=suffix,.f=function(x){ifelse(substr(x,start = 1,stop = 1)==".",paste(c("\\",x,"$"),collapse = ''),paste(c(x,"$"),collapse = ''))}))
    suffix <- paste(c("(",paste(suffix,collapse = "|"),")"),collapse = '')
  }
  
  # Set disable_bs
  if(missing(disable_bs)){
    disable_bs <- FALSE
  } else if(!is.logical(disable_bs)){
    disable_bs <- FALSE
  } else if(length(disable_bs)!=1){
    disable_bs <- FALSE
  }
  
  if(length(prefix)==0 & length(suffix)==0){
    tree_regex <- ''
  } else if(length(prefix)>0 & length(suffix)==0){
    tree_regex <- prefix
  } else if(length(prefix)==0 & length(suffix)>0){
    tree_regex <- suffix
  } else if(length(prefix)>0 & length(suffix)>0){
    tree_regex <- paste(paste(c(prefix,"(.*)",suffix),collapse = ""))
  }
  
  # Use dummy names for multiPhylo?
  if(missing(dummy_names)){
    dummy_names <- FALSE
  } else if(!is.logical(dummy_names)){
    dummy_names <- FALSE
  }
  
  # Figure out how many files are being read in
  if(length(to_root)==1){
    
    isFile <- file.exists(to_root) & !dir.exists(to_root)
    isDir <- dir.exists(to_root) & !isFile
    
    if(isFile){ # 'to_root' points to a single valid file
      
      tree_count <- 1
      to_root <- file_path_as_absolute(to_root)

    } else if(isDir){ # 'to_root' points to a valid directory
      
      if(has_error(silent=TRUE,list.files(path=to_root,pattern=tree_regex,full.names = TRUE,include.dirs = FALSE))){
        stop("Can't process file fetch. Check path or regex?")
      } else{
        
        to_root <- list.files(path=to_root,pattern=tree_regex,full.names = TRUE,include.dirs = FALSE)
        
        if(length(to_root)==0){
          stop("Directory found, but no files identified in 'to_root'. Check regex?")
        } else if(length(to_root)==1){
          tree_count <- 1
        } else{
          tree_count <- length(to_root)
          default_name <- purrr::map(to_root,.f = function(x){basename(x)}) %>% unlist() %>% as.character()
        }
      }
    } else{ stop("'to_root' points to neither a valid file or directory.") }
  
  } else{ # 'to_root' is a list of file paths
    
    file_check <- purrr::map(.x = to_root,.f=function(x){ file.exists(x) & !dir.exists(x)}) %>% unlist() %>% all() # Check that all paths in 'to_root' point to valid files
    
    if(!file_check){
      stop("At least one file from 'to_root' points to an invalid path.")
    } else{
      to_root <- purrr::map(.x=to_root,.f=function(x){file_path_as_absolute(x)}) %>% unlist()
      tree_count <- length(to_root)
      default_name <- purrr::map(to_root,.f = function(x){basename(x)}) %>% unlist() %>% as.character()
      }
  }
 
  # If a single tree path is provided, return a phylo
  if(tree_count == 1){
    tree <- Rboretum::readRooted_Worker(to_root,root_taxa,disable_bs)
    if(!Rboretum::isPhylo(tree)){
      stop("'to_root' cannot be rooted with 'root_taxa'")
    } else{
      return(tree)
    }
  } else if(tree_count > 1){ # If multiple tree paths are provided, return a named multiPhylo
    
    if(missing(tree_names)){
      tree_names <- default_name
    } else if(length(tree_names) != tree_count){
      print(paste(c("'tree_names' (",length(tree_names),") and number of trees (",tree_count,") do not match...using default names..."),collapse = ''))
      tree_names <- default_name
    }
    
    trees <- purrr::map(.x = to_root,.f = function(x){Rboretum::readRooted_Worker(x,root_taxa,disable_bs)})
    
    if(any(is.na(unlist(trees)))){
      stop("At least one tree from 'to_root' could not be rooted with 'root_taxa'")
    } else{
      class(trees) <- "multiPhylo"
      names(trees) <- tree_names
      
      if(dummy_names){
        print("Applying dummy tree names as requested...")
        trees <- Rboretum::treeNamer(trees)
      }
      return(trees)
    }
  }
}