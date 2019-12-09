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

readRooted <- function(to_root,root_taxa,tree_names,dummy_names,prefix,suffix){
  
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
      stop("At least on file from 'to_root' points to an invalid path.")
    } else{
      to_root <- purrr::map(.x=to_root,.f=function(x){file_path_as_absolute(x)}) %>% unlist()
      tree_count <- length(to_root)
      default_name <- purrr::map(to_root,.f = function(x){basename(x)}) %>% unlist() %>% as.character()
      }
  }
 
  # If a single tree path is provided, return a phylo
  if(tree_count == 1){
    tree <- Rboretum::readRooted_Worker(to_root,root_taxa)
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
    
    trees <- purrr::map(.x = to_root,.f = function(x){Rboretum::readRooted_Worker(x,root_taxa)})
    
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