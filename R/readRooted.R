#' Rboretum Rooted Tree Reader
#'
#' This function simultaenously reads in and roots one or more trees at a common root clade
#' @param to_root Where to find tree files. Options include:
#' \itemize{
#'   \item A character vector of one or more tree file paths
#'   \item A path to a single directory containing all tree files
#' }
#' @param root_taxa Character vector containing outgroup species IDs (Must be in tree(s) and monophyletic)
#' @param tree_regex OPTIONAL: If 'to_root' is a directory, provide a regex code to get all tree files. [See examples below; Default: Read all files in directory]
#' @param tree_names OPTIONAL: If multiple tree paths are provided, a character vector of names to assign to trees. Length must equal the number of trees. [Default: Trees will be autonamed based on the filename]
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
#' myTrees <- readRooted('/path/to/tree/dir/',root_taxa,tree_regex="\\.nwk$") # Trees will be named based off their filenames
#' 
#' # tree_regex examples
#' tree_regex = "^PREFIX" # Files that start with PREFIX
#' tree_regex = "\\.nwk$" # All files ending in .nwk
#' tree_regex = "\\.nwk$|\\.nex$" # All files ending in .nwk or .nex
#' @export

readRooted <- function(to_root,root_taxa,tree_regex,tree_names){
  
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
  
  if(length(to_root)==1){
    isFile <- file.exists(to_root) & !dir.exists(to_root)
    isDir <- dir.exists(to_root) & !isFile
    
    if(!isDir & !isFile){
      stop("'to_root' points to neither a valid file or directory.")
    } else{
      if(isFile){ # If 'to_root' points to a valid file...
        tree <- Rboretum::readRooted_Worker(to_root,root_taxa)
        if(!Rboretum::isPhylo(tree)){
          stop("'to_root' cannot be rooted with 'root_taxa'")
        } else{
          return(tree)
        }
        } else{ # If 'to_root' points to a valid directory...
          if(missing(tree_regex)){
            tree_paths <- list.files(path=to_root,full.names = TRUE,include.dirs = FALSE)
            if(length(tree_paths)==0){
              stop("No files found in 'to_root'")
            } else if(length(tree_paths)==1){
              tree <- Rboretum::readRooted_Worker(tree_paths,root_taxa)
              if(!Rboretum::isPhylo(tree)){
                stop("'to_root' cannot be rooted with 'root_taxa'")
              } else{
                return(tree)
              } 
            } else{ #If more than one file is returned...
              tree <- purrr::map(.x = tree_paths,.f = function(x){Rboretum::readRooted_Worker(x,root_taxa)})
              if(any(is.na(unlist(tree)))){
                stop("At least one tree from 'to_root' could not be rooted with 'root_taxa'")
              } else{
                class(tree) <- "multiPhylo"
                default_tree_names <- list.files(path=to_root,full.names = TRUE,include.dirs = FALSE) %>% basename()
              }
              if(missing(tree_names)){
                names(tree) <- default_tree_names
                return(tree)
              } else if(length(tree_names) != length(tree)){
                print("Not enough tree names supplied, using default naming based on filenames...")
                names(tree)  <- default_tree_names
                return(tree)
              }
              else{
                names(tree) <- tree_names
                return(tree)
              }
            }
            
          } else{
            if(has_error(silent=TRUE,list.files(path=to_root,pattern=tree_regex,full.names = TRUE,include.dirs = FALSE))){
              stop("Can't process file fetch. Check regex?")
            } else{
              tree_paths <- list.files(path=to_root,full.names = TRUE,include.dirs = FALSE)
              if(length(tree_paths)==0){
                stop("No files found in 'to_root'")
              } else if(length(tree_paths)==1){
                tree <- Rboretum::readRooted_Worker(tree_paths,root_taxa)
                if(!Rboretum::isPhylo(tree)){
                  stop("'to_root' cannot be rooted with 'root_taxa'")
                } else{
                  return(tree)
                } 
              } else{ #If more than one file is returned...
                tree <- purrr::map(.x = tree_paths,.f = function(x){Rboretum::readRooted_Worker(x,root_taxa)})
                if(any(is.na(unlist(tree)))){
                  stop("At least one tree from 'to_root' could not be rooted with 'root_taxa'")
                } else{
                  class(tree) <- "multiPhylo"
                  default_tree_names <- list.files(path=to_root,pattern=tree_regex,full.names = TRUE,include.dirs = FALSE) %>% basename()
                }
                if(missing(tree_names)){
                  names(tree) <- default_tree_names
                  return(tree)
                } else if(length(tree_names) != length(tree)){
                  print("Not enough tree names supplied, using default naming based on filenames...")
                  names(tree)  <- default_tree_names
                  return(tree)
                }
                else{
                  names(tree) <- tree_names
                  return(tree)
                }
              } 
            }
          }
      }
    }
  } else{ # If 'to_root' is a list of file paths...
    file_check <- purrr::map(.x = to_root,.f=function(x){ file.exists(to_root) & !dir.exists(to_root)}) %>% unlist() # Check that all paths in 'to_root' point to valid files
    if(!all(file_check)){
      stop("At least on file from 'to_root' points to an invalid path.")
    } else{ # If all paths are to valid files...
      tree <- purrr::map(.x = to_root,.f = function(x){Rboretum::readRooted_Worker(x,root_taxa)})
      if(any(is.na(unlist(tree)))){
        stop("At least one tree from 'to_root' could not be rooted with 'root_taxa'")
      } else{
        class(tree) <- "multiPhylo"
        default_tree_names <- list.files(path=to_root,full.names = TRUE,include.dirs = FALSE) %>% basename()
      }
      if(missing(tree_names)){
        names(tree) <- default_tree_names
        return(tree)
      } else if(length(tree_names) != length(tree)){
        print("Not enough tree names supplied, using default naming based on filenames...")
        names(tree)  <- default_tree_names
        return(tree)
      }
      else{
        names(tree) <- tree_names
        return(tree)
      }
    }
  }
}