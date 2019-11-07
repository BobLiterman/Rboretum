#' Rboretum Multiphylo Comparability Checker
#'
#' This function takes a multiphylo object and returns a dataframe that indicates which topologies can be compared (>= 3 species shared, with unique topologies). If so, the common species set is also returned.
#' @param trees Multiphylo object
#' @param tree_names OPTIONAL: Vector of tree names
#' @param return_only_comparable OPTIONAL: If TRUE, returns only comparable trees; (Default: FALSE)
#' @return Dataframe with information about possible comparisons
#' @export
#' @examples
#'
#' trees <- c(tree_1,tree_2,tree_3)
#' names <- c('CDS_Tree','Intron_Tree','UTR_Tree')
#'
#' # Return all comparisons
#' canCompare(trees,names)
#'
#' # Return only comparable trees
#' canCompare(trees,names,return_only_comparable=TRUE)
#'

canCompare <- function(trees,tree_names,return_only_comparable){
  if(missing(tree_names)){
    # Check that input is multiphylo and has at least 2 trees
    if(has_error(unlist(attributes(trees)$class))){ 
      stop("'trees' argument should be a multiPhylo object")
    } else if(!"multiPhylo" %in% unlist(attributes(trees)$class)){
      stop("'trees' argument should be a multiPhylo object")
    } else if(length(trees)<2){
      stop("At least two trees are required for comparison. For a single tree, use checkTreeTaxa()")
    }
    
    tree_count <- length(trees)
    tree_names <- unlist(lapply(X = 1:tree_count,function(x) paste(c("Tree",x),collapse = "_")))
    
  } else{
    if(has_error(unlist(attributes(trees)$class))){ 
      stop("'trees' argument should be a multiPhylo object")
    } else if(!"multiPhylo" %in% unlist(attributes(trees)$class)){
      stop("'trees' argument should be a multiPhylo object")
    } else if(length(trees)<2){
      stop("At least two trees are required for comparison. For a single tree, use checkTreeTaxa()")
    }
    
    tree_count <- length(trees)
    
    if(tree_count!=length(tree_names)){
      
      print("Multiphylo object and name vector are not the same length. Each tree must be named. Using numbers based on multiphlyo order")
      tree_names <- unlist(lapply(X = 1:tree_count,function(x) paste(c("Tree",x),collapse = "_")))
      
    } else{
      tree_names <- as.character(tree_names)
    }
  }
  
  if(missing(return_only_comparable)){
    return_only_comparable <- FALSE
  } else{ 
    if(!is.logical(return_only_comparable)){
      return_only_comparable <- FALSE
    }
  }
  
  names_1 <- c()
  names_2 <- c()
  comparables <- c()
  species_sets <- c()

  for(i in 1:(tree_count-1)){
    for(j in (i+1):tree_count){
      tree_1 <- trees[[i]]
      names_1 <- c(names_1,tree_names[i])

      tree_2 <- trees[[j]]
      names_2 <- c(names_2,tree_names[j])
    
      if(has_error(Rboretum::sameTopology(tree_1,tree_2))){
        comparable <- FALSE
      } else{
        if(!Rboretum::sameTopology(tree_1,tree_2)){
        comparable <- TRUE
      } else{ comparable <- FALSE }
      }
      
      comparables <- c(comparables,comparable)

      if(comparable){
        species_sets <- c(species_sets,paste(sort(Rboretum::getSharedSpecies(c(tree_1,tree_2))),collapse = ";"))
      } else{
        species_sets <- c(species_sets,NA)
      }
    }
  }

  compareTable <- data.frame("Tree_1"=as.character(names_1),"Tree_2"=as.character(names_2),'Comparable'=as.logical(comparables),'Species'=as.character(species_sets))

  if(return_only_comparable==FALSE){
    return(compareTable)
  }
  else{
      compareTable <- compareTable %>% filter(Comparable == TRUE)
      if(nrow(compareTable )>0){
        return(compareTable)
      } else{
          stop("No trees are comparable")
      }
    }
}
