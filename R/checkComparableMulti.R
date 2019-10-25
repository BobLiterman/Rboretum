#' Check Which Tree Topologies from Multiphylo Can Be Compared
#'
#' This function takes a multiphylo object and returns a dataframe that indicates which topologies can be compared (>= 3 species shared, with unique topologies). If so, the common species set is also returned.
#' @param trees Multiphylo object
#' @param names Vector of tree identifiers (must have identifier for each tree)
#' @param return_only_comparable OPTIONAL: If TRUE, returns only comparable trees; Anything else, including missing, returns all comparisons (Default: FALSE)
#' @return Dataframe with information about possible comparisons
#' @export
#' @examples
#'
#' trees <- c(tree_1,tree_2,tree_3)
#' names <- c('CDS_Tree','Intron_Tree','UTR_Tree')
#'
#' # Return all comparisons
#' checkComparableMulti(trees,names)
#'
#' # Return only comparable trees
#' checkComparableMulti(trees,names,return_only_comparable=TRUE)
#'

checkComparableMulti <- function(trees,names,return_only_comparable){

  if(missing(return_only_comparable)){
    return_only_comparable <- FALSE
  }
  if(return_only_comparable != TRUE){
    return_only_comparable <- FALSE
  }
  if(missing(names)){
    names <- unlist(lapply(X = 1:length(trees),function(x) paste(c("tree",x),collapse = "_")))
  }
  if(length(trees)!=length(names)){
    print("Multiphylo object and name vector are not the same length. Each tree must be named. Using numbers based on multiphlyo order")
    names <- unlist(lapply(X = 1:length(trees),function(x) paste(c("tree",x),collapse = "_")))
  }
  if(length(trees)<2){
    stop("At least two trees are required for comparison.")
  }
  compareTable <- data.frame("Tree_1"=character(),"Tree_2"=character(),'Comparable'=logical(),'Species'=character())

  names_1 <- c()
  names_2 <- c()
  comparables <- c()
  species_sets <- c()

  tree_count <- length(trees)
  for(i in 1:(tree_count-1)){
    for(j in (i+1):tree_count){
      tree_1 <- trees[[i]]
      names_1 <- c(names_1,names[i])

      tree_2 <- trees[[j]]
      names_2 <- c(names_2,names[j])

      if(has_error(Rboretum::checkComparable(tree_1,tree_2))){
        comparable <- FALSE
      } else{
        comparable <- Rboretum::checkComparable(tree_1,tree_2)
      }
      comparables <- c(comparables,comparable)

      if(comparable){
        species_sets <- c(species_sets,paste(sort(Rboretum::getSharedSpecies(tree_1,tree_2)),collapse = ";"))
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
