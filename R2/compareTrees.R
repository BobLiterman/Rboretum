#' Rboretum Tree Comparability Checker
#'
#' This function takes a multiPhylo object and returns a dataframe that indicates which topologies can be compared (>= 3 species shared, with unique topologies). If so, the common species set is also returned.
#' @param trees Named multiPhylo object [Set tree names by using names(trees) <- c('Tree1_Name','Tree2_Name',etc.)]
#' @return Dataframe with information about possible comparisons
#' @export

compareTrees <- function(trees){
  
  if(!Rboretum::isMultiPhylo(trees)){
    stop("'trees' does not appear to be a valid multiPhylo object with 2+ trees")
  } else if(is.null(names(trees))){
    stop("Trees in multiPhylo must be named for compareTrees. Name trees via names(trees) <- c('Tree1_Name','Tree2_Name',etc.)")
  }
  
  tree_names <- names(trees)
  tree_count <- length(trees)
  
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
      
      comparable <- Rboretum::checkComparable(tree_1,tree_2)
      
      comparables <- c(comparables,comparable)

      if(comparable){
        species_sets <- c(species_sets,paste(sort(Rboretum::getSharedTaxa(c(tree_1,tree_2))),collapse = ";"))
      } else{
        species_sets <- c(species_sets,NA) 
      }
    }
  }

  compareTable <- data.frame("Tree_1"=as.character(names_1),"Tree_2"=as.character(names_2),'Comparable'=as.logical(comparables),'Species'=as.character(species_sets))
  return(compareTable)
}
