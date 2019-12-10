#' Rboretum Alignment Signal Checker
#'
#' This function returns TRUE if the passed object is the output of getAlignmentSignal; otherwise, FALSE
#' @param test_object R object to check
#' @param tree OPTIONAL: Check that passed support object contains information about a specific set of trees. Options include:
#' \itemize{
#'   \item A rooted phylo object; or,
#'   \item A named, rooted multiPhylo object where all trees share 3+ taxa
#' }
#' @return TRUE if output of getTreeSupport(); otherwise, FALSE
#' @export

isTreeSupport <- function(test_object,tree){
  
  # Ensure dataframe columns match expected and that data exists
  if(!is.data.frame(test_object) & !is.list(test_object)){
    return(FALSE)
  } else if(!Rboretum::isPhylo(tree,check_rooted = TRUE) & !Rboretum::isMultiPhylo(tree,check_rooted = TRUE,check_named = TRUE,check_three_taxa = TRUE)){
    stop("'tree' is not a valid rooted phylo or named, rooted multiPhylo object")
  }
  
  # Get tree clades to test
  if(Rboretum::isPhylo(tree)){
    tree_clades <- Rboretum::getTreeClades(tree)
    tree_count <- 1
  } else if(Rboretum::isMultiPhylo(tree,check_all_equal = TRUE)){
    tree <- as.phylo(tree[[1]])
    tree_clades <- Rboretum::getTreeClades(tree)
    tree_count <- 1
  } else if(Rboretum::isMultiPhylo(tree,check_all_unique = TRUE)){
    tree_names <- names(tree)
    tree_clades <- purrr::map(.x=tree,.f=function(x){Rboretum::getTreeClades(x)})
    names(tree_clades) <- tree_names
    tree_count <- length(tree)
  } else if(Rboretum::isMultiPhylo(tree,check_some_unique = TRUE)){
    tree <- Rboretum::getUniqueTopologies(tree)
    tree_names <- names(tree)
    tree_clades <- purrr::map(.x=tree,.f=function(x){Rboretum::getTreeClades(x)})
    names(tree_clades) <- tree_names
    tree_count <- length(tree)
  }
  
  if(is.data.frame(test_object)){
    if(tree_count != 1){ # Single dataframe support must be for one tree
      return(FALSE)
    } else if(has_error(silent=TRUE,expr=!all(names(test_object)[1:3] == c("Clade","Mirror_Clade","Split_Node")))){
      return(FALSE)
    } else if(!all(names(test_object)[1:3] == c("Clade","Mirror_Clade","Split_Node"))){ # Test that cols 1:3 are split cols
      return(FALSE)
    } else if(nrow(test_object) != length(tree_clades)){ # Ensure data
      return(FALSE)
    } else if(any(is.na(test_object$Split_Node))){ # Tree support lacks root information
      return(FALSE)
    } else if(has_error(silent = TRUE,expr = !all(sort(tree_clades)==sort(as.character(test_object$Clade))))){
      return(FALSE)
    } else if(!all(sort(tree_clades)==sort(as.character(test_object$Clade)))){ # Ensure clades match between tree and support
      return(FALSE)
    } else{
      return(TRUE)
    }
  }
  
  if(is.list(test_object)){
    
    # Ensure results list is the same length as trees
    if(length(test_object) != tree_count){ 
      return(FALSE)
    } 
    
    # Ensure support list has the same name as trees
    if(!all(names(test_object) == tree_names)){ 
      return(FALSE)
    }
      
    # Test that cols 1:3 are split cols
    err_name <- purrr::map(.x=names(test_object),.f=function(x){has_error(silent=TRUE,expr=all(names(test_object[[x]])[1:3] == c("Clade","Mirror_Clade","Split_Node")))}) %>% unlist()
    if(any(err_name)){
      return(FALSE)
    } 
    
    name_check <- purrr::map(.x=names(test_object),.f=function(x){all(names(test_object[[x]])[1:3] == c("Clade","Mirror_Clade","Split_Node"))}) %>% unlist()
    if(!all(name_check)){
      return(FALSE)
    }
    
    # Ensure data
    data_check <- purrr::map(.x=names(test_object),.f=function(x){nrow(test_object[[x]]) == length(Rboretum::getTreeClades(tree[[x]]))}) %>% unlist()
    if(!all(data_check)){
      return(FALSE)
    }
    
    # Tree support lacks root information
    root_check <- purrr::map(.x=names(test_object),.f=function(x){any(is.na(test_object[[x]]$Split_Node))}) %>% unlist()
    if(any(root_check)){
      return(FALSE)
    }
    
    # Ensure clades match between tree and support
    err_clade <- purrr::map(.x=names(test_object),.f=function(x){has_error(silent=TRUE,expr=all(sort(test_object[[x]]$Clade) == Rboretum::getTreeClades(tree[[x]])))}) %>% unlist()
    if(any(err_name)){
      return(FALSE)
    } 
    
    clade_check <- purrr::map(.x=names(test_object),.f=function(x){all(sort(test_object[[x]]$Clade) == Rboretum::getTreeClades(tree[[x]]))}) %>% unlist()
    if(!all(name_check)){
      return(FALSE)
    }

    # Ensure each tree has data for all alignments
    col_nums <- purrr::map(.x=names(test_object),.f=function(x){ncol(test_object[[x]])}) %>% unlist() %>% as.integer()
    if(length(unique(col_nums))!=1){
      return(FALSE)
    }
    
    # Ensure each tree has data for all alignments
    col_names <- purrr::map(.x=names(test_object),.f=function(x){paste(names(test_object[[x]]),collapse = ";")}) %>% unlist() %>% as.character()
    if(length(unique(col_names))!=1){
      return(FALSE)
    }
    else{
      return(TRUE)
    }
  }
}