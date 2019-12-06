#' Rboretum MultiPhylo Checker
#'
#' This function returns TRUE if the passed object is of class multiPhylo and has 2+ trees; Otherwise, FALSE
#' @param test_object R object to check
#' @param check_named OPTIONAL: If TRUE, also check if trees are named. [Default: FALSE, don't check names]
#' @param check_rooted OPTIONAL: If TRUE, also check if all trees are rooted. [Default: FALSE, don't check rootedness]
#' @param check_three_taxa OPTIONAL: If TRUE, also check if all trees share at least three taxa. [Default: FALSE, don't check taxa]
#' @param check_all_taxa OPTIONAL: If TRUE, also check if all trees share all taxa. [Default: FALSE, don't check taxa]
#' @param check_all_equal OPTIONAL: If TRUE, also check whether all trees have one, shared topology after pruning to a common set of taxa [Default: FALSE, don't check topology]
#' @param check_some_equal OPTIONAL: If TRUE, also check whether some, but not all, trees share a topology after pruning to a common set of taxa [Default: FALSE, don't check topology]
#' @param check_all_unique OPTIONAL: If TRUE, also check whether all trees have a unique topology after pruning to a common set of taxa [Default: FALSE, don't check topology]
#' @return TRUE if 'test_object' is a multiPhylo with 2+ trees; otherwise, FALSE
#' @examples 
#' isMultiPhylo(myTrees) # Check if 'myTrees' is a valid multiPhylo object with 2+ trees
#' isMultiPhylo(myTrees,check_named=TRUE) # Check if 'myTrees' is a valid multiPhylo object with 2+ trees, and that all trees have names
#' isMultiPhylo(myTrees,check_rooted=TRUE) # Check if 'myTrees' is a valid multiPhylo object with 2+ trees, and that all trees are rooted
#' isMultiPhylo(myTrees,check_three_taxa) # Check if 'myTrees' is a valid multiPhylo object with 2+ trees, and that all trees share at least three taxa
#' isMultiPhylo(myTrees,check_all_taxa) # Check if 'myTrees' is a valid multiPhylo object with 2+ trees, and that all trees share all taxa
#' @export

isMultiPhylo <- function(test_object,check_named,check_rooted,check_three_taxa,check_all_taxa,check_all_equal,check_some_equal,check_all_unique){

  # Check if trees are named?
  if(missing(check_named)){
    check_named <- FALSE
  } else if(!is.logical(check_named)){
    check_named <- FALSE
  }
  
  # Check if all trees are rooted?
  if(missing(check_rooted)){
    check_rooted <- FALSE
  } else if(!is.logical(check_rooted)){
    check_rooted <- FALSE
  }
  
  # Check if all trees share at least three taxa?
  if(missing(check_three_taxa)){
    check_three_taxa <- FALSE
  } else if(!is.logical(check_three_taxa)){
    check_three_taxa <- FALSE
  }
  
  # Check if all trees share all taxa?
  if(missing(check_all_taxa)){
    check_all_taxa <- FALSE
  } else if(!is.logical(check_all_taxa)){
    check_all_taxa <- FALSE
  }
  
  # Check if all trees share one topology?
  if(missing(check_all_equal)){
    check_all_equal <- FALSE
  } else if(!is.logical(check_all_equal)){
    check_all_equal <- FALSE
  }
  
  # Check if all trees have a unique topology?
  if(missing(check_all_unique)){
    check_all_unique <- FALSE
  } else if(!is.logical(check_all_unique)){
    check_all_unique <- FALSE
  }

  # Check if some, but not all trees share a topology?
  if(missing(check_some_equal)){
    check_some_equal <- FALSE
  } else if(!is.logical(check_some_equal)){
    check_some_equal <- FALSE
  }
  
  if(check_all_equal & check_all_unique){
    stop("Cannot 'check_all_equal' and 'check_all_unique' in the same run.")
  } else if(check_all_equal & check_some_equal){
    stop("Cannot 'check_all_equal' and 'check_some_equal' in the same run.")
  } else if(check_all_unique & check_some_equal){
    stop("Cannot 'check_all_unique' and 'check_some_equal' in the same run.")
  }

  if(has_error(silent=TRUE,expr=unlist(attributes(test_object)))){ # Can attributes be unlisted?
    print("Can't access object attributes...")
    return(FALSE) # Object attributes can't be unlisted --> FALSE
  } else if(!'multiPhylo' %in% unlist(attributes(test_object)$class)){
      return(FALSE) # 'multiPhylo' not in object class --> FALSE
    } else if(length(test_object) < 2){
      print("'multiPhylo' contains fewer than two trees...")
      return(FALSE) # 'test_object' has length < 2 --> FALSE
    } else if(!purrr::map(.x = test_object,.f = function(x){Rboretum::isPhylo(x)}) %>% unlist() %>% all()){
      print("'multiPhylo' contains non-phylo objects...")
      return(FALSE) # 'test_object' not made up of phylo objects --> FALSE
    } else{
      
      # Summarize taxa in multiPhylo 'test_object'
      tree_taxa <- purrr::map(.x = test_object,.f = function(x){x$tip.label}) %>% unlist() %>% unique() %>% sort() # Get all unique tip labels among 'trees'
      taxa_count <- length(tree_taxa)
      tree_count <- length(test_object)
      
      # Get shared species
      tip_table  <- purrr::map(.x=test_object,.f=function(x){x$tip.label}) %>% unlist() %>% table() # Tally tip counts
      shared_species <- tree_taxa[tip_table==tree_count] # Find tips that occur in all 'trees'
      shared_count <- length(shared_species)
      
      has_names <- !is.null(names(test_object))
      
      name_error <- FALSE
      
      if(has_names){
        tree_names <- names(test_object)
        name_length <- length(names(test_object))
        name_error <- any(is.na(tree_names)) | any(tree_names == "")| any(is.null(tree_names)) | name_length != tree_count | any(duplicated(tree_names))
      }
      
      if(check_rooted & !all(unlist(purrr::map(.x = test_object,.f = ape::is.rooted)))){
        print('At least one tree is unrooted...')
        return(FALSE) # 'test_object' is a valid multiPhylo but trees are not rooted --> FALSE
      } else if(check_three_taxa & shared_count < 3){
        print('Trees in multiPhylo share fewer than three taxa...')
        return(FALSE) # 'test_object' is a valid multiPhylo but trees share fewer than three taxa --> FALSE
      } else if(check_all_taxa & shared_count != taxa_count){
        print('Trees in multiPhylo share do not have identical taxa...')
        return(FALSE) # 'test_object' is a valid multiPhylo but trees do not have identical taxa --> FALSE
      } else if(check_named & (!has_names | name_error)){
        print('Everything about the multiPhylo is fine, but the trees are not named, or are not named stably. You can add placeholder tree names by running Rboretum::treeNamer(trees)')
        return(FALSE) # 'test_object' is a valid multiPhylo but trees are not named --> FALSE
      } else if(check_all_equal | check_all_unique | check_some_equal){
        
        if(shared_count < 3){
          stop("Cannot test topologies when trees share fewer than three taxa.")
        }
        
        # Prune 'test_object' down to common set of taxa
        pruned_tree <- purrr::map(.x = test_object,.f = function(x){ape::drop.tip(x,x$tip.label[-match(shared_species, x$tip.label)])})
        class(pruned_tree) <- "multiPhylo"
      
        # Compare all tree topologies
        top_check <- c()
        
        for(i in 1:(tree_count-1)){
          for(j in (i+1):tree_count){
            top_check <- c(top_check,ape::all.equal.phylo(pruned_tree[[i]],pruned_tree[[j]],use.edge.length = FALSE))
          }
        }
        
        if(check_some_equal){
          if(!all(top_check) & any(top_check)){ # Not all trees share the same topology, but some do --> TRUE
            return(TRUE)
          } else{
            print("Trees either all share the same topolgy, or all trees are unique.")
            return(FALSE)
          }
        } else if(check_all_equal){
          if(all(top_check)){ # All trees share the same topology --> TRUE
            return(TRUE)
          } else{
            print("Trees do not all share a single topology.")
            return(FALSE)
          }
        } else if(check_all_unique){
          if(!any(top_check)){ # No trees share the same topology --> TRUE
            return(TRUE)
          } else{
            print("Trees do not all have a unique topology.")
            return(FALSE)
          }
        }
      } else{
        return(TRUE) # All checks passed --> TRUE
      }
    }
}
