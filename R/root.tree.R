#' Rboretum Tree Rooter
#'
#' This ape wrapper accepts a phylo object and a character vector of new root taxa and returns a rooted phylo object, if possible
#' @param tree Phylo object
#' @param root_taxa Character vector containing outgroup species IDs (Must be in tree and monophyletic)
#' @param root_node Root by node value rather than root taxa. NOTE: node setting overrides root_taxa
#' @param resolve_root OPTIONAL: Set ape::resolve.root [Default: TRUE]
#' @return A phylo object rooted at specified taxa
#' @export root.tree
#' @examples
#' root.tree(birdTree,c('Alligator','Turtle'))
#' 
root.tree <- function(tree,root_taxa,root_node,resolve_root){
  
  if(missing(resolve_root)){
    resolve_root <- TRUE
  } else if(!is.logical(resolve_root)){
    resolve_root <- TRUE
  }
  
  if(missing(root_taxa)){
    
    if(missing(root_node)){ stop("Must set either root_taxa or root_node.") } 
    
    else if(has_error(ape::root.phylo(tree,node = node,edgelabel = TRUE,resolve.root = resolve_root))){ stop("Ape cannot root tree at the specified node.") } 
    
    else{ return(ape::root.phylo(tree,node = node,edgelabel = TRUE,resolve.root = resolve_root)) }
  
    } else if(missing(root_node)){
      
      if(has_error(ape::root.phylo(tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = resolve_root))){ stop("Ape cannot root tree on these taxa.") }
      
      else{ return(ape::root.phylo(tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = resolve_root)) }
    
    } else if(has_error(ape::root.phylo(tree,node = node,edgelabel = TRUE,resolve.root = resolve_root))){ stop("Ape cannot root tree at the specified node.") } 
  
  else{ return(ape::root.phylo(tree,node = node,edgelabel = TRUE,resolve.root = resolve_root))}
  
}