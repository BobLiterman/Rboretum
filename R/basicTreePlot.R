#' Basic Tree Plotter
#'
#' This function plots a tree without branch lengths, and nodes can be labeled with bootstrap or node labels
#' @param tree_1 Phylo object
#' @param label Either 'node' [print node labels] or 'bs' [print bootstrap labels]. Any other entry (or missing) results in no node labels
#' @param xmax Optional: Provide a numeric value to extend plot axis (i.e. if right side is cut off)
#' @export
#' @examples
#' # Print tree with no labels
#' basicTreePlot(tree)
#'
#' # Print tree with bootstrap labels
#' basicTreePlot(tree,'bs')
#'
#' # Print tree with node labels
#' basicTreePlot(tree,'node')
#'

basicTreePlot <- function(tree,label,xmax){

  if(missing(xmax)){
    if(missing(label)){
      return(ggtree(tree,branch.length = 'none')+geom_tiplab())
    }
    else{
      if(label=="node"){
        return(ggtree(tree,branch.length = 'none')+geom_tiplab()+geom_nodelab(aes(label=node)))
      }
      else if(label=="bs"){
        return(ggtree(tree,branch.length = 'none')+geom_tiplab()+geom_nodelab())
      }
      else if(label != "node" & label != 'bs'){
        return(ggtree(tree,branch.length = 'none')+geom_tiplab())
      }
    }
  }
  else{
    if(is.numeric(xmax)){
      if(missing(label)){
        return(ggtree(tree,branch.length = 'none')+geom_tiplab() + xlim(0,xmax))
      }
      else{
        if(label=="node"){
          return(ggtree(tree,branch.length = 'none')+geom_tiplab()+geom_nodelab(aes(label=node))+ xlim(0,xmax))
        }
        else if(label=="bs"){
          return(ggtree(tree,branch.length = 'none')+geom_tiplab()+geom_nodelab()+ xlim(0,xmax))
        }
        else if(label != "node" & label != 'bs'){
          return(ggtree(tree,branch.length = 'none')+geom_tiplab()+ xlim(0,xmax))
        }
      }
    }
    else{
      stop("xmax argument must be numeric")
    }
  }

}
