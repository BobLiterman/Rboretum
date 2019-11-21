#' Rboretum Basic Tree Plotter
#' 
#' Given a phylo object (tree), this ggtree wrapper returns a ggtree plot object, adjusted with possible arguments
#' @usage myPlot <- treePlotter(tree=myTree,...)
#' @param tree Phylo object
#' @param node_label OPTIONAL: Node label choices include:
#' \itemize{
#'   \item "none": No node labels
#'   \item "node": Node number
#'   \item "bs": Bootstrap value [Default]
#' }
#' @param node_size OPTIONAL: Set ggtree node label size
#' @param node_nudge OPTIONAL: Set ggtree node label nudge_x 
#' @param branch_length OPTIONAL: TRUE [plot tree with branch lengths]; FALSE [Default: plot cladogram]
#' @param branch_weight OPTIONAL: Set ggtree branch thickness
#' @param taxa_size OPTIONAL: Set ggtree tip label size
#' @param taxa_italic OPTIONAL: TRUE [italic tip labels]; FALSE [Default: non-italic tip labels]
#' @param taxa_align OPTIONAL: 'left' or 'right' tip label alignment [Default: No alignment]
#' @param taxa_offset OPTIONAL: Set ggtree tip label offset
#' @param xmax OPTIONAL: Set ggplot xlim upper limit (e.g if long tip labels run off plot)
#' @param reverse_x OPTIONAL: TRUE [plot tree with tips on left]; FALSE [Default: plot tree with tips on right]
#' @param rename_tips OPTIONAL: Dataframe where column 1 contains tip labels for tree, and column 2 contains new, desired tip labels
#' @param to_color OPTIONAL: Color tips or clades via:
#' \itemize{
#'   \item Character vector of taxa, all to be colored the same color
#'   \item List of groups of taxa, each of which will have their own color. List can be named for use with a legend (set color_legend == TRUE)
#' }
#' @param color_legend TRUE [Include a legend for colored taxa/clades]; False [Default: No legend]
#' @return ggtree object
#' @export
#' @examples
#' # Print tree with bootstrap labels
#' basic.treePlot(tree)
#' 
#' # Print tree with no node labels
#' basic.treePlot(tree,node_label='none')
#'
#' # Print tree with thicker branches
#' basic.treePlot(tree,branch_weight=3)
#' 
#' # Print tree with buffer on right side to accomodate longer tip labels
#' basic.treePlot(tree,xmax=10)

basic.treePlot2 <- function(tree,branch_length,branch_weight,node_label,node_size,node_nudge,taxa_size,taxa_italic,taxa_align,taxa_offset,xmax,reverse_x,rename_tips,to_color,color_legend){
  
  if(missing(tree)){
    stop("No tree provided.")
  }
  
  if(missing(branch_length)){
    branch_length <- FALSE
  } else{
    if(branch_length != TRUE){
      branch_length <- FALSE
    }
  }
  
  if(missing(branch_weight)){
    bWeight <- FALSE
  } else{
    if(!is.numeric(branch_weight)){
      bWeight<-FALSE
    } else{ bWeight <- TRUE }
    }
  
  if(missing(node_label)){
    node_label <- 'bs'
  } else if(!any(node_label ==  c('bs','node','none'))){
    node_label <- 'bs'
  }
  
  if(missing(node_size)){
    nSize <- FALSE
  } else{
    if(!is.numeric(node_size) | node_label == "none"){
      nSize<-FALSE
    } else{ nSize <- TRUE }
  }
  
  if(missing(node_nudge)){
    nNudge <- FALSE
  } else{
    if(!is.numeric(node_nudge)){
      nNudge<-FALSE
    } else{ nNudge <- TRUE }
  }
  
  if(missing(taxa_size)){
    tSize <- FALSE
  } else{
    if(!is.numeric(taxa_size)){
      tSize<-FALSE
    } else{ tSize <- TRUE }
    }
  
  if(missing(taxa_italic)){
    taxa_italic <- FALSE
  } else{
    if(taxa_italic != TRUE){
      taxa_italic <- FALSE
    }
  }
  
  if(missing(taxa_align)){
    tAlign <- FALSE
  } else if(!any(taxa_align ==  c('left','right'))){
    tAlign <- FALSE
  } else{ tAlign <- TRUE }
  
  if(missing(taxa_offset)){
    tOffset <- FALSE
  } else{
    if(!is.numeric(taxa_offset)){
      tOffset<-FALSE
    } else{ tOffset <- TRUE }
  }
  
  if(missing(xmax)){
    extendX <- FALSE
  } else{
    if(!is.numeric(xmax)){
      extendX<-FALSE
    } else{ extendX <- TRUE }
  }
  
  if(missing(reverse_x)){
    reverseX <- FALSE
  } else if((reverse_x != TRUE)){
    reverseX <- FALSE
  } else{
    reverseX <- TRUE
  }
  
  if(missing(rename_tips)){
    tree <- tree
  } else if((is.data.frame(rename_tips) | is_tibble(rename_tips)) & length(names(rename_tips) >= 2)){
    old_id <- names(rename_tips)[1]
    new_id <- names(rename_tips)[2]
    if(has_error(Rboretum::convert.tips(tree,rename_tips,old_id,new_id))){
      tree <- tree
    } else{
      tree <- Rboretum::convert.tips(tree,rename_tips,old_id,new_id)
      name_dict <- list(new_id)
      names(name_dict) <- old_id
    }
  } else{
    tree <- tree
  }
  
  if(missing(to_color)){
    colorTips <- FALSE
  } else if(is.character(to_color)){
    if(check.tip(tree,to_color)){
      colorTips <- TRUE
    } else{
      print("Some taxa in 'to_color' not found in tree.")
      colorTips <- FALSE
    }
  } else if(is.list(to_color)){
    if(check.tip(tree,unlist(to_color))){
      colorTips <- TRUE
    } else{
      colorTips <- FALSE
      print("Some taxa in 'to_color' not found in tree.")
    }
  } else{
    print("'to_color' must be a character vector of taxa, or a list of taxa/clades")
    colorTips <- FALSE
  }
  
  if(missing(color_legend)){
    legendColor <- FALSE
  } else if(!is.logical(color_legend)){
    legendColor <- FALSE
  }
  
  # Create base tree
  if(colorTips){
    tree <- ggtree::groupOTU(tree,to_color)
    
    if(is.character(to_color)){
      if(bWeight & branch_length){
        return_tree <- ggtree(tree,size=branch_weight,aes(color=group))
      } else if(bWeight & !branch_length){
        return_tree <- ggtree(tree,size=branch_weight,branch.length = 'none',aes(color=group))
      } else if(!bWeight & branch_length){
        return_tree <- ggtree(tree,aes(color=group))
      } else if(!bWeight & !branch_length){
        return_tree <- ggtree(tree,branch.length = 'none',aes(color=group))
      }
    } else{
      if(bWeight & branch_length){
        return_tree <- ggtree(tree,size=branch_weight)
      } else if(bWeight & !branch_length){
        return_tree <- ggtree(tree,size=branch_weight,branch.length = 'none')
      } else if(!bWeight & branch_length){
        return_tree <- ggtree(tree)
      } else if(!bWeight & !branch_length){
        return_tree <- ggtree(tree,branch.length = 'none')
      }
    }
  } else{
    if(bWeight & branch_length){
      return_tree <- ggtree(tree,size=branch_weight)
    } else if(bWeight & !branch_length){
      return_tree <- ggtree(tree,size=branch_weight,branch.length = 'none')
    } else if(!bWeight & branch_length){
      return_tree <- ggtree(tree)
    } else if(!bWeight & !branch_length){
      return_tree <- ggtree(tree,branch.length = 'none')
    }
  }
  if(extendX){
    if(!reverseX){
    return_tree <- return_tree + ggplot2::xlim(0,xmax)
    }
  }
  
  if(reverseX){
    return_tree <- return_tree + scale_x_reverse()
  }
  
  # Process node labels
  
  if(node_label == "none"){
    return_tree <- return_tree
  } else{
    if(nNudge){
      if(nSize){
        if(node_label == "node"){
          return_tree <- return_tree + geom_nodelab(nudge_x = node_nudge,size=node_size,aes(label=node))
        } else{
          return_tree <- return_tree + geom_nodelab(nudge_x = node_nudge,size=node_size)
        }
      } else{
        if(node_label == "node"){
          return_tree <- return_tree + geom_nodelab(nudge_x = node_nudge,aes(label=node))
          } else{
            return_tree <- return_tree + geom_nodelab(nudge_x = node_nudge)
          }
      }
    } else{
      if(nSize){
        if(node_label == "node"){
          return_tree <- return_tree + geom_nodelab(size=node_size,aes(label=node))
        } else{
          return_tree <- return_tree + geom_nodelab(size=node_size)
        }
      } else{
        if(node_label == "node"){
          return_tree <- return_tree + geom_nodelab(aes(label=node))
        } else{
          return_tree <- return_tree + geom_nodelab()
        }
      }
    }
  }
  
  # Process tip labels
  if(!tAlign & !tOffset){
    if(tSize & taxa_italic){
      return_tree <- return_tree + geom_tiplab(size=taxa_size,fontface='italic')
    } else if(tSize & !taxa_italic){
      return_tree <- return_tree + geom_tiplab(size=taxa_size)
    } else if(!tSize & taxa_italic){
      return_tree <- return_tree + geom_tiplab(fontface='italic')
    } else{
      return_tree <- return_tree + geom_tiplab()
    }
  } else if(!tAlign & tOffset){
    if(tSize & taxa_italic){
      return_tree <- return_tree + geom_tiplab(size=taxa_size,fontface='italic',offset=taxa_offset)
    } else if(tSize & !taxa_italic){
      return_tree <- return_tree + geom_tiplab(size=taxa_size,offset=taxa_offset)
    } else if(!tSize & taxa_italic){
      return_tree <- return_tree + geom_tiplab(fontface='italic',offset=taxa_offset)
    } else{
      return_tree <- return_tree + geom_tiplab(offset=taxa_offset)
    }
  } else if(tAlign & !tOffset){
    if(taxa_align == 'right'){
      if(branch_legnth){
        if(tSize & taxa_italic){
          return_tree <- return_tree + geom_tiplab(size=taxa_size,fontface='italic',hjust=1,align=TRUE)
        } else if(tSize & !taxa_italic){
          return_tree <- return_tree + geom_tiplab(size=taxa_size,hjust=1,align=TRUE)
        } else if(!tSize & taxa_italic){
          return_tree <- return_tree + geom_tiplab(fontface='italic',hjust=1,align=TRUE)
        } else{
          return_tree <- return_tree + geom_tiplab(hjust=1,align=TRUE)
        }
      } else{
        if(tSize & taxa_italic){
          return_tree <- return_tree + geom_tiplab(size=taxa_size,fontface='italic',hjust=1)
        } else if(tSize & !taxa_italic){
          return_tree <- return_tree + geom_tiplab(size=taxa_size,hjust=1)
        } else if(!tSize & taxa_italic){
          return_tree <- return_tree + geom_tiplab(fontface='italic',hjust=1)
        } else{
          return_tree <- return_tree + geom_tiplab(hjust=1)
        }
      }
    } else{
      if(tSize & taxa_italic){
        return_tree <- return_tree + geom_tiplab(size=taxa_size,fontface='italic',hjust=0,align = TRUE)
      } else if(tSize & !taxa_italic){
        return_tree <- return_tree + geom_tiplab(size=taxa_size,hjust=0,align = TRUE)
      } else if(!tSize & taxa_italic){
        return_tree <- return_tree + geom_tiplab(fontface='italic',hjust=0,align = TRUE)
      } else{
        return_tree <- return_tree + geom_tiplab(hjust=0,align = TRUE)
      }
    }
  } else if(tAlign & tOffset){
    if(taxa_align == 'right'){
      if(branch_length){
        if(tSize & taxa_italic){
          return_tree <- return_tree + geom_tiplab(size=taxa_size,fontface='italic',hjust=1,offset=taxa_offset,align=TRUE)
        } else if(tSize & !taxa_italic){
          return_tree <- return_tree + geom_tiplab(size=taxa_size,hjust=1,offset=taxa_offset,align=TRUE)
        } else if(!tSize & taxa_italic){
          return_tree <- return_tree + geom_tiplab(fontface='italic',hjust=1,offset=taxa_offset,align=TRUE)
        } else{
          return_tree <- return_tree + geom_tiplab(hjust=1,offset=taxa_offset,align=TRUE)
        }
      } else{
        if(tSize & taxa_italic){
          return_tree <- return_tree + geom_tiplab(size=taxa_size,fontface='italic',hjust=1,offset=taxa_offset)
        } else if(tSize & !taxa_italic){
          return_tree <- return_tree + geom_tiplab(size=taxa_size,hjust=1,offset=taxa_offset)
        } else if(!tSize & taxa_italic){
          return_tree <- return_tree + geom_tiplab(fontface='italic',hjust=1,offset=taxa_offset)
        } else{
          return_tree <- return_tree + geom_tiplab(hjust=1,offset=taxa_offset)
        }
        }
      } else{
      if(tSize & taxa_italic){
        return_tree <- return_tree + geom_tiplab(size=taxa_size,fontface='italic',hjust=0,offset=taxa_offset,align = TRUE)
      } else if(tSize & !taxa_italic){
        return_tree <- return_tree + geom_tiplab(size=taxa_size,hjust=0,offset=taxa_offset,align = TRUE)
      } else if(!tSize & taxa_italic){
        return_tree <- return_tree + geom_tiplab(fontface='italic',hjust=0,offset=taxa_offset,align = TRUE)
      } else{
        return_tree <- return_tree + geom_tiplab(hjust=0,offset=taxa_offset,align = TRUE)
      }
    }
  }

  return(return_tree)
}