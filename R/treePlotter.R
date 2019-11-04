#' Rboretum Tree Plotter
#' 
#' Given a phylo object (tree), this ggtree wrapper plots the tree with following optional arguments:
#' - node_label: 'none' [no node labels], 'node' [print node number labels] or 'bs' [print bootstrap labels] (Default: Plot with bootstrap labels, if present)
#' - node_size: Set ggtree node label size (Default: ggtree default)
#' - node_nudge: Set ggtree nudge_x for node label
#' - branch_length: TRUE or FALSE; Plot tree with or without branch lengths (Default: FALSE; Plot as cladogram) 
#' - branch_weight: Set  ggtree branch thickness (Default: ggTree default)
#' - taxa_size: Set ggtree tip label size (Default: ggTree default)
#' - taxa_italic: TRUE of FALSE; Print tip labels in italic (Default: False)
#' - taxa_align: 'left' or 'right' alignment of tip labels
#' - taxa_offset: Set ggtree offset tip label value (Default: No offset)
#' - xmax: numeric extension value for ggplot2::xlim(0,xmax) (e.g. if tip labels are cut-off)
#'
#' @param tree Phylo object
#' @param node_label OPTIONAL: 'none', 'node', or 'bs' [Default: bs (print bootstrap)]
#' @param node_size OPTIONAL: Override default ggtree node label size
#' @param node_nudge OPTIONAL: Override default ggtree node label placement 
#' @param branch_length OPTIONAL: TRUE or FALSE [Default: FALSE (print cladogram)]
#' @param branch_weight OPTIONAL: Override default ggtree branch thickness
#' @param taxa_size OPTIONAL: Override default ggtree tip label size
#' @param taxa_italic OPTIONAL: TRUE or FALSE [Default: FALSE]
#' @param taxa_align OPTIONAL: 'left' or 'right' tip label alignment [Default: No alignment]
#' @param taxa_offset OPTIONAL: Override default ggtree tip label placement
#' @param xmax OPTIONAL: Override default ggplot xlim upper limit (e.g + xlim(0,max))
#' @export
#' @examples
#' # Print tree with bootstrap labels
#' treePlotter(tree)
#' 
#' # Print tree with no node labels
#' treePlotter(tree,node_label='none')
#'
#' # Print tree with thicker branches
#' treePlotter(tree,branch_weight=3)
#' 
#' # Print tree with buffer on right side to accomodate longer tip labels
#' treePlotter(tree,xmax=10)

treePlotter <- function(tree,branch_length,branch_weight,node_label,node_size,node_nudge,taxa_size,taxa_italic,taxa_align,taxa_offset,xmax){
  
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
  
  # Create base tree
  
  if(bWeight & branch_length){
    return_tree <- ggtree(tree,size=branch_weight)
  } else if(bWeight & !branch_length){
    return_tree <- ggtree(tree,size=branch_weight,branch.length = 'none')
  } else if(!bWeight & branch_length){
    return_tree <- ggtree(tree)
  } else if(!bWeight & !branch_length){
    return_tree <- ggtree(tree,branch.length = 'none')
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
  
  # Return with or without extended x-axis
  
  if(extendX){
    return(return_tree + ggplot2::xlim(0,xmax))
  } else{
    return(return_tree)
  }
}