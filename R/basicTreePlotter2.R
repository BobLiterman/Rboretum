#' Rboretum Basic Tree Plotter
#' 
#' Given a phylo or multiPhylo object (tree), this ggtree wrapper returns a ggtree plot object, adjusted with possible arguments
#' @param tree Tree(s) to plot. Options include:
#' \itemize{
#'   \item A single, rooted phylo object [Clade support disabled]; or, 
#'   \item A named, rooted multiPhylo object [Only unique topologies plotted]
#' }
#' @param branch_length OPTIONAL: TRUE [plot tree(s) with branch lengths]; FALSE [Default: plot cladogram(s)]
#' @param branch_weight OPTIONAL: Set ggtree branch thickness
#' @param node_label OPTIONAL: Node label choices include:
#' \itemize{
#'   \item "none": No node labels
#'   \item "node": Node number
#'   \item "bs": Bootstrap value [Default]
#' }
#' @param node_size OPTIONAL: Set ggtree node label size
#' @param node_nudge OPTIONAL: Set ggtree node label nudge_x 
#' @param taxa_size OPTIONAL: Set ggtree tip label size
#' @param taxa_italic OPTIONAL: TRUE [italic tip labels]; FALSE [Default: non-italic tip labels]
#' @param taxa_align OPTIONAL: 'left' or 'right' tip label alignment [Default: No alignment]
#' @param taxa_offset OPTIONAL: Set ggtree tip label offset
#' @param xmax OPTIONAL: Set ggplot xlim upper limit (e.g if long tip labels run off plot)
#' @param reverse_x OPTIONAL: TRUE [plot tree with tips on left]; FALSE [Default: plot tree with tips on right]
#' @param to_color OPTIONAL: Color tips or clades via:
#' \itemize{
#'   \item Character vector of taxa, all to be colored the same color
#'   \item List of groups of taxa, each of which will have their own color. List can be named for use with a legend (set color_legend == TRUE)
#' }
#' @param colors OPTIONAL: Colors for clade highlighting. Must be hex or valid R colors. Provide a color for each group (1 if character vector, 1 for each group if named list) or default colors will be used.
#' @param color_legend OPTIONAL: TRUE [Include a legend for colored taxa/clades]; False [Default: No legend]
#' @return ggtree object
#' @export

basicTreePlotter2 <- function(tree,branch_length,branch_weight,node_label,node_size,node_nudge,taxa_size,taxa_italic,taxa_align,taxa_offset,xmax,reverse_x,to_color,colors,color_legend){
  
  if(!Rboretum::isMultiPhylo(tree,check_rooted = TRUE,check_names = TRUE) & !Rboretum::isPhylo(tree,check_rooted = TRUE)){
    stop("'tree' must be either a rooted phylo object or a named, rooted, mulitPhlyo object basicTreePlotter")
  }
  
  if(Rboretum::isPhylo(tree)){
    
    tree_taxa <- sort(tree$tip.label)
    tree_count <- 1
    
    titlePlot <- FALSE
    
  } else if(Rboretum::checkSameTopology(tree)){
    
    print("'tree' contains trees that all have the same topology. Processing first tree...")
    tree <- as.phylo(tree[[1]])
    tree_taxa <- sort(tree$tip.label)
    tree_count <- 1
    
    titlePlot <- FALSE

  } else if(!Rboretum::isMultiPhylo(tree,check_unique = TRUE)){
    
    print("'tree' contains trees that are not unique. Reducing to unique topologies...")

    print(Rboretum::getUniqueTopologies(tree,return_table = TRUE))
    tree <- Rboretum::getUniqueTopologies(tree)
    
    tree_taxa <- Rboretum::getSharedTaxa(tree)
    tree_count <- length(tree)
    tree_names <- names(tree)
    
    plot_title <- names(tree)
    titlePlot <- TRUE
    
    if(any(duplicated(tree_names))){
      stop("'tree' multiPhlyo contains trees with identical names.")
    }
  } else{
    
    tree_taxa <- Rboretum::getSharedTaxa(tree)
    tree_count <- length(tree)
    tree_names <- names(tree)
    
    plot_title <- names(tree)
    titlePlot <- TRUE
    
    if(any(duplicated(tree_names))){
      stop("'tree' multiPhlyo contains trees with identical names.")
    }
  }
  
  if(missing(branch_length)){
    branch_length <- FALSE
  } else if(!is.logical(branch_length)){
      branch_length <- FALSE
    }
  
  if(missing(branch_weight)){
    bWeight <- FALSE
  } else if(!is.numeric(branch_weight)){
      bWeight<-FALSE
    } else{ bWeight <- TRUE }
  
  if(missing(node_label)){
    node_label <- 'bs'
  } else if(!(node_label %in% c('bs','node','none'))){
    node_label <- 'bs'
  }
  
  if(missing(node_size)){
    nSize <- FALSE
  } else if(!is.numeric(node_size) | node_label == "none"){
    nSize<-FALSE
    } else{ nSize <- TRUE }
  
  if(missing(node_nudge)){
    nNudge <- FALSE
  } else if(!is.numeric(node_nudge)){
    nNudge<-FALSE
    } else{ nNudge <- TRUE }
  
  if(missing(taxa_size)){
    tSize <- FALSE
  } else if(!is.numeric(taxa_size)){
      tSize<-FALSE
    } else{ tSize <- TRUE }
  
  if(missing(taxa_italic)){
    taxa_italic <- FALSE
  } else if(!is.logical(taxa_italic)){
      taxa_italic <- FALSE
    }
  
  if(missing(taxa_align)){
    tAlign <- FALSE
  } else if(!(taxa_align %in%  c('left','right'))){
    tAlign <- FALSE
  } else{ tAlign <- TRUE }
  
  if(missing(taxa_offset)){
    tOffset <- FALSE
  } else if(!is.numeric(taxa_offset)){
      tOffset<-FALSE
    } else{ tOffset <- TRUE }
  
  if(missing(xmax)){
    extendX <- FALSE
  } else if(!is.numeric(xmax)){
      extendX<-FALSE
    } else{ extendX <- TRUE }
  
  if(missing(reverse_x)){
    reverseX <- FALSE
  } else if(!is.logical(reverse_x)){
    reverseX <- FALSE
  } else{ reverseX <- reverse_x }

  if(missing(to_color)){
    colorTips <- FALSE
  } else if(is.character(to_color)){
    colorTips <- TRUE
    group_count <- 1
  } else if(is.list(to_color)){
    colorTips <- TRUE
    group_count <- length(to_color)
  } else{
    print("'to_color' must be a character vector of taxa, or a list of taxa/clades")
    colorTips <- FALSE
  }
  
  if(colorTips){
    if(missing(colors)){
      if(group_count == 1){
        colors <- c('black','red')
      }
      else if(group_count <= 8){
        all_colors <- c('#800000','#4363d8','#f58231','#e6beff','#000075','#a9a9a9','#fabebe','#ffe119') 
        colors <- c('black',all_colors[1:group_count])
      } else{
        colors <- c('black',viridisLite::viridis(group_count))
      }
    } else{
      if(group_count == 1){
        
        if(length(colors)>1){
          print('Too many highlight colors provided. Using first color in list...')
          colors <- colors[1]
        }
        
        if(has_error(silent=TRUE,expr=grDevices::col2rgb(colors))){
          print('Invalid color choice, using defaults...')
          colors <- c('black','red')
        } else{colors <- c('black',colors)}
      
        } else if(group_count > 1){
          
          if(length(colors) != group_count | any(has_error(silent=TRUE,expr=grDevices::col2rgb(colors)))){
          print("Invalid color choice. Using defaults.")
          
          if(group_count <= 8){
            all_colors <- c('#800000','#4363d8','#f58231','#e6beff','#000075','#a9a9a9','#fabebe','#ffe119') 
            colors <- c('black',all_colors[1:group_count])
          } else{
            colors <- c('black',viridisLite::viridis(group_count))
          }
        } else{colors <- c('black',colors)}
      }
    }
  }
  
  if(missing(color_legend)){
    color_legend <- FALSE
  } else if(!is.logical(color_legend)){
    color_legend <- FALSE
  }
  
  plotList <- list()
  
  for(i in 1:tree_count){
  
    # Create base tree
    if(colorTips){
      tree <- ggtree::groupOTU(tree,to_color)
      
      if(is.character(to_color)){
        if(bWeight & branch_length){
          return_tree <- ggtree(tree,size=branch_weight,aes(color=group)) + scale_color_manual(values = colors)
        } else if(bWeight & !branch_length){
          return_tree <- ggtree(tree,size=branch_weight,branch.length = 'none',aes(color=group)) + scale_color_manual(values = colors)
        } else if(!bWeight & branch_length){
          return_tree <- ggtree(tree,aes(color=group)) + scale_color_manual(values = colors)
        } else if(!bWeight & !branch_length){
          return_tree <- ggtree(tree,branch.length = 'none',aes(color=group)) + scale_color_manual(values = colors)
        }
      } else{
        if(bWeight & branch_length){
          return_tree <- ggtree(tree,size=branch_weight,aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors)
        } else if(bWeight & !branch_length){
          return_tree <- ggtree(tree,size=branch_weight,branch.length = 'none',aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors)
        } else if(!bWeight & branch_length){
          return_tree <- ggtree(tree,aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors)
        } else if(!bWeight & !branch_length){
          return_tree <- ggtree(tree,branch.length = 'none',aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors)
        }
        
        if(color_legend){
          return_tree <- return_tree + labs(color = "Focal Clades") + theme(legend.position = 'right')
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
    
    if(titlePlot){
      if(colorTips & color_legend){
      return_tree <- return_tree + ggplot2::ggtitle(names(trees)[i]) + theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))
      } else{
        return_tree <- return_tree + ggplot2::ggtitle(names(trees)[i]) + theme(plot.title = element_text(hjust = 0.5))
      }
    }
    
    if(tree_count > 1){
      plotList[[i]] <- return_tree
    } else{
      return(return_tree)
    }
  }
  
  tandemPlotter(plotList)
  # return(plotList)
}