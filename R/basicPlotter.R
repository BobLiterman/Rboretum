#' Rboretum Basic Tree Plot Generator
#'
#' Given a phylo or multiPhylo object (tree), this ggtree wrapper returns a ggtree plot object, adjusted with possible arguments
#' @param tree Tree(s) to plot. Options include:
#' \itemize{
#'   \item A single, rooted phylo object; or,
#'   \item A named, rooted multiPhylo object [Only unique topologies plotted, from most to least common]
#' }
#' @param branch_length OPTIONAL: TRUE [plot tree(s) with branch lengths]; FALSE [Default: plot cladogram(s)]
#' @param branch_weight OPTIONAL: Set ggtree branch thickness
#' @param node_label OPTIONAL: Node label choices include:
#' \itemize{
#'   \item "none": No node labels
#'   \item "node": Node number
#'   \item "bs": Bootstrap value/node labels [Default]
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
#'   \item List of groups of taxa, each of which will have their own color. List can be named for use with a legend (set highlight_legend == TRUE)
#' }
#' @param colors OPTIONAL: Colors for clade highlighting. Must be hex or valid R colors. Provide a color for each group (1 if character vector, 1 for each group if named list) or default colors will be used.
#' @param highlight_legend OPTIONAL: TRUE [Include a legend for colored tips, given a list]; False [Default: No legend]
#' @param color_branches OPTIONAL: If TRUE and coloring taxa or clades, color the branches rather than the tip labels [Default: FALSE, colorize tip labels]
#' @param plot_title OPTIONAL: Character vector containing plot titles (1 per tree) [Default: No title for phylo, tree name for multiPhylo]
#' @return ggtree object or list of ggtree objects
#' @export

basicPlotter <- function(tree,branch_length,branch_weight,node_label,node_size,node_nudge,taxa_size,taxa_italic,taxa_align,taxa_offset,xmax,reverse_x,to_color,colors,highlight_legend,color_branches,plot_title){
  
  # Ensure tree is valid for plotter
  if(!Rboretum::isMultiPhylo(tree,check_rooted = TRUE,check_named = TRUE,check_three_taxa = TRUE) & !Rboretum::isPhylo(tree,check_rooted = TRUE)){
    stop("'tree' must be either a rooted phylo object or a named, rooted, mulitPhlyo object basicPlotter")
  }
  
  # If one tree is provided...
  if(Rboretum::isPhylo(tree)){
    
    tree_taxa <- sort(tree$tip.label)
    tree_count <- 1
    
    if(missing(plot_title)){
      titlePlot <- FALSE
    } else if(!is.character(plot_title)){
      titlePlot <- FALSE
    } else if(length(plot_title)!=1){
      titlePlot <- FALSE
    } else{
      titlePlot <- TRUE
    }
    
    # If a multiPhylo is provided, reduce to unique topologies...
  } else if(Rboretum::isMultiPhylo(tree,check_all_equal = TRUE)){
    
    tree_taxa <- Rboretum::getSharedTaxa(tree)
    tree <- Rboretum::treeTrimmer(tree[[1]],tree_taxa)
    tree_count <- 1
    
    if(missing(plot_title)){
      titlePlot <- FALSE
    } else if(!is.character(plot_title)){
      titlePlot <- FALSE
    } else if(length(plot_title)!=1){
      titlePlot <- FALSE
    } else{
      titlePlot <- TRUE
    }
    
  } else if(Rboretum::isMultiPhylo(tree,check_all_unique = TRUE)){
    
    tree_taxa <- Rboretum::getSharedTaxa(tree)
    tree <- Rboretum::treeTrimmer(tree,tree_taxa)
    tree_count <- length(tree)
    tree_names <- names(tree)
    
    if(any(duplicated(tree_names))){
      stop("'tree' multiPhlyo contains trees with identical names.")
    }
    
    titlePlot <- TRUE
    
    if(missing(plot_title)){
      plot_title <- tree_names
    } else if(!is.character(plot_title)){
      plot_title <- tree_names
    } else if(length(plot_title)!=tree_count){
      plot_title <- tree_names
    }
    
  } else if(Rboretum::isMultiPhylo(tree,check_some_equal = TRUE)){
    
    tree_taxa <- Rboretum::getSharedTaxa(tree)
    tree_table <- Rboretum::getUniqueTopologies(tree,return_table = TRUE)
    tree <- Rboretum::getUniqueTopologies(tree)
    tree_count <- length(tree)
    tree_names <- as.character(tree_table$Trees_with_Topology)
    
    if(any(duplicated(tree_names))){
      stop("'tree' multiPhlyo contains trees with identical names.")
    }
    
    titlePlot <- TRUE
    
    if(missing(plot_title)){
      plot_title <- tree_names
    } else if(!is.character(plot_title)){
      plot_title <- tree_names
    } else if(length(plot_title)!=tree_count){
      plot_title <- tree_names
    }
  }
  
  # Get other arguments and establish defaults
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
          colors <- c('black',viridisLite::viridis(group_count))
        } else{colors <- c('black',colors)
        }
      }
    }
  }
  
  if(missing(highlight_legend)){
    highlight_legend <- FALSE
  } else if(!is.logical(highlight_legend)){
    highlight_legend <- FALSE
  }
  
  if(missing(color_branches)){
    color_branches <- FALSE
  } else if(!is.logical(color_branches)){
    color_branches <- FALSE
  } else if(!colorTips){
    color_branches <- FALSE
  }
  
  # Create empty plot list
  plotList <- list()
  
  # Create dummy multiPhylo for simple handling
  if(tree_count == 1){
    tree <- c(tree,tree)
  }
  
  for(i in 1:tree_count){
    temp_tree <- tree[[i]]
    
    temp_highlight_legend <- highlight_legend
    
    if(i!=tree_count){
      temp_highlight_legend <- FALSE
    }
    
    # Create base tree
    
    if(colorTips){
      temp_tree <- ggtree::groupOTU(temp_tree,to_color)
    }
    
    if(colorTips & color_branches){
      
      if(is.character(to_color)){
        if(bWeight & branch_length){
          return_tree <- ggtree(temp_tree,size=branch_weight,aes(color=group)) + scale_color_manual(values = colors)
        } else if(bWeight & !branch_length){
          return_tree <- ggtree(temp_tree,size=branch_weight,branch.length = 'none',aes(color=group)) + scale_color_manual(values = colors)
        } else if(!bWeight & branch_length){
          return_tree <- ggtree(temp_tree,aes(color=group)) + scale_color_manual(values = colors)
        } else if(!bWeight & !branch_length){
          return_tree <- ggtree(temp_tree,branch.length = 'none',aes(color=group)) + scale_color_manual(values = colors)
        } #+ ggplot2::ggtitle(plot_title[i]) + labs(color = "Focal Clades") + theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))
      } else{
        if(titlePlot){
          if(temp_highlight_legend){
            if(bWeight & branch_length){
              return_tree <- ggtree(temp_tree,size=branch_weight,aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                ggplot2::ggtitle(plot_title[i]) + 
                labs(color = "Focal Clades") + 
                theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))
            } else if(bWeight & !branch_length){
              return_tree <- ggtree(temp_tree,size=branch_weight,branch.length = 'none',aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                ggplot2::ggtitle(plot_title[i]) + 
                labs(color = "Focal Clades") + 
                theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))
            } else if(!bWeight & branch_length){
              return_tree <- ggtree(temp_tree,aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                ggplot2::ggtitle(plot_title[i]) + 
                labs(color = "Focal Clades") + 
                theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))
            } else if(!bWeight & !branch_length){
              return_tree <- ggtree(temp_tree,branch.length = 'none',aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                ggplot2::ggtitle(plot_title[i]) + 
                labs(color = "Focal Clades") + 
                theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))
            }
          } else{
            if(bWeight & branch_length){
              return_tree <- ggtree(temp_tree,size=branch_weight,aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                ggplot2::ggtitle(plot_title[i]) + theme(plot.title = element_text(hjust = 0.5))

                
            } else if(bWeight & !branch_length){
              return_tree <- ggtree(temp_tree,size=branch_weight,branch.length = 'none',aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                ggplot2::ggtitle(plot_title[i]) + theme(plot.title = element_text(hjust = 0.5))
            } else if(!bWeight & branch_length){
              return_tree <- ggtree(temp_tree,aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                ggplot2::ggtitle(plot_title[i]) + theme(plot.title = element_text(hjust = 0.5))
            } else if(!bWeight & !branch_length){
              return_tree <- ggtree(temp_tree,branch.length = 'none',aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                ggplot2::ggtitle(plot_title[i]) + theme(plot.title = element_text(hjust = 0.5))
            }
          }
        } else{
          if(temp_highlight_legend){
            if(bWeight & branch_length){
              return_tree <- ggtree(temp_tree,size=branch_weight,aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                labs(color = "Focal Clades") + 
                theme(legend.position = 'right')
            } else if(bWeight & !branch_length){
              return_tree <- ggtree(temp_tree,size=branch_weight,branch.length = 'none',aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                labs(color = "Focal Clades") + 
                theme(legend.position = 'right')
            } else if(!bWeight & branch_length){
              return_tree <- ggtree(temp_tree,aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                labs(color = "Focal Clades") + 
                theme(legend.position = 'right')
            } else if(!bWeight & !branch_length){
              return_tree <- ggtree(temp_tree,branch.length = 'none',aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                labs(color = "Focal Clades") + 
                theme(legend.position = 'right')
            }
          } else{
            if(bWeight & branch_length){
              return_tree <- ggtree(temp_tree,size=branch_weight,aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                labs(color = "Focal Clades") + 
                theme(legend.position = 'right')
            } else if(bWeight & !branch_length){
              return_tree <- ggtree(temp_tree,size=branch_weight,branch.length = 'none',aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                labs(color = "Focal Clades") + 
                theme(legend.position = 'right')
            } else if(!bWeight & branch_length){
              return_tree <- ggtree(temp_tree,aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                labs(color = "Focal Clades") + 
                theme(legend.position = 'right')
            } else if(!bWeight & !branch_length){
              return_tree <- ggtree(temp_tree,branch.length = 'none',aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors) + 
                labs(color = "Focal Clades") + 
                theme(legend.position = 'right')
            }
          }
        }
      }
    } else{
      if(bWeight & branch_length){
        return_tree <- ggtree(temp_tree,size=branch_weight)
      } else if(bWeight & !branch_length){
        return_tree <- ggtree(temp_tree,size=branch_weight,branch.length = 'none')
      } else if(!bWeight & branch_length){
        return_tree <- ggtree(temp_tree)
      } else if(!bWeight & !branch_length){
        return_tree <- ggtree(temp_tree,branch.length = 'none')
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
    if(colorTips & !color_branches){
      
      if(is.character(to_color)){
        if(!tAlign & !tOffset){
          if(tSize & taxa_italic){
            return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic') + scale_color_manual(values = colors)
          } else if(tSize & !taxa_italic){
            return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size) + scale_color_manual(values = colors)
          } else if(!tSize & taxa_italic){
            return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic') + scale_color_manual(values = colors)
          } else{
            return_tree <- return_tree + geom_tiplab(aes(color=group)) + scale_color_manual(values = colors)
          }
        } else if(!tAlign & tOffset){
          if(tSize & taxa_italic){
            return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic',offset=taxa_offset) + scale_color_manual(values = colors)
          } else if(tSize & !taxa_italic){
            return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,offset=taxa_offset) + scale_color_manual(values = colors)
          } else if(!tSize & taxa_italic){
            return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic',offset=taxa_offset) + scale_color_manual(values = colors)
          } else{
            return_tree <- return_tree + geom_tiplab(aes(color=group),offset=taxa_offset) + scale_color_manual(values = colors)
          }
        } else if(tAlign & !tOffset){
          if(taxa_align == 'right'){
            if(branch_legnth){
              if(tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic',hjust=1,align=TRUE) + scale_color_manual(values = colors)
              } else if(tSize & !taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,hjust=1,align=TRUE) + scale_color_manual(values = colors)
              } else if(!tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic',hjust=1,align=TRUE) + scale_color_manual(values = colors)
              } else{
                return_tree <- return_tree + geom_tiplab(aes(color=group),hjust=1,align=TRUE) + scale_color_manual(values = colors)
              }
            } else{
              if(tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic',hjust=1) + scale_color_manual(values = colors)
              } else if(tSize & !taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,hjust=1) + scale_color_manual(values = colors)
              } else if(!tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic',hjust=1) + scale_color_manual(values = colors)
              } else{
                return_tree <- return_tree + geom_tiplab(aes(color=group),hjust=1) + scale_color_manual(values = colors)
              }
            }
          } else{
            if(tSize & taxa_italic){
              return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic',hjust=0,align = TRUE) + scale_color_manual(values = colors)
            } else if(tSize & !taxa_italic){
              return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,hjust=0,align = TRUE) + scale_color_manual(values = colors)
            } else if(!tSize & taxa_italic){
              return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic',hjust=0,align = TRUE) + scale_color_manual(values = colors)
            } else{
              return_tree <- return_tree + geom_tiplab(aes(color=group),hjust=0,align = TRUE) + scale_color_manual(values = colors)
            }
          }
        } else if(tAlign & tOffset){
          if(taxa_align == 'right'){
            if(branch_length){
              if(tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic',hjust=1,offset=taxa_offset,align=TRUE) + scale_color_manual(values = colors)
              } else if(tSize & !taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,hjust=1,offset=taxa_offset,align=TRUE) + scale_color_manual(values = colors)
              } else if(!tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic',hjust=1,offset=taxa_offset,align=TRUE) + scale_color_manual(values = colors)
              } else{
                return_tree <- return_tree + geom_tiplab(aes(color=group),hjust=1,offset=taxa_offset,align=TRUE) + scale_color_manual(values = colors)
              }
            } else{
              if(tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic',hjust=1,offset=taxa_offset) + scale_color_manual(values = colors)
              } else if(tSize & !taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,hjust=1,offset=taxa_offset) + scale_color_manual(values = colors)
              } else if(!tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic',hjust=1,offset=taxa_offset) + scale_color_manual(values = colors)
              } else{
                return_tree <- return_tree + geom_tiplab(aes(color=group),hjust=1,offset=taxa_offset) + scale_color_manual(values = colors)
              }
            }
          } else{
            if(tSize & taxa_italic){
              return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic',hjust=0,offset=taxa_offset,align = TRUE) + scale_color_manual(values = colors)
            } else if(tSize & !taxa_italic){
              return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,hjust=0,offset=taxa_offset,align = TRUE) + scale_color_manual(values = colors)
            } else if(!tSize & taxa_italic){
              return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic',hjust=0,offset=taxa_offset,align = TRUE) + scale_color_manual(values = colors)
            } else{
              return_tree <- return_tree + geom_tiplab(aes(color=group),hjust=0,offset=taxa_offset,align = TRUE) + scale_color_manual(values = colors)
            }
          }
        }
        
      } else{
        
        if(!tAlign & !tOffset){
          if(tSize & taxa_italic){
            return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic') + scale_color_manual(breaks = names(to_color),values = colors)
          } else if(tSize & !taxa_italic){
            return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size) + scale_color_manual(breaks = names(to_color),values = colors)
          } else if(!tSize & taxa_italic){
            return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic') + scale_color_manual(breaks = names(to_color),values = colors)
          } else{
            return_tree <- return_tree + geom_tiplab(aes(color=group)) + scale_color_manual(breaks = names(to_color),values = colors)
          }
        } else if(!tAlign & tOffset){
          if(tSize & taxa_italic){
            return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic',offset=taxa_offset) + scale_color_manual(breaks = names(to_color),values = colors)
          } else if(tSize & !taxa_italic){
            return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,offset=taxa_offset) + scale_color_manual(breaks = names(to_color),values = colors)
          } else if(!tSize & taxa_italic){
            return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic',offset=taxa_offset) + scale_color_manual(breaks = names(to_color),values = colors)
          } else{
            return_tree <- return_tree + geom_tiplab(aes(color=group),offset=taxa_offset) + scale_color_manual(breaks = names(to_color),values = colors)
          }
        } else if(tAlign & !tOffset){
          if(taxa_align == 'right'){
            if(branch_legnth){
              if(tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic',hjust=1,align=TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
              } else if(tSize & !taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,hjust=1,align=TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
              } else if(!tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic',hjust=1,align=TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
              } else{
                return_tree <- return_tree + geom_tiplab(aes(color=group),hjust=1,align=TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
              }
            } else{
              if(tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic',hjust=1) + scale_color_manual(breaks = names(to_color),values = colors)
              } else if(tSize & !taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,hjust=1) + scale_color_manual(breaks = names(to_color),values = colors)
              } else if(!tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic',hjust=1) + scale_color_manual(breaks = names(to_color),values = colors)
              } else{
                return_tree <- return_tree + geom_tiplab(aes(color=group),hjust=1) + scale_color_manual(breaks = names(to_color),values = colors)
              }
            }
          } else{
            if(tSize & taxa_italic){
              return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic',hjust=0,align = TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
            } else if(tSize & !taxa_italic){
              return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,hjust=0,align = TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
            } else if(!tSize & taxa_italic){
              return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic',hjust=0,align = TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
            } else{
              return_tree <- return_tree + geom_tiplab(aes(color=group),hjust=0,align = TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
            }
          }
        } else if(tAlign & tOffset){
          if(taxa_align == 'right'){
            if(branch_length){
              if(tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic',hjust=1,offset=taxa_offset,align=TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
              } else if(tSize & !taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,hjust=1,offset=taxa_offset,align=TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
              } else if(!tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic',hjust=1,offset=taxa_offset,align=TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
              } else{
                return_tree <- return_tree + geom_tiplab(aes(color=group),hjust=1,offset=taxa_offset,align=TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
              }
            } else{
              if(tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic',hjust=1,offset=taxa_offset) + scale_color_manual(breaks = names(to_color),values = colors)
              } else if(tSize & !taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,hjust=1,offset=taxa_offset) + scale_color_manual(breaks = names(to_color),values = colors)
              } else if(!tSize & taxa_italic){
                return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic',hjust=1,offset=taxa_offset) + scale_color_manual(breaks = names(to_color),values = colors)
              } else{
                return_tree <- return_tree + geom_tiplab(aes(color=group),hjust=1,offset=taxa_offset) + scale_color_manual(breaks = names(to_color),values = colors)
              }
            }
          } else{
            if(tSize & taxa_italic){
              return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,fontface='italic',hjust=0,offset=taxa_offset,align = TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
            } else if(tSize & !taxa_italic){
              return_tree <- return_tree + geom_tiplab(aes(color=group),size=taxa_size,hjust=0,offset=taxa_offset,align = TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
            } else if(!tSize & taxa_italic){
              return_tree <- return_tree + geom_tiplab(aes(color=group),fontface='italic',hjust=0,offset=taxa_offset,align = TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
            } else{
              return_tree <- return_tree + geom_tiplab(aes(color=group),hjust=0,offset=taxa_offset,align = TRUE) + scale_color_manual(breaks = names(to_color),values = colors)
            }
          }
        }
        
      }
      
    } else {
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
    }
    
    if(titlePlot){
      if(temp_highlight_legend){
        return_tree <- return_tree + ggplot2::ggtitle(plot_title[i]) + labs(color = "Focal Clades") + theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))
      } else{
        return_tree <- return_tree + ggplot2::ggtitle(plot_title[i]) + theme(plot.title = element_text(hjust = 0.5))
      }
    }
    
    if(tree_count > 1){
      plotList[[i]] <- return_tree
    } else{
      return(return_tree)
    }
  }
  return(plotList)
}