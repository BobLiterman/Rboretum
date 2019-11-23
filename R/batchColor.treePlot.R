#' Rboretum Basic Tree Plotter
#' 
#' Given a phylo object (tree), this ggtree wrapper returns a ggtree plot object, adjusted with possible arguments
#' @usage myPlot <- treePlotter(tree=myTree,...)
#' @param tree Phylo object
#' @param branch_length OPTIONAL: TRUE [plot tree with branch lengths]; FALSE [Default: plot cladogram]
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
#' @param color_legend TRUE [Include a legend for colored taxa/clades]; False [Default: No legend]
#' @param plot_titles OPTIONAL: Character vector of titles for plots
#' @return ggtree object
#' @export

batchColor.treePlot <- function(trees,to_color,branch_length,branch_weight,node_label,node_size,node_nudge,taxa_size,taxa_italic,taxa_align,taxa_offset,xmax,reverse_x,colors,color_legend,plot_titles){
  
  if(!Rboretum::is.multiPhylo(trees)){
    stop("'trees' does not appear to be a multiPhylo object. Use basic.treePlot() for single trees.")
  }
  
  if(missing(to_color)){
    stop("'to_color' argument required for batchColor.treePlot") 
  }
  else if(!is.character(to_color) && !is.list(to_color)){
    stop("'to_color' must be a character vector of taxa, or a list of taxa/clades")
  }
   
  if(missing(branch_length)){
    branch_length <- NA
  }
  
  if(missing(branch_weight)){
    branch_weight <- NA
  }
  
  if(missing(node_label)){
    node_label <- 'bs'
  }
  
  if(missing(node_size)){
    node_size <- NA
  }
  
  if(missing(node_nudge)){
    node_nudge <- NA
  }
    
  if(missing(taxa_size)){
    taxa_size <- NA
  }
  
  if(missing(taxa_italic)){
    taxa_italic <- NA
  }
  
  if(missing(taxa_align)){
    taxa_align <- NA
  }
  
  if(missing(taxa_offset)){
    taxa_offset <- NA
  }
  
  if(missing(xmax)){
    xmax <- NA
  }
  
  if(missing(reverse_x)){
    reverse_x <- NA
  }

  if(missing(colors)){
    colors <- NA
  }
  
  if(missing(color_legend)){
    color_legend <- NA
  }
  
  if(missing(plot_titles)){
    plot_titles <- rep(NA,length(trees))
  } else if(length(plot_titles) != length(trees)){
    print("Not enought plot titles for each tree, no titles passed.")
    plot_titles <- rep(NA,length(trees))
  }
  
  plotList <- list()
  
  for(i in 1:length(trees)){
    plotList[[i]] <- basic.treePlot(tree = trees[[i]],branch_length,branch_weight,node_label,node_size,node_nudge,taxa_size,taxa_italic,taxa_align,taxa_offset,xmax,reverse_x,to_color,colors,color_legend,plot_title=plot_titles[[i]])
  }
  
  tandem.treePlot(plotList)
}