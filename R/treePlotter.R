#' Rboretum Tree Plot Generator
#'
#' Given a phylo or multiPhylo object (tree), and some optional metadata, this ggtree wrapper returns a ggtree plot object, adjusted with possible arguments
#' @param tree Tree(s) to plot. Options include:
#' \itemize{
#'   \item A single, rooted phylo object; or,
#'   \item A rooted multiPhylo object where all trees share 3+ taxa [Only unique topologies plotted, from most to least common]
#' }
#' @param basic_plot OPTIONAL: If providing a multiPhylo and TRUE, plot all trees as stored and disable support functions [Default: FALSE, plot unique topologies after pruning to common taxa]
#' @param tree_support OPTIONAL: Output of getAlignmentSupport, including data from all clades in 'tree'
#' @param plot_root_support OPTIONAL: If TRUE and plotting node support via 'tree_support', include root nodes [Default: FALSE, do not include support for root nodes]
#' @param clade_support OPTIONAL: Output of getTreeClades(return_counts=TRUE), including data from all clades in 'tree'
#' @param geom_size OPTIONAL: If plotting support via tree_support, how should geom_nodepoint (or pies) be sized? Options include:	
#' \itemize{
#'   \item Single numeric value: All geom_nodepoint geoms will be this size [Default: 4]
#'   \item c(min,max): If used with tree_support, geoms will be sized based on the total support across datasets and rescaled between min and max.
#'   \item "log": If used with tree_support, geoms will be sized based on the log-transformed total support across datasets
#'   \item "none": Do not print geom_nodepoints
#' }
#' @param scale_range OPTIONAL: If plotting tree_support, supply a min and max of values to be scaled. Values outside this range will be displayed as the otherwise min/max geom size. [Default: scale all data]
#' @param use_pies OPTIONAL: If TRUE and tree_support contains inforomation from 2+ datasets, data from tree_support will be displayed as an inset pie chart rather than a geom_nodepoint [Default: FALSE, use geom_nodepoint; DEACTIVATED WHEN clade_support IS PROVIDED]	
#' @param pie_colors OPTIONAL If plotting with pies, a vector of colors with one item per support column [Default: ggplot default colors]
#' @param pie_scaler OPTIONAL: If plotting with pies, use pie_scaler to set a relative downward adjustment for geom size (ie. 0.1, 0.5, etc). If pie_scaler > 1, pies will be plotted without adjustment
#' @param pie_xnudge OPTIONAL: Set ggtree pie label hjust [Default: 0]	
#' @param pie_ynudge OPTIONAL: Set ggtree pie label yjust [Default: 0]	
#' @param pie_legend_position OPTIONAL: Numerical vector of length four specifying legend xmin,xmax,ymin,ymax respectively. [Default = c(1,1,1,1); Necessary oddity of generating a legend for nodepie data]
#' @param branch_length OPTIONAL: TRUE [plot tree(s) with branch lengths]; FALSE [Default: plot cladogram(s)]
#' @param branch_weight OPTIONAL: Set ggtree branch thickness [Default: 1]
#' @param node_label OPTIONAL: Choice of node labels include:
#' \itemize{
#'   \item "bs": Bootstrap value [Default]
#'   \item "none": No node labels
#'   \item "node": Node number
#'   \item "support": Raw total support count (Only with tree_support)
#'   \item "clade": Percent of trees in clade_support that support this clade (Only with clade_support)
#' }
#' @param node_label_font_size OPTIONAL: Set ggtree node label font size [Default: 5]
#' @param node_label_fontface OPTIONAL: Set node label fontface. Options include:
#' \itemize{
#'   \item "plain" [Default]
#'   \item "bold"
#'   \item "italic"
#'   \item "bold.italic"
#' }
#' @param node_label_color OPTIONAL: Set ggtree node label color [Default: 'black']
#' @param node_label_box OPTIONAL: Set ggtree node label as a label with a white background [Default: TRUE]
#' @param node_label_nudge OPTIONAL: Set ggtree node label nudge_x [Default: 0]
#' @param taxa_font_size OPTIONAL: Set ggtree tip label font size [Default: 5]
#' @param taxa_fontface OPTIONAL: Set tip label fontface. Options include:
#' \itemize{
#'   \item "plain" [Default]
#'   \item "bold"
#'   \item "italic"
#'   \item "bold.italic"
#' }
#' @param taxa_offset OPTIONAL: Set ggtree tip label offset [Default: 0]
#' @param xmax OPTIONAL: Set ggplot xlim upper limit (e.g if long tip labels run off plot)
#' @param xmin OPTIONAL: Set ggplot xlim lower limit (Must be used along with xmax [Default: 0])
#' @param reverse_x OPTIONAL: TRUE [plot tree with tips on left]; FALSE [Default: plot tree with tips on right]
#' @param to_color OPTIONAL: Color tips or clades via:
#' \itemize{
#'   \item Character vector of taxa, all to be colored the same color
#'   \item List of groups of taxa, each of which will have their own color. List can be named for use with a legend (set highlight_legend == TRUE)
#' }
#' @param colors OPTIONAL: Colors for clade highlighting. Must be hex or valid R colors. Provide a color for each group (1 if character vector, 1 for each group if named list) or default colors will be used.
#' @param highlight_legend OPTIONAL: Include a legend for colored tips, given a list (disabled if providing clade_support); [Default: False: No highlight legend]
#' @param color_branches OPTIONAL: If TRUE and coloring taxa or clades, color the branches rather than the tip labels [Default: FALSE, colorize tip labels; DEACTIVATED IF clade_support is provided]
#' @param plot_title OPTIONAL: Character vector containing plot titles (1 per tree) [Default: No title for phylo, tree name for multiPhylo]
#' @param plot_title_size OPTIONAL: Set ggplot title font size [Default: 14]
#' @param plot_title_fontface OPTIONAL: Set ggplot title fontface
#' \itemize{
#'   \item "plain"
#'   \item "bold" [Default]
#'   \item "italic"
#'   \item "bold.italic"
#' }
#' @param legend_shape_size OPTIONAL: ggplot2 size for legend icons [Default: 5]
#' @param legend_font_size OPTIONAL: ggplot2 size for legend font [Default: 10]
#' @param legend_title_size OPTIONAL: ggplot size for legend title [Default: 10]
#' @param geom_alpha OPTIONAL: ggplot2 alpha value for geom_nodepoint (or pies) [Default: 0.9]
#' @param geom_color OPTIONAL: ggplot2 color value for geom_nodepoint if clade_support not provided [Default: 'red']
#' @return ggtree object or list of ggtree objects
#' @export

treePlotter <- function(tree,basic_plot,tree_support,plot_root_support,clade_support,geom_size,scale_range,use_pies,pie_colors,pie_scaler,pie_xnudge,pie_ynudge,pie_legend_position,branch_length,branch_weight,node_label,node_label_font_size,node_label_fontface,node_label_color,node_label_box,node_label_nudge,taxa_font_size,taxa_fontface,taxa_offset,xmax,xmin,reverse_x,to_color,colors,highlight_legend,color_branches,plot_title,plot_title_size,plot_title_fontface,legend_shape_size,legend_font_size,legend_title_size,geom_alpha,geom_color){  
  
  # Ensure tree is valid for plotter
  if(!Rboretum::isMultiPhylo(tree) & !Rboretum::isPhylo(tree)){
    stop("'tree' must be either a phylo object or a mulitPhlyo object")
  } 
  
  # Check basic_plot
  if(missing(basic_plot)){
    basic_plot <- FALSE
  } else if(!is.logical(basic_plot)){
    basic_plot <- FALSE
  } else if(length(basic_plot)!=1){
    basic_plot <- FALSE
  } else if(Rboretum::isPhylo(tree)){
    basic_plot <- FALSE
  }
  
  if(!basic_plot & Rboretum::isMultiPhylo(tree) & !Rboretum::isMultiPhylo(tree,check_three_taxa=TRUE)){
    stop("multiPhylo objects must share 3+ taxa unless 'basic_plot' is set to TRUE...")
  }
  
  # Get tree root status and tree taxa
  if(Rboretum::isPhylo(tree)){ # If one tree is provided...
    
    root_status <- Rboretum::isPhylo(tree,check_rooted = TRUE)
    tree_taxa <- naturalsort(tree$tip.label)
    
  } else{ # If a multiPhylo is provided...
    
    # Add dummy tree names if necessary
    if(!Rboretum::isMultiPhylo(tree,check_named = TRUE)){
      tree <- Rboretum::treeNamer(tree)
    }
    
    raw_tree_names <- names(tree)
    root_status <- purrr::map(.x=tree,.f=function(x){Rboretum::isPhylo(x,check_rooted = TRUE)}) %>% unlist() %>% all()
    
    # If not doing basic plotting, reduce multiPhylo to common taxa
    if(!basic_plot){
      tree_taxa <- Rboretum::getSharedTaxa(tree)
      
      # If trimming of some trees is required, trim
      if(!Rboretum::isMultiPhylo(tree,check_all_taxa = TRUE)){
        tree <- Rboretum::treeTrimmer(tree)
      }
      names(tree) <-  raw_tree_names
    }
  }
  
  # Process trees and plot titles...
  if(Rboretum::isPhylo(tree)){ # If one tree is provided...
    
    tree_count <- 1
    
    if(root_status){
      tree_clades <- Rboretum::getTreeClades(tree)
    }
    
    # Default: No title for single trees
    if(missing(plot_title)){
      titlePlot <- FALSE
    } else if(!is.character(plot_title)){
      titlePlot <- FALSE
    } else if(length(plot_title)!=1){
      titlePlot <- FALSE
    } else{
      titlePlot <- TRUE
    }
    
  } else{ # If a multiPhylo is provided...
    
    # If simply plotting all trees, or if trees contain unrooted trees, collect relevant data
    if(basic_plot | !root_status){
      tree_count <- length(tree)
      tree_names <- names(tree)
      
      if(any(duplicated(tree_names))){
        stop("'tree' multiPhylo contains trees with identical names.")
      }
      
      titlePlot <- TRUE
      
      if(missing(plot_title)){
        plot_title <- tree_names
      } else if(!is.character(plot_title)){
        plot_title <- tree_names
      } else if(length(plot_title)!=tree_count){
        plot_title <- tree_names
      }
    } else{
      
      # If not plotting all trees, reduce multiPhylo to unique topologies after pruning to common taxa
      if(Rboretum::isMultiPhylo(tree,check_all_equal = TRUE)){
        
        tree_count <- 1
        tree <- Rboretum::getUniqueTopologies(tree,tree_names = TRUE)
        tree_clades <- Rboretum::getTreeClades(tree)
        
        # Default: No title for single trees
        if(missing(plot_title)){
          titlePlot <- FALSE
        } else if(!is.character(plot_title)){
          titlePlot <- FALSE
        } else if(length(plot_title)!=1){
          titlePlot <- FALSE
        } else{
          titlePlot <- TRUE
        }
        
      } else if(Rboretum::isMultiPhylo(tree,check_all_unique = TRUE)){ # If all trees are unique...
        
        tree_clades <- Rboretum::getTreeClades(tree)
        
        tree_count <- length(tree)
        tree_names <- names(tree)
        
        if(any(duplicated(tree_names))){
          stop("'tree' multiPhylo contains trees with identical names.")
        }
        
        titlePlot <- TRUE
        
        if(missing(plot_title)){
          plot_title <- tree_names
        } else if(!is.character(plot_title)){
          plot_title <- tree_names
        } else if(length(plot_title)!=tree_count){
          plot_title <- tree_names
        }
        
      } else if(Rboretum::isMultiPhylo(tree,check_some_equal = TRUE)){ # If trees contains duplicate topologies...
        
        tree_table <- Rboretum::getUniqueTopologies(tree,return_table = TRUE,tree_names = TRUE)
        tree <- Rboretum::getUniqueTopologies(tree,tree_names = TRUE)
        
        tree_clades <- Rboretum::getTreeClades(tree)
        
        tree_count <- length(tree)
        tree_names <- names(tree)
        
        if(any(duplicated(tree_names))){
          stop("'tree' multiPhylo contains trees with identical names.")
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
    }
  }
  
  # Is clade prevalence being provided?
  if(missing(clade_support) | basic_plot){
    
    cladeSupport <- FALSE
    
  } else if(!Rboretum::isCladeSupport(clade_support,tree,partial = TRUE)){
    stop("'clade_support' does not contain information about all the clades in 'tree'")
  } else if(!root_status){ # Don't process signal if unrooted trees are present
    cladeSupport <- FALSE
  } else{
    
    cladeSupport <- TRUE
    
    # Get breaks for legend
    clade_numbers <- unique(clade_support$Count) %>% sort()
    
    # Calculate percents for potential labels
    clade_tree_count <-  clade_support$Trees %>% semiVector() %>% unique() %>% length()
    clade_percent_char <- round(((as.numeric(clade_support$Count)/clade_tree_count)*100),digits = 1) %>% lapply(.,function(x) paste(c(x,"%"),collapse = '')) %>% unlist() %>% as.character()
    
    clade_support <- clade_support %>%
      select(Clade,Count) %>%
      rename(clade_count = 'Count') %>%
      mutate(Clade=as.character(Clade),clade_percent = as.character(clade_percent_char))
  }
  
  # Is alignment signal being mapped onto these trees?
  if(missing(tree_support) | basic_plot){
    
    treeSupport <- FALSE
    piePossible <- FALSE
  } else if(!Rboretum::isAlignmentSupport(tree_support,tree,partial = TRUE)){
    stop("'tree_support' does not contain information about all the clades in 'tree'")
  }  else if(!root_status){ # Don't process signal if unrooted trees are present
    treeSupport <- FALSE
    piePossible <- FALSE
  } else{
    
    treeSupport <- TRUE
    tree_support <- tree_support %>%
      filter(Clade %in% tree_clades)
    
    # Pie charts are only possible when tree_support contains 2+ columns of support	
    if(ncol(tree_support)>2 & !cladeSupport){	
      piePossible <- TRUE	
    } else{	
      piePossible <- FALSE # TODO: Fix for node labels	
    }
  }
  
  # Plot root support?
  if(missing(plot_root_support)){	
    plot_root_support <- FALSE	
  } else if(!is.logical(plot_root_support)){	
    plot_root_support <- FALSE	
  } else if(!plot_root_support){	
    plot_root_support <- FALSE	
  }	
  
  # Print pie charts on nodes?	
  if(missing(use_pies)){	
    use_pies <- FALSE	
  } else if(!is.logical(use_pies)){	
    use_pies <- FALSE	
  } else if(!piePossible){	
    use_pies <- FALSE	
  }	
  
  # Scale pies?	
  if(missing(pie_scaler)){	
    pie_scaler <- 1	
  } else if(!is.numeric(pie_scaler)){	
    pie_scaler <- 1
  }	else if(pie_scaler > 1){
    pie_scaler <- 1
  }
  
  # Nudge pies?	
  if(missing(pie_xnudge)){	
    pie_xnudge <- 0	
  } else if(!is.numeric(pie_xnudge)){	
    pie_xnudge <- 0	
  }	
  
  if(missing(pie_ynudge)){	
    pie_ynudge <- 0	
  } else if(!is.numeric(pie_ynudge)){	
    pie_ynudge <- 0	
  }
  
  # Only add nodepoint labels if tree_support or clade_support is provided
  if(!treeSupport & !cladeSupport){
    geom_size <- "none"
    nodePoint <- FALSE
  } else if(use_pies){
    nodePoint <- FALSE	
  } else{
    
    nodePoint <-  TRUE
    
    if(missing(geom_size)){
      geom_size <- 4
    } else if(length(geom_size)==1){
      if(is.character(geom_size)){
        if(!geom_size %in% c('log','none')){
          stop("If 'geom_size' is a character, only 'log' or 'none' are allowed")
        } else if(!geom_size == 'none'){
          if(treeSupport){
            geom_size <- 'log'
          } else{
            geom_size <- 4
          }
        } else{
          geom_size <- "none"
          nodePoint <- FALSE
        }
      } else if(is.numeric(geom_size)){
        geom_size <- as.numeric(geom_size)
      }
    } else if(length(geom_size)==2){
      
      if(!is.numeric(geom_size)){
        stop("'geom_size' is length 2, but not c(min,max)")
      }
      
      min_geom <- geom_size[1]
      max_geom <- geom_size[2]
      if(min_geom < max_geom){
        
        if(treeSupport){
          geom_size <- as.numeric(geom_size)
        } else{
          geom_size <- 4
        }
        
      } else{
        stop("'geom_size' is length 2, but not c(min,max)")
      }
    } else{
      stop("If 'geom_size' is numeric, it should be one or two items long.")
    }
  }
  
  # Print tree with branch lengths?
  if(missing(branch_length)){
    branch_length <- FALSE
  } else if(!is.logical(branch_length)){
    branch_length <- FALSE
  }
  
  # Set branch thickness
  if(missing(branch_weight)){
    branch_weight <- 1
  } else if(!is.numeric(branch_weight)){
    branch_weight <- 1
  }
  
  # Choose node label
  if(missing(node_label)){
    node_label <- 'bs'
  } else if(!is.character(node_label)){
    node_label <- 'bs'
  } else if(!(node_label %in% c('bs','node','none','support','clade'))){
    node_label <- 'bs'
  } else if(!treeSupport & node_label == 'support'){
    node_label <- 'bs'
  } else if(!cladeSupport & node_label == 'clade'){
    node_label <- 'bs'
  }
  
  # Choose node label font size
  if(missing(node_label_font_size)){
    node_label_font_size <- 5
  } else if(!is.numeric(node_label_font_size)){
    node_label_font_size <- 5
  }
  
  # Choose node label fontface
  if(missing(node_label_fontface)){
    node_label_fontface <- 'plain'
  } else if(!is.character(node_label_fontface)){
    node_label_fontface <- 'plain'
  } else if(!(node_label_fontface %in% c('plain','bold','italic','bold.italic'))){
    node_label_fontface <- 'plain'
  }
  
  # Color node label?
  if(missing(node_label_color)){
    node_label_color <- 'black'
  } else if(length(node_label_color)!= 1){
    node_label_color <- 'black'
  } else if(has_error(silent=TRUE,expr=grDevices::col2rgb(node_label_color))){
    node_label_color <- 'black'
  }
  
  # Print geom_nodelab in a box?
  if(missing(node_label_box)){
    node_label_box <- TRUE
  } else if(!is.logical(node_label_box)){
    node_label_box <- TRUE
  }
  
  # Nudge node label?
  if(missing(node_label_nudge)){
    node_label_nudge <- 0
  } else if(!is.numeric(node_label_nudge)){
    node_label_nudge <- 0
  }
  
  # Choose tip label font size
  if(missing(taxa_font_size)){
    taxa_font_size <- 5
  } else if(!is.numeric(taxa_font_size)){
    taxa_font_size <- 5
  }
  
  # Choose tip label fontface
  if(missing(taxa_fontface)){
    taxa_fontface <- 'plain'
  } else if(!is.character(taxa_fontface)){
    taxa_fontface <- 'plain'
  } else if(!(taxa_fontface %in% c('plain','bold','italic','bold.italic'))){
    taxa_fontface <- 'plain'
  }
  
  # Offset tip labels?
  if(missing(taxa_offset)){
    taxa_offset <- 0
  } else if(!is.numeric(taxa_offset)){
    taxa_offset <- 0
  }
  
  # Extend X-axis?
  if(missing(xmax) & missing(xmin)){
    extendX <- FALSE
  } else if(missing(xmax) & !missing(xmin)){
    extendX <- FALSE
  } else if(!is.numeric(xmax)){
    extendX <- FALSE
  } else if(missing(xmin)){
    extendX <- TRUE
    xmin <- 0
  } else if(!is.numeric(xmin)){
    extendX <- FALSE
  }  else if(xmin >= xmax){
    extendX <- FALSE
  } else{ extendX <- TRUE }
  
  # Reverse X-axis?
  if(missing(reverse_x)){
    reverse_x <- FALSE
  } else if(!is.logical(reverse_x)){
    reverse_x <- FALSE
  }
  
  # Set up distinct color list
  distinct_colors <- c('#e6194B','#4363d8','#9A6324','#f032e6','#800000','#000075','#f58231','#469990','#a9a9a9','#3cb44b','#42d4f4','#fabebe','#e6beff')
  
  # Color tips?
  if(missing(to_color)){
    to_color <- NULL
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
        colors <- c('#000000','#e6194B')
      } else{
        colors <- c('#000000',distinct_colors[1:group_count])
      }
    } else{
      if(group_count == 1){
        
        if(length(colors)>1){
          print('Too many highlight colors provided. Using first color in list...')
          colors <- colors[1]
        }
        
        if(has_error(silent=TRUE,expr=grDevices::col2rgb(colors))){
          print('Invalid color choice, using defaults...')
          colors <- c('#000000','#e6194B')
        } else{colors <- c('#000000',colors)}
        
      } else if(group_count > 1){
        if(length(colors) != group_count | any(has_error(silent=TRUE,expr=grDevices::col2rgb(colors)))){
          print("Invalid color choice. Using defaults.")
          colors <- c('#000000',distinct_colors[1:group_count])
        } else{colors <- c('#000000',colors)
        }
      }
    }
  }
  
  # Add color to tip labels or branches?
  if(missing(color_branches)){
    color_branches <- FALSE
  } else if(!is.logical(color_branches)){
    color_branches <- FALSE
  } else if(!colorTips){
    color_branches <- FALSE
  } else if(cladeSupport){
    color_branches <- FALSE
  }
  
  # Set legend icon shape size
  if(missing(legend_shape_size)){
    legend_shape_size <- 5
  } else if(!is.numeric(legend_shape_size)){
    legend_shape_size <- 5
  }
  
  # Set legend font size
  if(missing(legend_font_size)){
    legend_font_size <- 10
  } else if(!is.numeric(legend_font_size)){
    legend_font_size <- 10
  }
  
  # Set legend title size
  if(missing(legend_title_size)){
    legend_title_size <- 10
  } else if(!is.numeric(legend_title_size)){
    legend_title_size <- 10
  }
  
  # Set geom alpha value
  if(missing(geom_alpha)){
    geom_alpha <- 0.9
  } else if(!is.numeric(geom_alpha)){
    geom_alpha <- 0.9
  } else if(geom_alpha < 0 | geom_alpha > 1){
    geom_alpha <- 0.9
  }
  
  # Set geom color if tree support without clade support
  if(missing(geom_color)){
    geom_color <- '#e6194B'
  } else if(length(geom_color)>1){
    geom_color <- '#e6194B'
  } else if(has_error(silent=TRUE,expr=grDevices::col2rgb(geom_color))){
    geom_color <- '#e6194B'
  }
  
  # Handle legend arguments
  
  # Pie legend hack	
  if(missing(pie_legend_position)){	
    pie_legend_position <- c(1,1,1,1)	
  } else if(!is.numeric(pie_legend_position)){	
    pie_legend_position <- c(1,1,1,1)	
  } else if(length(pie_legend_position) != 4){	
    pie_legend_position <- c(1,1,1,1)	
  }	
  
  # Include legend if highlight from list of taxa?
  if(missing(highlight_legend)){
    highlight_legend <- FALSE
  } else if(!is.logical(highlight_legend)){
    highlight_legend <- FALSE
  } else if(!colorTips){
    highlight_legend <- FALSE
  } else if(group_count < 2){
    highlight_legend <- FALSE
  }
  
  # Rescale tree support
  if(treeSupport){
    raw_totals <- Rboretum::rescaleTreeSupport(tree_support,return_total = TRUE)
    
    if(missing(scale_range)){
      scaled_totals  <- Rboretum::rescaleTreeSupport(tree_support,scale = geom_size)
      tree_support_summary <- data.frame(Clade=as.character(tree_support$Clade),total_sites = as.integer(raw_totals),scaled_total = as.numeric(scaled_totals),pie_scales = pie_scaler*(scaled_totals/max(geom_size)),stringsAsFactors = FALSE)
    } else{
      scaled_totals  <- Rboretum::rescaleTreeSupport(tree_support,scale = geom_size,scale_range=scale_range)
      tree_support_summary <- data.frame(Clade=as.character(tree_support$Clade),total_sites = as.integer(raw_totals),scaled_total = as.numeric(scaled_totals),pie_scales = pie_scaler*(scaled_totals/max(geom_size)),stringsAsFactors = FALSE)
    }
  }
  
  # If plotting clade support, tip highlight legend is disabled
  if(cladeSupport){
    highlight_legend <- FALSE
  }
  
  # Create dummy multiPhylo for simple handling
  if(tree_count == 1){
    tree <- c(tree,tree)
  }
  
  # Create empty plot list
  plotList <- list()
  
  for(i in 1:tree_count){
    
    # Is this the last tree?
    last_tree <- ifelse(i==tree_count,TRUE,FALSE)
    
    # Pull out tree
    temp_tree <- tree[[i]]
    
    # Group for coloring tips/branches
    if(colorTips){
      temp_tree <- ggtree::groupOTU(temp_tree,to_color)
    }
    
    # Create ggtree_df if generating geom_nodepoint labels
    
    # Unrooted trees are plotted simply
    if(!root_status | basic_plot){
      ggtree_df <- tibble(node=integer(),Clade=character())
    } else{
      
      ggtree_df <- Rboretum::getTreeSplits(temp_tree) %>%
        select(Split_Node,Clade,Root) %>%
        rename(node=Split_Node) %>%
        mutate(node=as.integer(node), Clade = as.character(Clade))
      
      if(!plot_root_support){
        ggtree_df <- ggtree_df %>% filter(!Root) %>% select(-Root)
      } else{
        ggtree_df <- ggtree_df %>% select(-Root)
      }
    }
    
    if(cladeSupport & !treeSupport){
      
      ggtree_df <- ggtree_df %>%
        left_join(clade_support,by='Clade') %>%
        select(node,clade_count,clade_percent) %>%
        mutate(clade_count = as.factor(clade_count))
      
    } else if(treeSupport & !cladeSupport){
      
      ggtree_df <- ggtree_df %>%
        left_join(tree_support_summary,by='Clade') %>%
        select(node,total_sites,scaled_total) %>%
        rename(size_geom = 'scaled_total')
      
    } else if(treeSupport & cladeSupport){
      
      ggtree_df <- ggtree_df %>%
        left_join(clade_support,by='Clade') %>%
        left_join(tree_support_summary,by='Clade') %>%
        select(node,clade_count,clade_percent,total_sites,scaled_total) %>%
        rename(size_geom = 'scaled_total') %>%
        mutate(clade_count = as.factor(clade_count))
    }
    
    if(use_pies){	
      
      # Generate pie chart from tree support and make dummy legend	
      pie_df <- Rboretum::getTreeSplits(temp_tree) %>%	
        select(Split_Node,Clade,Root) %>%	
        rename(node=Split_Node) %>%	
        mutate(node=as.integer(node), Clade = as.character(Clade)) %>%	
        left_join(tree_support,by='Clade') %>% 	
        left_join(tree_support_summary,by='Clade') %>%	
        select(-Clade)
      
      if(!plot_root_support){
        pie_df <- pie_df %>% filter(!Root) %>% select(-Root)
      } else{
        pie_df <- pie_df %>% select(-Root)
      }
      
      support_counts <- length(2:ncol(tree_support))
      
      if(missing(pie_colors)){
        pies <- nodepie(pie_df,cols=2:ncol(tree_support),alpha = geom_alpha)
        
        pie_color_data <- cbind(tree_support[2:ncol(tree_support)]) %>%	
          gather(key = 'Dataset',value = 'Count') %>%	
          ggplot(aes(x=Dataset, y=Count,fill=Dataset)) +	
          geom_bar(stat="identity", width=1) +      	
          theme(legend.position="right",	
                legend.title=element_text(size=legend_title_size),	
                legend.text=element_text(size=legend_font_size)) +	
          guides(fill = guide_legend(override.aes = list(size = legend_shape_size)))	
      } else if(length(pie_colors) != support_counts){
        pies <- nodepie(pie_df,cols=2:ncol(tree_support),alpha = geom_alpha)
        
        pie_color_data <- cbind(tree_support[2:ncol(tree_support)]) %>%	
          gather(key = 'Dataset',value = 'Count') %>%	
          ggplot(aes(x=Dataset, y=Count,fill=Dataset)) +	
          geom_bar(stat="identity", width=1) +      	
          theme(legend.position="right",	
                legend.title=element_text(size=legend_title_size),	
                legend.text=element_text(size=legend_font_size)) +	
          guides(fill = guide_legend(override.aes = list(size = legend_shape_size)))
      } else if(any(has_error(silent=TRUE,expr=grDevices::col2rgb(pie_colors)))){
        pies <- nodepie(pie_df,cols=2:ncol(tree_support),alpha = geom_alpha)
        
        pie_color_data <- cbind(tree_support[2:ncol(tree_support)]) %>%	
          gather(key = 'Dataset',value = 'Count') %>%	
          ggplot(aes(x=Dataset, y=Count,fill=Dataset)) +	
          geom_bar(stat="identity", width=1) +      	
          theme(legend.position="right",	
                legend.title=element_text(size=legend_title_size),	
                legend.text=element_text(size=legend_font_size)) +	
          guides(fill = guide_legend(override.aes = list(size = legend_shape_size)))
      } else{
        pies <- nodepie(pie_df,cols=2:ncol(tree_support),alpha = geom_alpha,color=pie_colors)
        
        pie_color_data <- cbind(tree_support[2:ncol(tree_support)]) %>%	
          gather(key = 'Dataset',value = 'Count') %>%	
          ggplot(aes(x=Dataset, y=Count,fill=Dataset)) +	
          scale_fill_manual(values = pie_colors) +
          geom_bar(stat="identity", width=1) +      	
          theme(legend.position="right",	
                legend.title=element_text(size=legend_title_size),	
                legend.text=element_text(size=legend_font_size)) +	
          guides(fill = guide_legend(override.aes = list(size = legend_shape_size)))
      }
      
      legend_for_pie <- gtable::gtable_filter(ggplot_gtable(ggplot_build(pie_color_data)), "guide-box") 	
    }
    
    # Build base tree (Branch lengths? Colored branches?)
    
    # If group count = tip count, don't add black color label
    if(colorTips | color_branches){
      if(group_count == length(temp_tree$tip.label)){
        colors <- colors[-1]
      }
    }
    
    if(branch_length){
      if(colorTips & color_branches){
        
        if(is.character(to_color)){
          
          return_tree <- ggtree(temp_tree,size=branch_weight,aes(color=group)) %<+% ggtree_df
          return_tree <- return_tree + 
            scale_color_manual(values = colors,guide = (last_tree && highlight_legend))
          
        } else if(is.list(to_color)){
          return_tree <- ggtree(temp_tree,size=branch_weight,aes(color=group)) %<+% ggtree_df
          return_tree <- return_tree + 
            scale_color_manual("Focal Clades",breaks = names(to_color),values = colors,guide = (last_tree && highlight_legend))
        }
      } else{
        return_tree <- ggtree(temp_tree,size=branch_weight) %<+% ggtree_df
      }
    } else{
      if(colorTips & color_branches){
        
        if(is.character(to_color)){
          
          return_tree <- ggtree(temp_tree,branch.length = 'none',size=branch_weight,aes(color=group)) %<+% ggtree_df
          return_tree <- return_tree + 
            scale_color_manual(values = colors,guide = (last_tree && highlight_legend))
          
        } else if(is.list(to_color)){
          return_tree <- ggtree(temp_tree,branch.length = 'none',size=branch_weight,aes(color=group)) %<+% ggtree_df
          return_tree <- return_tree + 
            scale_color_manual("Focal Clades",breaks = names(to_color),values = colors,guide = (last_tree && highlight_legend))
        }
      } else{
        return_tree <- ggtree(temp_tree,branch.length = 'none',size=branch_weight) %<+% ggtree_df
      }
    }
    
    # Process tip labels
    if(!colorTips | color_branches){
      if(reverse_x){
        return_tree <- return_tree + geom_tiplab(size=taxa_font_size,fontface=taxa_fontface,offset = -taxa_offset,align=TRUE,linetype=NA,hjust=1)
      } else{
        return_tree <- return_tree + geom_tiplab(size=taxa_font_size,fontface=taxa_fontface,offset = taxa_offset)
      }
    } else{
      if(is.character(to_color)){
        if(reverse_x){
          return_tree <- return_tree + geom_tiplab(size=taxa_font_size,fontface=taxa_fontface,offset = -taxa_offset,aes(color=group),align=TRUE,linetype=NA,hjust=1) +
            scale_color_manual(values = colors,guide = (last_tree && highlight_legend))
        } else{
          return_tree <- return_tree + geom_tiplab(size=taxa_font_size,fontface=taxa_fontface,offset = taxa_offset,aes(color=group)) +
            scale_color_manual(values = colors,guide = (last_tree && highlight_legend))          
        }
        
      } else if(is.list(to_color)){
        if(reverse_x){
          return_tree <- return_tree + geom_tiplab(size=taxa_font_size,fontface=taxa_fontface,offset = -taxa_offset,aes(color=group),align=TRUE,linetype=NA,hjust=1) +
            scale_color_manual("Focal Clades",breaks = names(to_color),values = colors,guide = (last_tree && highlight_legend))
        } else{
          return_tree <- return_tree + geom_tiplab(size=taxa_font_size,fontface=taxa_fontface,offset = taxa_offset,aes(color=group)) +
            scale_color_manual("Focal Clades",breaks = names(to_color),values = colors,guide = (last_tree && highlight_legend))            
        }
        
      }
    }
    
    # Process geom_nodepoint
    if(nodePoint){
      if(treeSupport & !cladeSupport){
        return_tree <- return_tree + geom_nodepoint(color=geom_color, alpha=geom_alpha,aes(size=size_geom)) +
          scale_size_identity()
      } else if(cladeSupport & !treeSupport){
        if(colorTips){
          return_tree <- return_tree + 
            ggnewscale::new_scale_color() + 
            geom_nodepoint(size=geom_size, alpha=geom_alpha,aes(color=clade_count)) +
            scale_discrete_manual(aesthetics = c('color'),limits = factor(clade_numbers),values = viridisLite::viridis(length(clade_numbers)),name = "Trees with Split",guide = last_tree)
        } else{
          return_tree <- return_tree + geom_nodepoint(size=geom_size, alpha=geom_alpha,aes(color=clade_count)) +
            scale_discrete_manual(aesthetics = c('color'),limits = factor(clade_numbers),values = viridisLite::viridis(length(clade_numbers)),name = "Trees with Split",guide = last_tree)
        }
      } else if(cladeSupport & treeSupport){
        if(colorTips){
          return_tree <- return_tree + 
            ggnewscale::new_scale_color() + 
            geom_nodepoint(alpha=geom_alpha,aes(size=size_geom,color=clade_count)) +
            scale_size_identity() + 
            scale_discrete_manual(aesthetics = c('color'),limits = factor(clade_numbers),values = viridisLite::viridis(length(clade_numbers)),name = "Trees with Split",guide = last_tree)
        } else{
          return_tree <- return_tree + 
            geom_nodepoint(alpha=geom_alpha,aes(size=size_geom,color=clade_count)) +
            scale_size_identity() + 
            scale_discrete_manual(aesthetics = c('color'),limits = factor(clade_numbers),values = viridisLite::viridis(length(clade_numbers)),name = "Trees with Split",guide = last_tree)
        }
      }
    }
    
    # Reverse X-axis? Extend X-axis?
    
    if(reverse_x & extendX){
      return_tree <- return_tree + scale_x_reverse(limits = c(xmax,xmin))
    } else if(reverse_x & !extendX){
      return_tree <- return_tree + scale_x_reverse()
    } else if(!reverse_x & extendX){
      return_tree <- return_tree + ggplot2::xlim(xmin,xmax)
    }
    
    # Add titles
    
    # Set title size
    if(missing(plot_title_size)){
      plot_title_size <- 14
    } else if(!is.numeric(plot_title_size)){
      plot_title_size <- 14
    }
    
    # Set title fontface
    if(missing(plot_title_fontface)){
      plot_title_fontface <- 'bold'
    } else if(!is.character(plot_title_fontface)){
      plot_title_fontface <- 'bold'
    } else if(!(plot_title_fontface %in% c('plain','bold','italic','bold.italic'))){
      plot_title_fontface <- 'bold'
    }
    
    
    if(titlePlot){
      return_tree <- return_tree + ggplot2::ggtitle(plot_title[i])
    }
    
    # Adjust theme
    if(!last_tree){
      if(titlePlot){
        return_tree <- return_tree + theme(plot.title = element_text(hjust = 0.5,face = plot_title_fontface,size = plot_title_size))
      }
    } else{
      if(titlePlot){
        return_tree <- return_tree +
          theme(legend.position="right",
                legend.title=element_text(size=legend_title_size),
                legend.text=element_text(size=legend_font_size),
                plot.title = element_text(hjust = 0.5,face = plot_title_fontface,size = plot_title_size))
      } else{
        return_tree <- return_tree +
          theme(legend.position="right",
                legend.title=element_text(size=legend_title_size),
                legend.text=element_text(size=legend_font_size))
      }
      
      if(highlight_legend | cladeSupport){
        return_tree <- return_tree +
          guides(color = guide_legend(override.aes = list(size = legend_shape_size)))
      } else if(!cladeSupport & !highlight_legend){
        return_tree <- return_tree +
          guides(color = FALSE)
      }
    }
    
    # Add pies to all plots, but only add legend to last plot
    if(use_pies){	
      if(i!=tree_count){
        if(reverse_x){	
          return_tree <- inset(return_tree,pies,height=pie_df$pie_scales,width=pie_df$pie_scales,hjust=pie_xnudge,vjust=pie_ynudge,reverse_x = TRUE)	
        } else{	
          return_tree <- inset(return_tree,pies,height=pie_df$pie_scales,width=pie_df$pie_scales,hjust=pie_xnudge,vjust=pie_ynudge)	
        }	
      } else{
        if(reverse_x){	
          return_tree <- inset(return_tree,pies,height=pie_df$pie_scales,width=pie_df$pie_scales,hjust=pie_xnudge,vjust=pie_ynudge,reverse_x = TRUE) +
            ggplot2::annotation_custom(grob = legend_for_pie,xmin = pie_legend_position[1],xmax = pie_legend_position[2],ymin = pie_legend_position[3],ymax = pie_legend_position[4])	
          
        } else{	
          return_tree <- inset(return_tree,pies,height=pie_df$pie_scales,width=pie_df$pie_scales,hjust=pie_xnudge,vjust=pie_ynudge)	+
            ggplot2::annotation_custom(grob = legend_for_pie,xmin = pie_legend_position[1],xmax = pie_legend_position[2],ymin = pie_legend_position[3],ymax = pie_legend_position[4])
        }	
      }
    }
    
    # Process node labels
    if(node_label == "none"){
      return_tree <- return_tree
    } else if(node_label == "node"){
      if(node_label_box){
        return_tree <- return_tree + geom_nodelab(geom='label',fill='white',color=node_label_color,nudge_x = node_label_nudge,size=node_label_font_size,fontface=node_label_fontface,aes(label=node))
      } else{
        return_tree <- return_tree + geom_nodelab(color=node_label_color,nudge_x = node_label_nudge,size=node_label_font_size,fontface=node_label_fontface,aes(label=node))
      }
    } else if(node_label == "bs"){
      if(node_label_box){
        return_tree <- return_tree + geom_nodelab(geom='label',fill='white',color=node_label_color,nudge_x = node_label_nudge,size=node_label_font_size,fontface=node_label_fontface)
      } else{
        return_tree <- return_tree + geom_nodelab(color=node_label_color,nudge_x = node_label_nudge,size=node_label_font_size,fontface=node_label_fontface)
      }
    } else if(node_label == "clade"){
      if(node_label_box){
        return_tree <- return_tree + geom_nodelab(geom='label',fill='white',color=node_label_color,nudge_x = node_label_nudge,size=node_label_font_size,fontface=node_label_fontface,aes(label=clade_percent))
      } else{
        return_tree <- return_tree + geom_nodelab(color=node_label_color,nudge_x = node_label_nudge,size=node_label_font_size,fontface=node_label_fontface,aes(label=clade_percent))
      }
    } else if(node_label == "support"){
      if(node_label_box){
        return_tree <- return_tree + geom_nodelab(geom='label',fill='white',color=node_label_color,nudge_x = node_label_nudge,size=node_label_font_size,fontface=node_label_fontface,aes(label=total_sites))
      } else{
        return_tree <- return_tree + geom_nodelab(color=node_label_color,nudge_x = node_label_nudge,size=node_label_font_size,fontface=node_label_fontface,aes(label=total_sites))
      }
    }
    
    ###### FINISHED TREE; If a single tree is passed, return it; otherwise, add to plot list and continue...
    if(tree_count > 1){
      plotList[[i]] <- return_tree
    } else{
      return(return_tree)
    }
  }
  return(Rboretum::tandemPlotter(plotList))
}
