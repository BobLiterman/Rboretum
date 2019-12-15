#' Rboretum Tree Plot Generator
#'
#' Given a phylo or multiPhylo object (tree), and some optional metadata, this ggtree wrapper returns a ggtree plot object, adjusted with possible arguments
#' @param tree Tree(s) to plot. Options include:
#' \itemize{
#'   \item A single, rooted phylo object; or,
#'   \item A rooted multiPhylo object where all trees share 3+ taxa [Only unique topologies plotted, from most to least common]
#' }
#' @param clade_support OPTIONAL: Output of getTreeClades(return_counts=TRUE), including data from all clades in 'tree'
#' @param tree_support OPTIONAL: Output of getTreeSupport, including data from all clades in 'tree'
#' @param geom_size OPTIONAL: How should geom_nodepoint (or pies) be sized? Options include:
#' \itemize{
#'   \item Single numeric value: All geom_nodepoint geoms will be this size [Default: 4]
#'   \item c(min,max): If used with tree_support, geoms will be sized based on the total support across datasets and rescaled between min and max.
#'   \item "log": If used with tree_support, geoms will be sized based on the log-transformed total support across datasets
#' }
#' @param use_pies OPTIONAL: If TRUE and tree_support contains inforomation from 2+ datasets, data from tree_support will be displayed as an inset pie chart rather than a geom_nodepoint [Default: FALSE, use geom_nodepoint; DEACTIVATED WHEN clade_support IS PROVIDED]
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
#' @param node_label_nudge OPTIONAL: Set ggtree node label nudge_x [Default: 0]
#' @param taxa_font_size OPTIONAL: Set ggtree tip label font size [Default: 5]
#' @param taxa_fontface OPTIONAL: Set tip label fontface. Options include:
#' \itemize{
#'   \item "plain" [Default]
#'   \item "bold"
#'   \item "italic"
#'   \item "bold.italic"
#' }
#' @param taxa_offset OPTIONAL: Set ggtree tip label offset
#' @param xmax OPTIONAL: Set ggplot xlim upper limit (e.g if long tip labels run off plot)
#' @param reverse_x OPTIONAL: TRUE [plot tree with tips on left]; FALSE [Default: plot tree with tips on right]
#' @param to_color OPTIONAL: Color tips or clades via:
#' \itemize{
#'   \item Character vector of taxa, all to be colored the same color
#'   \item List of groups of taxa, each of which will have their own color. List can be named for use with a legend (set highlight_legend == TRUE)
#' }
#' @param colors OPTIONAL: Colors for clade highlighting. Must be hex or valid R colors. Provide a color for each group (1 if character vector, 1 for each group if named list) or default colors will be used.
#' @param highlight_legend OPTIONAL: Include a legend for colored tips, given a list; [Default: False: No highlight legend]
#' @param color_branches OPTIONAL: If TRUE and coloring taxa or clades, color the branches rather than the tip labels [Default: FALSE, colorize tip labels; DEACTIVATED IF clade_support is provided]
#' @param plot_title OPTIONAL: Character vector containing plot titles (1 per tree) [Default: No title for phylo, tree name for multiPhylo]
#' @param legend_shape_size OPTIONAL: ggplot2 size for legend icons [Default: 5]
#' @param legend_font_size OPTIONAL: ggplot2 size for legend font [Default: 10]
#' @param legend_title_size OPTIONAL: ggplot size for legend title [Default: 10]
#' @param geom_alpha OPTIONAL: ggplot2 alpha value for geom_nodepoint (or pies) [Default: 0.9]
#' @param geom_color OPTIONAL: ggplot2 color value for geom_nodepoint if clade_support not provided [Default: 'red']
#' @return ggtree object or list of ggtree objects
#' @export

treePlotter <- function(tree,clade_support,tree_support,geom_size,use_pies,pie_xnudge,pie_ynudge,pie_legend_position,branch_length,branch_weight,node_label,node_label_font_size,node_label_fontface,node_label_nudge,taxa_font_size,taxa_fontface,taxa_offset,xmax,reverse_x,to_color,colors,highlight_legend,color_branches,plot_title,legend_shape_size,legend_font_size,legend_title_size,geom_alpha,geom_color){
  
  # Ensure tree is valid for plotter
  if(!Rboretum::isMultiPhylo(tree,check_rooted = TRUE,check_three_taxa = TRUE) & !Rboretum::isPhylo(tree,check_rooted = TRUE)){
    stop("'tree' must be either a rooted phylo object or a rooted, mulitPhlyo object basicPlotter")
  }
  
  # Process trees and plot titles...
  if(Rboretum::isPhylo(tree)){ # If one tree is provided...
    
    tree_taxa <- sort(tree$tip.label)
    tree_clades <- Rboretum::getTreeClades(tree)
    tree_count <- 1
    
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
    
    # Add dummy tree names if necessary
    if(!Rboretum::isMultiPhylo(tree,check_named = TRUE)){
      tree <- Rboretum::treeNamer(tree)
    }
    
    if(Rboretum::isMultiPhylo(tree,check_all_equal = TRUE)){
      
      tree_taxa <- Rboretum::getSharedTaxa(tree)
      tree <- Rboretum::treeTrimmer(tree[[1]],tree_taxa)
      tree_clades <- Rboretum::getTreeClades(tree)
      tree_count <- 1
      
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
      
      tree_taxa <- Rboretum::getSharedTaxa(tree)
      tree <- Rboretum::treeTrimmer(tree,tree_taxa)
      tree_clades <- Rboretum::getTreeClades(tree)
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
      
    } else if(Rboretum::isMultiPhylo(tree,check_some_equal = TRUE)){ # If trees contains duplicate topologies...
      
      tree_taxa <- Rboretum::getSharedTaxa(tree)
      tree_table <- Rboretum::getUniqueTopologies(tree,return_table = TRUE)
      tree <- Rboretum::getUniqueTopologies(tree)
      tree_clades <- Rboretum::getTreeClades(tree)
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
  }
  
  # Is clade prevalence being provided?
  if(missing(clade_support)){
    
    cladeSupport <- FALSE
    
  } else if(!Rboretum::isCladeSupport(clade_support,tree,partial = TRUE)){
    stop("'clade_support' does not contain information about all the clades in 'tree'")
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
  if(missing(tree_support)){
    
    treeSupport <- FALSE
    piePossible <- FALSE
    
  } else if(!Rboretum::isTreeSupport(tree_support,tree_clades,partial = TRUE)){
    stop("'tree_support' does not contain information about all the clades in 'tree'")
  } else{
    
    treeSupport <- TRUE
    tree_support <- tree_support %>%
      filter(Clade %in% tree_clades)
    
    # Pie charts are only possible when tree_support contains 2+ columns of support
    if(ncol(tree_support)>2 & !cladeSupport){
      piePossible <- TRUE
    } else{
      piePossible <- FALSE
    }
  }
  
  # Print pie charts on nodes?
  if(missing(use_pies)){
    use_pies <- FALSE
  } else if(!is.logical(use_pies)){
    use_pies <- FALSE
  } else if(!piePossible){
    use_pies <- FALSE
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
    nodePoint <- FALSE
  } else if(use_pies){
    nodePoint <- FALSE
  } else{
    nodePoint <-  TRUE
  }
  
  # Process geom_size
  if(missing(geom_size)){
    geom_size <- 4
    
  } else if(is.character(geom_size)){
    
    if(length(geom_size)!=1){
      stop("If 'geom_size' is a character, only 'log' is allowed")
    } else if(geom_size != 'log'){
      stop("If 'geom_size' is a character, only 'log' is allowed")
    } else{
      
      if(treeSupport){
        geom_size <- 'log'
      } else{
        geom_size <- 4
      }
    }
    
  } else if(is.numeric(geom_size)){
    
    if(length(geom_size)==1){
      geom_size <- as.numeric(geom_size)
      
    } else if(length(geom_size)==2){
      
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
  
  # Print tree with branch lenghts?
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
    tOffset <- FALSE
  } else if(!is.numeric(taxa_offset)){
    tOffset<-FALSE
  } else{ tOffset <- TRUE }
  
  # Extend X-axis?
  if(missing(xmax)){
    extendX <- FALSE
  } else if(!is.numeric(xmax)){
    extendX<-FALSE
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
    scaled_totals  <- Rboretum::rescaleTreeSupport(tree_support,scale = geom_size)
    tree_support_summary <- data.frame(Clade=as.character(tree_support$Clade),total_sites = as.integer(raw_totals),scaled_total = as.numeric(scaled_totals)) %>%
      mutate(Clade=as.character(Clade))
  }
  
  # Create dummy multiPhylo for simple handling
  if(tree_count == 1){
    tree <- c(tree,tree)
  }
  
  # Create empty plot list
  plotList <- list()
  
  print(colors)
  for(i in 1:tree_count){
    
    # Pull out tree
    temp_tree <- tree[[i]]
    
    # Group for coloring tips/branches
    if(colorTips){
      temp_tree <- ggtree::groupOTU(temp_tree,to_color)
    }
    
    # Create ggtree_df if generating geom_nodepoint labels
    ggtree_df <- Rboretum::getTreeSplits(temp_tree) %>%
      filter(!is.na(Split_Node)) %>%
      select(Split_Node,Clade) %>%
      rename(node=Split_Node) %>%
      mutate(node=as.integer(node), Clade = as.character(Clade))
    
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
        filter(!is.na(Split_Node)) %>%
        select(Split_Node,Clade) %>%
        rename(node=Split_Node) %>%
        mutate(node=as.integer(node), Clade = as.character(Clade)) %>%
        left_join(tree_support,by='Clade') %>% 
        left_join(tree_support_summary,by='Clade') %>%
        select(-Clade)
      
      pies <- nodepie(pie_df,cols=2:ncol(tree_support),alpha = geom_alpha)
      
      pie_color_data <- cbind(tree_support[2:ncol(tree_support)]) %>%
        gather(key = 'Dataset',value = 'Count') %>%
        ggplot(aes(x=Dataset, y=Count,fill=Dataset)) +
        geom_bar(stat="identity", width=1) +      
        theme(legend.position="right",
              legend.title=element_text(size=legend_title_size),
              legend.text=element_text(size=legend_font_size)) +
        guides(fill = guide_legend(override.aes = list(size = legend_shape_size)))
      
      legend_for_pie <- gtable::gtable_filter(ggplot_gtable(ggplot_build(pie_color_data)), "guide-box") 
    }
    
    # Build base tree (Branch lengths? Colored branches?)
    
    if(branch_length){
      if(colorTips & color_branches){
        
        if(is.character(to_color)){
          
          return_tree <- ggtree(temp_tree,size=branch_weight,aes(color=group),show.legend=FALSE) %<+% ggtree_df
          return_tree <- return_tree + 
            scale_color_manual(values = colors)
          
        } else if(is.list(to_color)){
          
          return_tree <- ggtree(temp_tree,size=branch_weight,aes(color=group)) %<+% ggtree_df
          return_tree <- return_tree + 
            scale_color_manual("Focal Clades",breaks = names(to_color),values = colors)
        }
      } else{
        return_tree <- ggtree(temp_tree,size=branch_weight) %<+% ggtree_df
      }
    } else{
      if(colorTips & color_branches){
        
        if(is.character(to_color)){
          
          return_tree <- ggtree(temp_tree,branch.length = 'none',size=branch_weight,aes(color=group),show.legend=FALSE) %<+% ggtree_df
          return_tree <- return_tree + 
            scale_color_manual(values = colors)
          
        } else if(is.list(to_color)){
          
          return_tree <- ggtree(temp_tree,branch.length = 'none',size=branch_weight,aes(color=group)) %<+% ggtree_df
          return_tree <- return_tree + 
            scale_color_manual("Focal Clades",breaks = names(to_color),values = colors)
        }
      } else{
        return_tree <- ggtree(temp_tree,branch.length = 'none',size=branch_weight) %<+% ggtree_df
      }
    }
    
    # Process tip labels
    if(!colorTips | color_branches){
      return_tree <- return_tree + geom_tiplab(size=taxa_font_size,fontface=taxa_fontface,offset = taxa_offset)
    } else{
      if(is.character(to_color)){
        return_tree <- return_tree + geom_tiplab(size=taxa_font_size,fontface=taxa_fontface,offset = taxa_offset,aes(color=group),show.legend=FALSE) +
          scale_color_manual(values = colors)
      } else if(is.list(to_color)){
        return_tree <- return_tree + geom_tiplab(size=taxa_font_size,fontface=taxa_fontface,offset = taxa_offset,aes(color=group)) +
          scale_color_manual("Focal Clades",breaks = names(to_color),values = colors)
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
            Rboretum::new_scale_color() + 
            geom_nodepoint(size=geom_size, alpha=geom_alpha,aes(color=clade_count)) +
            scale_discrete_manual(aesthetics = c('color'),limits = clade_numbers,values = viridisLite::viridis(length(clade_numbers)),name = "Trees with Split")
        } else{
          return_tree <- return_tree + geom_nodepoint(size=geom_size, alpha=geom_alpha,aes(color=clade_count)) +
            scale_discrete_manual(aesthetics = c('color'),limits = clade_numbers,values = viridisLite::viridis(length(clade_numbers)),name = "Trees with Split")
        }
      } else if(cladeSupport & treeSupport){
        if(colorTips){
          return_tree <- return_tree + 
            Rboretum::new_scale_color() + 
            geom_nodepoint(alpha=geom_alpha,aes(size=size_geom,color=clade_count)) +
            scale_size_identity() + 
            scale_discrete_manual(aesthetics = c('color'),limits = clade_numbers,values = viridisLite::viridis(length(clade_numbers)),name = "Trees with Split")
        } else{
          return_tree <- return_tree + 
            geom_nodepoint(alpha=geom_alpha,aes(size=size_geom,color=clade_count)) +
            scale_size_identity() + 
            scale_discrete_manual(aesthetics = c('color'),limits = clade_numbers,values = viridisLite::viridis(length(clade_numbers)),name = "Trees with Split")
        }
      }
    }
    
    # Reverse X-axis?
    if(reverse_x){
      return_tree <- return_tree + scale_x_reverse()
    }
    
    # Extend X-axis?
    if(extendX){
      if(!reverse_x){
        return_tree <- return_tree + ggplot2::xlim(0,xmax)
      }
    }
    
    # Add titles
    if(titlePlot){
      return_tree <- return_tree + ggplot2::ggtitle(plot_title[i])
    }
    
    # Adjust theme
    if(i!=tree_count){
      if(titlePlot){
        return_tree <- return_tree + theme(plot.title = element_text(hjust = 0.5))
      }
    } else{
      if(titlePlot){
        return_tree <- return_tree +
          theme(legend.position="right",
                legend.title=element_text(size=legend_title_size),
                legend.text=element_text(size=legend_font_size),
                plot.title = element_text(hjust = 0.5)) +
          guides(color = guide_legend(override.aes = list(size = legend_shape_size)))
        
      } else{
        return_tree <- return_tree +
          theme(legend.position="right",
                legend.title=element_text(size=legend_title_size),
                legend.text=element_text(size=legend_font_size)) +
          guides(color = guide_legend(override.aes = list(size = legend_shape_size)))
        
      }
    }
    
    if(use_pies){
      if(temp_pie_legend){
        return_tree <- return_tree + 
          ggplot2::annotation_custom(grob = legend_for_pie,xmin = pie_legend_position[1],xmax = pie_legend_position[2],ymin = pie_legend_position[3],ymax = pie_legend_position[4]) 
      }
      
      if(reverse_x){
        return_tree <- inset(return_tree,pies,height=pie_df$scaled_total,width=pie_df$scaled_total,hjust=pie_xnudge,vjust=pie_ynudge,reverse_x = TRUE)
      } else{
        return_tree <- inset(return_tree,pies,height=pie_df$scaled_total,width=pie_df$scaled_total,hjust=pie_xnudge,vjust=pie_ynudge)
      }
    }
    
    # Process node labels
    if(node_label == "none"){
      return_tree <- return_tree
    } else if(node_label == "node"){
      return_tree <- return_tree + geom_nodelab(nudge_x = node_label_nudge,size=node_label_font_size,fontface=node_label_fontface,aes(label=node))
    } else if(node_label == "bs"){
      return_tree <- return_tree + geom_nodelab(nudge_x = node_label_nudge,size=node_label_font_size,fontface=node_label_fontface)
    } else if(node_label == "clade"){
      return_tree <- return_tree + geom_nodelab(nudge_x = node_label_nudge,size=node_label_font_size,fontface=node_label_fontface,aes(label=clade_percent))
    } else if(node_label == "support"){
      return_tree <- return_tree + geom_nodelab(nudge_x = node_label_nudge,size=node_label_font_size,fontface=node_label_fontface,aes(label=total_sites))
    }
    
    ###### FINISHED TREE; If a single tree is passed, return it; otherwise, add to plot list and continue...
    if(tree_count > 1){
      plotList[[i]] <- return_tree
    } else{
      return(return_tree)
    }
  }
  
  return(plotList)
}