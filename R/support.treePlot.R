#' Rboretum Site and Node Support Plotter
#' Given a phylo object (tree), and the output from tree.support(), this function plots the tree [adjusted with possible arguments] along with node support information
#' @param tree Rooted phylo object
#' @param tree_support  OPTIONAL: Output trom tree.support() run on the same tree
#' @param clade_support OPTIONAL: Output from compare.clades(). Will colorize nodepoint labels based on how many trees in multiPhylo contain that split
#' @param support_scales OPTIONAL: Scaling factor for nodepoint labels. Options include:
#' \itemize{
#'   \item "log" - Log-converted support values
#'   \item Single numeric value: All nodepoint labels will be this size [Default: 10]
#'   \item c(x,y): Support values will be re-scaled from numeric values x-y
#' }
#' @param node_alpha OPTIONAL: ggplot2 alpha value for geom_nodepoint [Default: 0.9]
#' @param node_color OPTIONAL: ggplot2 fill value for geom_nodepoint [Default: "red"]
#' @param legend_shape_size OPTIONAL: ggplot2 size for legend icons [Default: 10]
#' @param legend_font_size OPTIONAL: ggplot2 size for legend font [Default: 10]
#' @param legend_title_size OPTIONAL: ggplot size for legend title [Default: 10]
#' @param branch_length OPTIONAL: TRUE [plot tree with branch lengths]; FALSE [Default: plot cladogram]
#' @param branch_weight OPTIONAL: Set ggtree branch thickness
#' @param node_label OPTIONAL: Choice of node labels include:
#' \itemize{
#'   \item "none": No node labels [Default]
#'   \item "node": Node number
#'   \item "bs": Bootstrap value
#'   \item "support": Raw total support count
#' }
#' @param node_size OPTIONAL: Set ggtree node label size
#' @param node_nudge OPTIONAL: Set ggtree node label nudge_x 
#' @param taxa_size OPTIONAL: Set ggtree tip label size
#' @param taxa_italic OPTIONAL: TRUE [italic tip labels]; FALSE [Default: non-italic tip labels]
#' @param taxa_align OPTIONAL: 'left' or 'right' tip label alignment [Default: No alignment]
#' @param taxa_offset OPTIONAL: Set ggtree tip label offset
#' @param xmax OPTIONAL: Set ggplot xlim upper limit (e.g if long tip labels run off plot)
#' @param reverse_x OPTIONAL: TRUE [plot tree with tips on left]; FALSE [Default: plot tree with tips on right]
#' @param rename_tips OPTIONAL: Dataframe where column 1 contains tip labels for tree, and column 2 contains new, desired tip labels
#' @param plot_title OPTIONAL: Title for plot
#' @return ggtree object
#' @export
#' @examples
#' # Print tree with bootstrap labels
#' basic.treePlot(tree)

support.treePlot <- function(tree,tree_support,clade_support,support_scales,node_alpha,node_color,legend_shape_size,legend_font_size,legend_title_size,branch_length,branch_weight,node_label,node_size,node_nudge,taxa_size,taxa_italic,taxa_align,taxa_offset,xmax,reverse_x,rename_tips,plot_title){
  
  if(has_error(ape::is.rooted(tree))){
    stop("Error in ape::is.rooted. Is 'tree' a phylo object?")
  } else if(!ape::is.rooted(tree)){
    stop("Tree must be rooted for support.treePlot")}
  
  if(missing(tree_support)){
    tree_support <- Rboretum::get.splits(tree) %>%
      mutate(RBORETUM_DUMMY = 1)
    dummy_col <- TRUE
  } else if(!all(names(tree_support)[1:4] == c('Clade','Mirror_Clade','Split_Node','Split_Bootstrap'))){
    stop("'tree_support' argument must be output from tree.support()")
  } else if(ncol(tree_support)==4){
    tree_support$RBORETUM_DUMMY <- 1
    dummy_col <- TRUE
  } else{
    support_clades <- tree_support %>% pull(Clade) %>% as.character() %>% sort()
    tree_clades <- Rboretum::get.splits(tree) %>% pull(Clade) %>% as.character() %>% sort()
    if(all(support_clades == tree_clades) & all(tree_clades == support_clades)){
      support_cols <- 5:ncol(tree_support)
      dummy_col <- FALSE
    } else{
      stop("'tree' and 'tree_support' arguments contain different split information.")
    }
  }
  
  if(missing(clade_support)){
    if(dummy_col){ 
      stop("'tree_support' contained no signal columns and 'clade_support' missing. No data to map onto tree!") 
    } else{
      clade_support <- FALSE 
      tree_support$Tree_Count <- NA
    }
  } else{
    if(all(names(clade_support)==c('Clade','Tree_Count','Clade_Size','Tree_Percent','Trees_with_Clade'))){
      tree_clades <- Rboretum::get.clades(tree) %>% sort()
      support_clades <- clade_support$Clade %>% as.character() %>% sort()
      
      if(all(tree_clades %in% support_clades)){
        tree_support <- tree_support %>%
          mutate(Clade = as.character(Clade))
        clade_support <- clade_support %>%
          mutate(Clade = as.character(Clade))
        tree_support <- left_join(tree_support,select(clade_support,Clade,Tree_Count),by='Clade')
        clade_support <- TRUE
      } else{ stop("Clades from tree absent from 'clade_support'.")}
    } else { stop("'clade_support' argument must be output from compare.clades(return_shared_only=FALSE)") }
  }
  
  if(missing(support_scales)){
    support_scales <-  10
  } else if(is.character(support_scales)){
    if(support_scales == "log"){
      if(dummy_col){
        support_scales <- 10
      }
    } else { stop("Invalid argument for 'support_scales'") }
  } else if(is.numeric(support_scales)){
    if(length(support_scales) == 1){
      support_scales <- support_scales
    } else if (length(support_scales) == 2 & (support_scales[1] <= support_scales[2])){
      if(dummy_col){
        support_scales <- max(support_scales)
      }
    } else{ stop("Invalid argument for 'support_scales'") } 
  } else{ stop("Invalid argument for 'support_scales'") }
  
  if(missing(node_alpha)){
    node_alpha <- 0.9
  } else if(!is.numeric(node_alpha)){
    node_alpha <- 0.9
  }
  
  if(missing(node_color)){
    node_color <- "red"
  } else if(has_error(grDevices::col2rgb(node_color))){
    node_color <- "red"
  }
  
  if(missing(legend_shape_size)){
    legend_shape_size <- 10
  } else if(!is.numeric(legend_shape_size)){
    legend_shape_size <- 10
  }
  
  if(missing(legend_font_size)){
    legend_font_size <- 10
  } else if(!is.numeric(legend_font_size)){
    legend_font_size <- 10
  }
  
  if(missing(legend_title_size)){
    legend_title_size <- 10
  } else if(!is.numeric(legend_title_size)){
    legend_title_size <- 10
  }
  
  if(missing(branch_length)){
    branch_length <- FALSE
  } else{
    if(!is.logical(branch_length)){
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
    node_label <- 'none'
  } else if(!node_label %in%  c('bs','node','none','support')){
    node_label <- 'none'
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
    if(!is.logical(taxa_italic)){
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
  } else if(!is.logical(reverse_x)){
    reverseX <- FALSE
  } else{
    reverseX <- TRUE
  }
  
  if(missing(rename_tips)){
    renameTip <- FALSE
  } else if((is.data.frame(rename_tips) | is_tibble(rename_tips)) & length(names(rename_tips) >= 2)){
    old_id <- names(rename_tips)[1]
    new_id <- names(rename_tips)[2]
    if(has_error(Rboretum::convert.tips(tree,rename_tips,old_id,new_id))){
      renameTip <- FALSE
    } else{
      renameTip <- TRUE
    }
  } else{
    renameTip <- FALSE
  }
  
  if(missing(plot_title)){
    titlePlot <- FALSE
  } else if(is.na(plot_title)){
    titlePlot <- FALSE
  } else{
    plot_title <- as.character(plot_title)
    titlePlot <- TRUE
  }
  
  # Create ggtree_df
  
  if(!dummy_col){
    
    if(length(support_cols) == 1){
      tree_support$support <- tree_support[,5]
    } else{
      tree_support$support <- rowSums(tree_support[,support_cols])
    }
    
    if(all(tree_support$support == 0)){
      stop("All support values = 0")
    }
    
    non_zero_support <- tree_support %>% filter(support != 0)
    zero_support <- tree_support %>% filter(support == 0)
    
    if(length(support_scales)==1){
      if(support_scales == "log"){
        
        non_zero_support$scaled_support <- log(non_zero_support$support)
        
        if(nrow(zero_support)>0){
          zero_support$scaled_support <- 0
        }
        
        tree_support <- non_zero_support %>%
          rbind(zero_support) %>%
          arrange(Split_Node)
        
      } else if(is.numeric(support_scales)){
        
        non_zero_support$scaled_support <- support_scales
        
        if(nrow(zero_support)>0){
          zero_support$scaled_support <- 0
        }        
        
        tree_support <- non_zero_support %>%
          rbind(zero_support) %>%
          arrange(Split_Node)
        
      } else{ stop("Invalid argument for 'support_scales'") }
    } else if(is.numeric(support_scales) & length(support_scales) == 2){
      
      non_zero_support$scaled_support <- scales::rescale(non_zero_support$support,to=support_scales)

      if(nrow(zero_support)>0){
        zero_support$scaled_support <- 0
      }
      
      tree_support <- non_zero_support %>%
        rbind(zero_support) %>%
        arrange(Split_Node)
      
    } else{ stop("Invalid argument for 'support_scales'") }
    
    ggtree_df <- tree_support %>%
      select(Split_Node,Split_Bootstrap,support,scaled_support,Tree_Count) %>%
      rename(node = 'Split_Node', bootstrap = 'Split_Bootstrap', tree_count = 'Tree_Count') %>%
      filter(!is.na(node))
    
  } else{
    ggtree_df <- tree_support %>%
      select(Split_Node,Split_Bootstrap,Tree_Count) %>%
      mutate(support = NA,scaled_support = support_scales) %>%
      rename(node = 'Split_Node', bootstrap = 'Split_Bootstrap', tree_count = 'Tree_Count') %>%
      filter(!is.na(node))
  }
  
  # Create base tree
  
  if(renameTip){
    tree <- Rboretum::convert.tips(tree,rename_tips,old_id,new_id)
  }
  
  if(bWeight & branch_length){
    return_tree <- ggtree(tree,size=branch_weight) %<+% ggtree_df
  } else if(bWeight & !branch_length){
    return_tree <- ggtree(tree,size=branch_weight,branch.length = 'none') %<+% ggtree_df
  } else if(!bWeight & branch_length){
    return_tree <- ggtree(tree) %<+% ggtree_df
  } else if(!bWeight & !branch_length){
    return_tree <- ggtree(tree,branch.length = 'none') %<+% ggtree_df
  }
  
  if(extendX){
    if(!reverseX){
      return_tree <- return_tree + ggplot2::xlim(0,xmax)
    }
  }
  
  if(reverseX){
    return_tree <- return_tree + scale_x_reverse()
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
  
  if(!clade_support){
    return_tree <- return_tree + 
      geom_nodepoint(alpha=node_alpha,color=node_color,aes(size=scaled_support)) + 
      scale_size_identity()
    
  } else{
    return_tree <- return_tree + 
      geom_nodepoint(alpha=node_alpha,aes(size=scaled_support,color=as.integer(tree_count))) + 
      scale_size_identity() +
      scale_color_viridis(breaks = sort(unique(ggtree_df$tree_count)),name = "Trees with Split") +
      theme(legend.position="right",
            legend.title=element_text(size=legend_title_size), 
            legend.text=element_text(size=legend_font_size)) +
      guides(color = guide_legend(override.aes = list(size = legend_shape_size)))
  }
  
  # Process node labels
  if(node_label == "none"){
    return_tree <- return_tree
  } else if(node_label == "bs"){
    if(nNudge){
      if(nSize){
        return_tree <- return_tree + geom_nodelab(nudge_x = node_nudge,size=node_size,aes(label=bootstrap))
      } else{
        return_tree <- return_tree + geom_nodelab(nudge_x = node_nudge,aes(label=bootstrap))
      }
    } else{
      if(nSize){
        return_tree <- return_tree + geom_nodelab(size=node_size,aes(label=bootstrap))
      } else{
        return_tree <- return_tree + geom_nodelab(aes(label=bootstrap))
      }
    }
  } else if(node_label == "node"){
    if(nNudge){
      if(nSize){
        return_tree <- return_tree + geom_nodelab(nudge_x = node_nudge,size=node_size,aes(label=node))
      } else{
        return_tree <- return_tree + geom_nodelab(nudge_x = node_nudge,aes(label=node))
      }
    } else{
      if(nSize){
        return_tree <- return_tree + geom_nodelab(size=node_size,aes(label=node))
      } else{
        return_tree <- return_tree + geom_nodelab(aes(label=node))
      }
    }
  } else if(node_label == "support"){
    if(nNudge){
      if(nSize){
        return_tree <- return_tree + geom_nodelab(nudge_x = node_nudge,size=node_size,aes(label=support))
      } else{
        return_tree <- return_tree + geom_nodelab(nudge_x = node_nudge,aes(label=support))
      }
    } else{
      if(nSize){
        return_tree <- return_tree + geom_nodelab(size=node_size,aes(label=support))
      } else{
        return_tree <- return_tree + geom_nodelab(aes(label=support))
      }
    }
  }
  
  if(titlePlot){
    if(clade_support){
      return_tree <- return_tree + ggplot2::ggtitle(plot_title) + theme(legend.position="right",
                                                                        legend.title=element_text(size=legend_title_size),
                                                                        legend.text=element_text(size=legend_font_size),
                                                                        plot.title = element_text(hjust = 0.5))
    } else{
      return_tree <- return_tree + ggplot2::ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5))
    }
  }
  
  return(return_tree)
}