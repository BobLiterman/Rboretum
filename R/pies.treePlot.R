#' Rboretum Pie Graph Node Support Plotter
#' 
#' Given a phylo object (tree), and the output from tree.support(), this function plots the tree [adjusted with possible arguments] along with node support information
#' @param tree Rooted phylo object
#' @param tree_support  Output trom tree.support() run on the same tree
#' @param support_scales OPTIONAL: Scaling factor for nodepoint labels. Options include:
#' \itemize{
#'   \item "log" - Log-converted support values
#'   \item Single numeric value: All nodepoint labels will be this size [Default: 3]
#'   \item c(x,y): Support values will be re-scaled from numeric values x-y
#' }
#' @param node_alpha OPTIONAL: ggplot2 alpha value for geom_nodepoint [Default: 0.9]
#' @param legend_shape_size OPTIONAL: ggplot2 size for legend icons [Default: 10]
#' @param legend_font_size OPTIONAL: ggplot2 size for legend font [Default: 10]
#' @param legend_title_size OPTIONAL: ggplot size for legend title [Default: 10]
#' @param legend_position OPTIONAL: Numerical vector of length four specifying legend xmin,xmax,ymin,ymax respectively. [Necessary oddity of generating a legend for nodepie data]
#' @param node_label OPTIONAL: Choice of node labels include:
#' \itemize{
#'   \item "none": No node labels [Default]
#'   \item "node": Node number
#'   \item "bs": Bootstrap value
#'   \item "support": Raw total support count
#' }
#' @param node_size OPTIONAL: Set ggtree node label size
#' @param node_nudge OPTIONAL: Set ggtree node label nudge_x 
#' @param pie_xnudge OPTIONAL: Set ggtree pie label hjust
#' @param pie_ynudge OPTIONAL: Set ggtree pie label yjust
#' @param branch_length OPTIONAL: TRUE [plot tree with branch lengths]; FALSE [Default: plot cladogram]
#' @param branch_weight OPTIONAL: Set ggtree branch thickness
#' @param taxa_size OPTIONAL: Set ggtree tip label size
#' @param taxa_italic OPTIONAL: TRUE [italic tip labels]; FALSE [Default: non-italic tip labels]
#' @param taxa_align OPTIONAL: 'left' or 'right' tip label alignment [Default: No alignment]
#' @param taxa_offset OPTIONAL: Set ggtree tip label offset
#' @param xmax OPTIONAL: Set ggplot xlim upper limit (e.g if long tip labels run off plot)
#' @param rename_tips OPTIONAL: Dataframe where column 1 contains tip labels for tree, and column 2 contains new, desired tip labels
#' @return ggtree object
#' @export
#' @examples
#' # Print tree with bootstrap labels
#' basic.treePlot(tree)

pies.treePlot <- function(tree,tree_support,support_scales,node_alpha,legend_shape_size,legend_font_size,legend_title_size,legend_position,branch_length,branch_weight,node_label,node_size,node_nudge,pie_xnudge,pie_ynudge,taxa_size,taxa_italic,taxa_align,taxa_offset,xmax,rename_tips){
  
  if(has_error(ape::is.rooted(tree))){
    stop("Error in ape::is.rooted. Is 'tree' a phylo object?")
  } else if(!ape::is.rooted(tree)){
    stop("Tree must be rooted for support.treePlot")}
  
  if(missing(tree_support)){
    stop("tree_support required for pie charts")
  } else if(!all(names(tree_support)[1:4] == c('Clade','Mirror_Clade','Split_Node','Split_Bootstrap'))){
    stop("'tree_support' argument must be output from tree.support()")
  } else if(ncol(tree_support)<6){
    stop("tree_support required for 2+ alignments to make a pie chart")
  } else{
    support_clades <- tree_support %>% pull(Clade) %>% as.character() %>% sort()
    tree_clades <- Rboretum::get.splits(tree) %>% pull(Clade) %>% as.character() %>% sort()
    if(all(support_clades == tree_clades) & all(tree_clades == support_clades)){
      support_cols <- 5:ncol(tree_support)
    } else{
      stop("'tree' and 'tree_support' arguments contain different split information.")
    }
  }
  
  if(missing(support_scales)){
    support_scales <-  3
  } else if(is.character(support_scales)){
    if(support_scales == "log"){
      support_scales <- support_scales
    } else { stop("Invalid argument for 'support_scales'") }
  } else if(is.numeric(support_scales)){
    if(length(support_scales) == 1){
      support_scales <- support_scales
    } else if (length(support_scales) == 2 & (support_scales[1] <= support_scales[2])){
      support_scales <- support_scales
    } else{ stop("Invalid argument for 'support_scales'") } 
  } else{ stop("Invalid argument for 'support_scales'") }
  
  if(missing(node_alpha)){
    node_alpha <- 0.9
  } else if(!is.numeric(node_alpha)){
    node_alpha <- 0.9
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
  
  if(missing(legend_position)){
    legend_position <- FALSE
  } else if(!is.numeric(legend_position)){
    legend_position <- FALSE
  } else if(length(legend_position) != 4){
    legend_position <- FALSE
  } else{
    leg_xmin <- legend_position[1]
    leg_xmax <- legend_position[2]
    leg_ymin <- legend_position[3]
    leg_ymax <- legend_position[4]
    legend_position <- TRUE
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
  
  if(missing(pie_xnudge)){
    pxNudge <- FALSE
  } else{
    if(!is.numeric(pie_xnudge)){
      pxNudge<-FALSE
    } else{ pxNudge <- TRUE }
  }
  
  if(missing(pie_ynudge)){
    pyNudge <- FALSE
  } else{
    if(!is.numeric(pie_ynudge)){
      pyNudge<-FALSE
    } else{ pyNudge <- TRUE }
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

  # Create ggtree_df
  
  tree_support$support <- rowSums(tree_support[,support_cols])
  
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
        arrange(Split_Node)} 
    
    else if(is.numeric(support_scales)){
      
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
    select(Split_Node,Split_Bootstrap,support,scaled_support) %>%
    rename(node = 'Split_Node', bootstrap = 'Split_Bootstrap') %>%
    filter(!is.na(node))
  
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
      return_tree <- return_tree + ggplot2::xlim(0,xmax)
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
  
  tree_support <- tree_support %>% 
    filter(!is.na(Split_Node)) %>%
    rename(node='Split_Node')
  
  pies <- nodepie(tree_support,cols=support_cols,alpha = node_alpha)
  
  pie_color_data <- cbind(tree_support[support_cols]) %>%
    gather(key = 'Dataset',value = 'Count') %>%
    ggplot(aes(x=Dataset, y=Count,fill=Dataset)) +
    geom_bar(stat="identity", width=1) +      
    theme(legend.title.align=0.5,legend.position="right",legend.title=element_text(size=legend_title_size),legend.text=element_text(size=legend_font_size)) +
    guides(color = guide_legend(override.aes = list(size = legend_shape_size)))
  
  pie_legend <- gtable::gtable_filter(ggplot_gtable(ggplot_build(pie_color_data)), "guide-box") 
  
  if(legend_position){
    return_tree <- return_tree + ggplot2::annotation_custom(grob = pie_legend,xmin = leg_xmin,xmax = leg_xmax,ymin = leg_ymin,ymax = leg_ymax) 
  }
  
  if(pxNudge){
    if(pyNudge){
      return(inset(return_tree,pies,height=tree_support$scaled_support,width = tree_support$scaled_support,hjust=pie_xnudge,vjust=pie_ynudge))
    } else{ return(inset(return_tree,pies,height=tree_support$scaled_support,width = tree_support$scaled_support,hjust=pie_xnudge)) }
  } else{
    if(pyNudge){
      return(inset(return_tree,pies,height=tree_support$scaled_support,width = tree_support$scaled_support,vjust=pie_ynudge))
    } else {return(inset(return_tree,pies,height=tree_support$scaled_support,width = tree_support$scaled_support)) }
  }
}