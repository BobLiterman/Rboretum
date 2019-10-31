#' Plot Pie Charts for Support from Different Datasets
#'
#' DESCRIPTION
#'
#' @param tree x
#' @param tree_support x
#' @param support_scales x
#' @param xmax x
#' @return x
#' @export
#' @examples
#' 
#'
plotTreeSupportBreakdown <- function(tree,tree_support,support_scales,xmax,plot_alpha,plot_type){
  
  if(missing(support_scales)){
    scale_values <- FALSE
  } else{
    if(is.numeric(support_scales)){
      if(length(support_scales) == 1){
        scale_values <- TRUE
        support_scales <- c(support_scales,support_scales)
      } else{
        if(length(support_scales) == 2){
          scale_values <- TRUE
          support_scales <- c(min(support_scales),max(support_scales))
        } else{
          stop('support_scales must be a one (fixed-size) or two (bounded limits) element numeric vector specifying minimum and maximum geom size')
        }
      } 
    } else{
      stop('support_scales must be a one (fixed-size) or two (bounded limits) element numeric vector specifying minimum and maximum geom size')
    }
  }
  
  if(missing(xmax)){
    rescale_x <- FALSE
  } else{
    if(is.numeric(xmax)){
      rescale_x <- TRUE
    } else{
      stop('xmax argument must be numeric')
    }
  }
  
  if(missing(plot_alpha)){
    plot_alpha <- 0.6
  }
  
  if(missing(plot_type)){
    plot_pies <- TRUE
    plot_bars <- FALSE
  } else{
    if(plot_type == 'pies'){
      plot_pies <- TRUE
      plot_bars <- FALSE 
    } else if (plot_type == 'bars'){
      plot_pies <- FALSE
      plot_bars <- TRUE
    } else{
      plot_pies <- TRUE
      plot_bars <- FALSE
    }
  }
  
  # Check that tree file matches support file
  if(any(is.na(tree_support$Split_Node))){
    tree_splits <- getAllSplits(rooted_tree = tree) %>% pull(Clade) %>% sort()
  } else{
    tree_splits <- getAllSplits(rooted_tree = tree,exclude_root = TRUE) %>% pull(Clade) %>% sort()
  }
  
  support_splits <- tree_support %>% pull(Clade) %>% sort()
  
  if(!(all(tree_splits == support_splits))){
    stop("Provided tree doesn't match support file.")
  }
  
  # Check that support file is from getTreeSupport
  tree_support_cols <- c("Clade","Mirror_Clade","Split_Node")
  multi_support_cols <- names(tree_support)
  
  if(!all(multi_support_cols[1:3] == tree_support_cols)){
    stop("'tree_support' argument doesn't appear to be the output from getTreeSupport (Check colnames")
  }
  if(ncol(tree_support)<=4){
    stop("This function requires multiple support columns. The supplied dataframe contains 1 or fewer support columns. To add support columns, re-run getTreeSupport with existing_splits argument.")
  }
  
  last_support <- ncol(tree_support)
  support_cols  <- 4:last_support
  dataset_count <- length(support_cols)
  support_ids <- multi_support_cols[support_cols]
  tree_support <- tree_support %>%
    mutate(Support = rowSums(.[4:last_support]))
  
  total_col <- last_support+1
  
  plot_df <- tree_support %>%
    rename(node = Split_Node) %>%
    filter(!is.na(node)) %>%
    arrange(node) %>%
    filter(Support != 0)
  
  if(rescale_x){
    blank_tree <- ggtree(tree,branch.length = "none") + geom_tiplab() + xlim(0,xmax)
  }
  else{
    blank_tree <- ggtree(tree,branch.length = "none") + geom_tiplab()
  }
  
  if(plot_pies){
    if(scale_values){
      plot_df <- plot_df %>%
        mutate(Support = rescale(Support,to=support_scales))
    }
    pies <- nodepie(plot_df,cols=support_cols,alpha = plot_alpha)
    inset(blank_tree,pies,height=plot_df$Support,width = plot_df$Support)
  } else{
    
    support_col <- plot_df %>% pull(Support) %>% as.numeric()

    for(x in support_cols){
      data_col <- plot_df %>% pull(x) %>% as.numeric()
      data_breakdown <- data_col/support_col
      data_median <- median(data_breakdown)
      data_expected <- support_col*data_median
      data_diff <- data_col-data_expected
      data_perc <- data_diff/data_expected
      plot_df[x] <- data_perc
    }
    
    bars <- nodebar(plot_df,cols=support_cols,alpha = plot_alpha,position="dodge")
    print(head(plot_df))
    inset(blank_tree,bars,height=3,width = 1)
  }
}