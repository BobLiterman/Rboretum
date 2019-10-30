#' Basic Tree Plotter
#'
#' This function plots a tree without branch lengths, and nodes can be labeled with bootstrap or node labels
#' @param tree Phylo object
#' @param tree_support Output from getTreeSupport
#' @param support_scales OPTIONAL: Vector of values to rescale support geoms
#' @export
#' @examples
#' # Set up data
#' tree< - myTree
#' signal <- getAlignmentSignal(alignment_path,species_info = myTree)
#' tree_support <- getTreeSupport(tree,signal)
#'
#' supportPlot(tree,tree_support)
#'

supportPlot <- function(tree,tree_support,support_scales,xmax){

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

  ggtree_df <- tree_support %>%
    select(-Clade,-Mirror_Clade) %>%
    `colnames<-`(c('node','Support')) %>%
    filter(!is.na(node)) %>%
    arrange(node)

  if(scale_values){
    non_zero_support <- ggtree_df %>% filter(Support != 0)
    zero_support <-  ggtree_df %>% filter(Support == 0)
    ggtree_df <- non_zero_support %>%
      mutate(Support = rescale(Support,to=support_scales)) %>%
      rbind(zero_support) %>%
      arrange(node)
  }

  if(rescale_x){
    blank_tree <- ggtree(tree,branch.length = "none") + geom_tiplab() + xlim(0,xmax)
  }
  else{
    blank_tree <- ggtree(tree,branch.length = "none") + geom_tiplab()
  }
  
  supportTree <- blank_tree  %<+% ggtree_df +
    geom_nodepoint(color="red",alpha=0.9,aes(size=Support)) +
    scale_size_identity()

  return(supportTree)
}
