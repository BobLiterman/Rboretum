#' Plot Alignment Support and Node Prevalence for Specific Topology
#'
#' DESCRIPTION
#'
#' @param comparison x
#' @param signal x
#' @param tree x
#' @param support_scales x
#' @param xmax x
#' @return Tree with nodepoint geoms scaled to support, and colored by number of trees containing that split
#' @export
#' @examples
#' plotShuffledSupport(comparison,signal,tree,support_scales,xmax)
#'

plotShuffledSupport <- function(comparison,signal,tree,support_scales,xmax){

  if(missing(support_scales)){
    scale_values <- FALSE
  } else{
    if(length(support_scales) == 2 & is.numeric(support_scales)){
      scale_values <- TRUE
    } else{
      stop('support_scales must be a two element numeric vector specifying minimum and maximum geom size')
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

  # Get taxa from signal analysis
  signal_taxa <- signal %>%
    filter(!is.na(Split_1)) %>%
    head(1) %>%
    select(starts_with('Split_')) %>%
    select_if(~ !any(is.na(.))) %>%
    unite(col = "Taxa",sep = ";") %>%
    pull() %>% as.character() %>% str_split(pattern = ";") %>% unlist() %>% sort()

  # Get taxa from multiphylo comparison
  comparison_taxa <- comparison %>%
    pull(Clade) %>%
    paste(collapse = ';') %>%
    semiVector() %>%
    unique() %>%
    sort()

  # Get taxa from tree
  tree_taxa <- tree$tip.label %>% sort()

  # Ensure all taxa sets match
  if(!all(sapply(list(comparison_taxa, signal_taxa), FUN = identical, tree_taxa))){
    print("Tree Taxa:")
    print(tree_taxa)
    print("Signal Taxa:")
    print(signal_taxa)
    print("Comparison Taxa:")
    print(comparison_taxa)
    stop("Taxa from arguments don't match")
  }

  # Get all tree splits (excluding root)
  tree_splits <- Rboretum::getAllSplits(rooted_tree = tree) %>% select(-Mirror_Clade) %>%
    mutate(Clade = as.character(Clade)) %>%
    filter(!is.na(Split_Node))

  informative_patterns <- c('biallelic','gap_biallelic','triallelic','gap_triallelic','quadallelic','gap_quadallelic','pentallelic','gap_pentallelic')

  signal_counts <- signal %>%
    filter(Site_Pattern %in% informative_patterns) %>% select(starts_with('Split_')) %>%
    unlist() %>% table()

  alignment_name <- signal$Alignment_Name[1]

  pruned_comparison_splits <- comparison %>%
    filter(Clade %in% tree_splits$Clade) %>%
    left_join(tree_splits,by = "Clade")

  support_list <- c()
  for(clade in tree_splits$Clade){
    support_list <- c(support_list,tableCount(signal_counts,clade))
  }

  pruned_comparison_splits$Support <- as.integer(support_list)

  ggtree_df <- pruned_comparison_splits %>%
    select(Split_Node,Support,Tree_Percent,Tree_Count) %>%
    rename(node = Split_Node)

  if(scale_values){
    non_zero_support <- ggtree_df %>% filter(Support != 0)
    zero_support <-  ggtree_df %>% filter(Support == 0)
    ggtree_df <- non_zero_support %>%
      mutate(Support = rescale(Support,to=support_scales)) %>%
      rbind(zero_support)
  }

  blank_tree <- ggtree(tree,branch.length = "none") + geom_tiplab()

  if(rescale_x){
    supportTree <- blank_tree  %<+% ggtree_df +
      geom_nodepoint(alpha=0.9,aes(color = Tree_Count,size=Support)) +
      scale_size_identity() +
      scale_color_viridis(name = "Trees with Split") +
      xlim(0,xmax) +
      theme(legend.position="right") +
      guides(color = guide_legend())

  } else{
    supportTree <- blank_tree  %<+% ggtree_df +
      geom_nodepoint(alpha=0.9,aes(color = Tree_Count,size=Support)) +
      scale_size_identity() +
      scale_color_viridis(name = "Trees with Split") +
      theme(legend.position="right") +
      guides(color = guide_legend())
  }

  return(supportTree)
}
