#' Basic Tree Plotter
#'
#' This function plots a tree without branch lengths, and nodes can be labeled with bootstrap or node labels
#' @param tree_1 Phylo object
#' @param label Either 'node' [print node labels] or 'bootstrap' [print bootstrap labels]. Any other entry (or missing) results in no node labels
#' @export
#' @examples
#' # Print tree with no labels
#' basicTreePlot(tree)
#'
#' # Print tree with bootstrap labels
#' basicTreePlot(tree,'bootstrap')
#'
#' # Print tree with node labels
#' basicTreePlot(tree,'node')
#'

# Plot trees with shared nodes highlighted
sharedNodePlotter <- function(tree,split_df){

  test_names <- c("Clade","Tree_Count","Clade_Size")
  if(!(all(names(split_df) == test_names))){
    stop("Split table provided is not the output of Rboretum script.")
  }

  tree_species_list <- sort(tree$tip.label)
  split_species_list <- split_df %>% pull(Clade) %>% c() %>% str_split(pattern = ";") %>% unlist() %>% unique() %>% sort()

  if(!all.equal(tree_species_list,split_species_list)){
    stop("Tree and split table contain different species.")
  }

  shared_splits <- shared_split_plot[[2]] %>% pull(Clade)
  other_splits <- shared_split_plot[[1]] %>% filter(Clade_Count >= min_count) %>% pull(Clade)

  shared_node_list <- c()
  for(f in shared_splits){
    clade <- unlist(str_split(f,pattern = ","))
    shared_node_list <- c(shared_node_list,getMRCA(tree,clade))
  }
  shared_color_list <- rep('red',each=length(shared_node_list))

  other_node_list <-  c()
  for(g in other_splits){
    clade <- unlist(str_split(g,pattern = ","))
    if(is.monophyletic(tree,clade)){
      other_node_list <- c(other_node_list,getMRCA(tree,clade))
    }
  }
  other_color_list <- rep('black',each=length(other_node_list))

  ggtree_df <- data.frame('node'=c(shared_node_list,other_node_list),'Color'=c(shared_color_list,other_color_list))

  plot(ggtree(tree,branch.length = 'none')  %<+% ggtree_df +
         geom_tiplab() +
         geom_nodelab() +
         geom_nodepoint(alpha=0.5,aes(color=Color,fill=Color,size=4)) +
         scale_color_identity() +
         scale_fill_identity())
}
