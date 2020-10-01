library(Rboretum)
sourceRboretum()

# Read in multiPhylo 
myMultiPhylo <- readRooted(to_root = rb_all_unrooted,root_taxa = c('Species_C','Species_H'))
myMultiPhylo
treePlotter(myMultiPhylo,basic_plot = TRUE,xmin=-0.5,xmax=10,taxa_font_size = 3,node_label_font_size = 3)

# Reduce to unique topologies, use tree names for new trees
myUniqueTrees_1 <- getUniqueTopologies(trees = myMultiPhylo,tree_names = TRUE,print_table = TRUE)
myUniqueTrees_1
treePlotter(myUniqueTrees_1,basic_plot = TRUE,xmin=-0.5,xmax=10,taxa_font_size = 3,node_label_font_size = 3)

# Reduce to unique topologies, use dummy names for new trees
myUniqueTrees_2 <- getUniqueTopologies(trees = myMultiPhylo,tree_names = FALSE,print_table = TRUE)
myUniqueTrees_2
treePlotter(myUniqueTrees_2,basic_plot = TRUE,xmin=-0.5,xmax=10,taxa_font_size = 3,node_label_font_size = 3)
