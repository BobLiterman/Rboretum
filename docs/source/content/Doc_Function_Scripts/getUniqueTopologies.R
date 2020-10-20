library(Rboretum)
sourceRboretum()

# Read in multiPhylo 
myMultiPhylo <- readRooted(to_root = rb_all_unrooted,root_taxa = c('Species_C','Species_H'))

# Reduce to unique topologies, use tree names for new trees; Print results and return unique trees
myUniqueTrees_1 <- getUniqueTopologies(trees = myMultiPhylo,tree_names = TRUE,print_table = TRUE)

# Reduce to unique topologies, use dummy names for new trees; Print results and return unique trees
myUniqueTrees_2 <- getUniqueTopologies(trees = myMultiPhylo,tree_names = FALSE,print_table = TRUE)

# Reduce to unique topologies, use tree names for new trees; Return results table
myUniqueTrees_3 <- getUniqueTopologies(trees = myMultiPhylo,tree_names = TRUE,return_table = TRUE)

myUniqueTrees_1
myUniqueTrees_2
myUniqueTrees_3