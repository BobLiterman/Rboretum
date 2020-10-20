library(Rboretum)
sourceRboretum()

# Read in trees
myTrees <- readRooted(rb_unroot_dir,root_taxa = c('Species_C','Species_H'))
summarizeMultiPhylo(myTrees)

# Simulate trees with different numbers of taxa
mySimTrees <- c(rtree(10),rtree(10),rtree(20))
summarizeMultiPhylo(mySimTrees)
