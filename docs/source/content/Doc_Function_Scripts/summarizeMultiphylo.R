library(Rboretum)
sourceRboretum()

# Read in trees with two topologies (1=2, 1!=3)
myTree_1 <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))
myTree_2 <- readRooted(rb_tree4_path,root_taxa = c('Species_C','Species_H'))

myTrees <- c(myTree_1,myTree_2)

summarizeMultiPhylo(myTrees)
