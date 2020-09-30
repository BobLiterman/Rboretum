library(Rboretum)
sourceRboretum()

# Read in trees with two topologies (1=2, 1!=3)
myTree_1 <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))
myTree_2 <- readRooted(rb_tree2_path,root_taxa = c('Species_C','Species_H'))
myTree_3 <- readRooted(rb_tree4_path,root_taxa = c('Species_C','Species_H'))

mySameTrees <- c(myTree_1,myTree_2)
myDifferentTrees <- c(myTree_1,myTree_3)

# Topology A
getTreeSplits(myTree_1)
getTreeSplits(myTree_2)

# Topology B
getTreeSplits(myTree_3)

# Get splits when trees are identical
getTreeSplits(mySameTrees)

# Get splits when trees are not identical
getTreeSplits(myDifferentTrees)
