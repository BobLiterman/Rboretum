library(Rboretum)
sourceRboretum()

# Read in trees with two topologies (1=2, 1!=3)
myTree_1 <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))
myTree_2 <- readRooted(rb_tree2_path,root_taxa = c('Species_C','Species_H'))
myTree_3 <- readRooted(rb_tree4_path,root_taxa = c('Species_C','Species_H'))

# Topology A
getTreeClades(myTree_1)
getTreeClades(myTree_2)

# Topology B
getTreeClades(myTree_3)

# Get splits from a multiPhylo
myMultiPhylo <- readRooted(c(rb_tree1_path,rb_tree2_path,rb_tree4_path),root_taxa = c('Species_C','Species_H'))

# Return all identified splits
getTreeClades(myMultiPhylo)
length(getTreeClades(myMultiPhylo))

# Return all identified splits, but exclude root split
getTreeClades(myMultiPhylo,include_root = FALSE)
length(getTreeClades(myMultiPhylo,include_root = FALSE))

# Return only splits present in all trees
getTreeClades(myMultiPhylo,return_shared = TRUE)
length(getTreeClades(myMultiPhylo,return_shared = TRUE))

# Return a table of results
getTreeClades(myMultiPhylo,return_counts = TRUE)
