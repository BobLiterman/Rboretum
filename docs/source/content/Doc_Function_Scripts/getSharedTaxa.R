library(Rboretum)
sourceRboretum()

# Create a multiPhylo where all trees share all taxa

tree_1 <- ape::rtree(25)
tree_2 <- ape::rtree(25)
tree_3 <- ape::rtree(25)

trees <- c(tree_1,tree_2,tree_3)

getSharedTaxa(trees)
length(getSharedTaxa(trees))

# Create a multiPhylo where all trees share 10 taxa

tree_1 <- ape::rtree(30)
tree_2 <- ape::rtree(20)
tree_3 <- ape::rtree(10)

getSharedTaxa(c(tree_1,tree_2,tree_3))
length(getSharedTaxa(c(tree_1,tree_2,tree_3)))