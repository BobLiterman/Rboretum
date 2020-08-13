library(Rboretum)

# Set test data directory
sourceRboretum()

# View raw tree
raw_tree <- ape::read.tree(rb_tree1_path)
treePlotter(tree=raw_tree,xmax=10)

# Read in a single tree, rooted at the clade of Species C + Species H
myTree <- readRooted(to_root = rb_tree1_path, root_taxa = c('Species_C','Species_H'))
treePlotter(tree=myTree,xmax=10)

# From a directory containing multiple trees, read in all '.nwk' files and root at the clade of Species C + Species H
myTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),suffix=".nwk")
treePlotter(tree=myTrees,xmax=10)

# Same as above, but add user-defined tree tree_names as opposed to tree file basenames
myTreeNames <- c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E')
myNamedTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),suffix=".nwk",tree_names=myTreeNames)
treePlotter(tree=myTreeNames,xmax=10)

# Same as above, but add placeholder tree_names ('Tree_1' - 'Tree_5') as opposed to tree file basenames
myDummyTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),suffix=".nwk",dummy_names=TRUE)
treePlotter(tree=myDummyTrees,xmax=10)

names(myTrees)
names(myNamedTrees)
names(myDummyTrees)
