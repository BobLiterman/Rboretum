library(Rboretum)

# Set test data directory
sourceRboretum()

# Read in a single tree, rooted at the clade of Species C + Species H
myTree <- readRooted(to_root = rb_tree1_path, root_taxa = c('Species_C','Species_H'))

# From a directory containing multiple trees, read in all '.nwk' files and root at the clade of Species C + Species H
myTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),suffix=".nwk")

# Same as above, but add user-defined tree tree_names as opposed to tree file basenames
myTreeNames <- c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E')
myNamedTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),suffix=".nwk",tree_names=myTreeNames)

# Same as above, but add placeholder tree_names ('Tree_1' - 'Tree_5') as opposed to tree file basenames
myDummyTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),suffix=".nwk",dummy_names=TRUE)

names(myTrees)
names(myNamedTrees)
names(myDummyTrees)
