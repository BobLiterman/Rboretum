library(Rboretum)

# Set test data directory
sourceRboretum()

# Read in a single tree rooted at the clade of Species C + Species H
myTree <- readRooted(to_root = rb_tree1_path, root_taxa = c('Species_C','Species_H'))

# Trim myTree down to Species A - I   
taxa_to_keep <- c("Species_A","Species_B","Species_C","Species_D","Species_E","Species_F","Species_G","Species_H","Species_I")
taxa_to_remove <- c("Species_J","Species_K","Species_L","Species_M","Species_N","Species_O")

# Trim myTree by supplying a list of taxa to keep
myTrimmedTree_keep <- treeTrimmer(tree = myTree,taxa = taxa_to_keep)

# Trim myTree by supplying a list of taxa to remove
myTrimmedTree_remove <- treeTrimmer(tree = myTree,taxa = taxa_to_remove,remove = TRUE)

# Check tip labels
naturalsort(myTrimmedTree_keep$tip.label)

naturalsort(myTrimmedTree_remove$tip.label)

# Read in a a multiPhylo of trees rooted at the clade of Species C + Species H
myTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),suffix=".nwk")

# Make a multiPhlyo where trees do not share all taxa
mixed_trees <- c(myTrees,myTrimmedTree_keep)

# Default behavior is to trim down to common taxa
myTrimmedTrees_mixed <- treeTrimmer(tree=mixed_trees)

# Check multiPhylo trimming
purrr::map(.x = myTrimmedTrees_mixed, .f = function(x){naturalsort(x$tip.label)})

# Trim multiPhylo given a list of taxa
myTrimmedTrees <- treeTrimmer(tree=myTrees,taxa=taxa_to_keep)

# Check multiPhylo trimming
purrr::map(.x = myTrimmedTrees, .f = function(x){naturalsort(x$tip.label)})
