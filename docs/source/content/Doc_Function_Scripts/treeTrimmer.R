library(Rboretum)
sourceRboretum()

# Read in a single tree rooted at the clade of Species C + Species H
is.rooted(ape::read.tree(rb_tree1_path))
myTree <- readRooted(to_root = rb_tree1_path, root_taxa = c('Species_C','Species_H'))
is.rooted(myTree)

# Trim myTree, keeping only Species A - I 
taxa_to_keep <- c("Species_A","Species_B","Species_C","Species_D","Species_E","Species_F","Species_G","Species_H","Species_I")
taxa_to_remove <- c("Species_J","Species_K","Species_L","Species_M","Species_N","Species_O")

# Trim myTree by supplying a list of taxa to keep
myTrimmedTree_keep <- treeTrimmer(tree = myTree,taxa = taxa_to_keep)

# Trim myTree by supplying a list of taxa to remove
myTrimmedTree_remove <- treeTrimmer(tree = myTree,taxa = taxa_to_remove,remove = TRUE)

# Check tip labels
naturalsort(myTrimmedTree_keep$tip.label)
naturalsort(myTrimmedTree_remove$tip.label)

# Read in all trees from a directory that start with "Gene" and end with ".nwk", and root at the clade of Species C + Species H
myTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),prefix="Gene",suffix=".nwk",dummy_names = TRUE)
purrr::map(.x=myTrees,.f=function(x){is.rooted(x)}) %>% unlist()

# Create a multiPhlyo where trees do not share all taxa (one tree only has Species_A - Species_I)
mixed_trees <- c(myTrees,myTrimmedTree_keep)

# Default behavior is to trim down to common taxa
myTrimmedTrees_mixed <- treeTrimmer(tree=mixed_trees)

# Check multiPhylo trimming
purrr::map(.x = myTrimmedTrees_mixed, .f = function(x){naturalsort(x$tip.label)})

# Trim all trees in the multiPhylo down to a given a list of taxa
myTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),prefix="Gene",suffix=".nwk",dummy_names = TRUE)
myTrimmedTrees <- treeTrimmer(tree=myTrees,taxa='Species_A;Species_B;Species_C;Species_H')

# Check multiPhylo trimming
purrr::map(.x = myTrimmedTrees, .f = function(x){naturalsort(x$tip.label)})
