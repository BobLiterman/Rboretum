library(Rboretum)
sourceRboretum()

# Read in a single tree and root at the clade of Species C + Species H
raw_tree <- ape::read.tree(rb_tree1_path)
ape::is.rooted(raw_tree)

myTree <- readRooted(to_root = rb_tree1_path, root_taxa = c('Species_C','Species_H'))
ape::is.rooted(myTree)

# Read in a multiple unrooted trees and root at the clade of Species C + Species H
purrr::map(.x=rb_all_unrooted,.f=function(x){ape::read.tree(x) %>% ape::is.rooted(.)}) %>% unlist()

myTrees <- readRooted(to_root = rb_all_unrooted, root_taxa = c('Species_C','Species_H'))
purrr::map(.x=myTrees,.f=function(x){ape::is.rooted(x)}) %>% unlist()

# From a directory containing multiple unrooted trees, read in all '.nwk' files and root at the clade of Species C + Species H
myTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),prefix="Gene",suffix=".nwk")
purrr::map(.x=myTrees,.f=function(x){ape::is.rooted(x)}) %>% unlist()

# Same as above, but add placeholder tree_names ('Tree_1' - 'Tree_5') as opposed to tree file basenames
myDummyTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),prefix="Gene",suffix=".nwk",dummy_names = TRUE)
purrr::map(.x=myDummyTrees,.f=function(x){ape::is.rooted(x)}) %>% unlist()

# Same as above, but add user-defined tree tree_names as opposed to tree file basenames
myTreeNames <- c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E')
myNamedTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),prefix="Gene",suffix=".nwk",tree_names=myTreeNames)
purrr::map(.x=myNamedTrees,.f=function(x){ape::is.rooted(x)}) %>% unlist()