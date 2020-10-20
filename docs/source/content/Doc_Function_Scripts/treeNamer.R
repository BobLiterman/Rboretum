library(Rboretum)
sourceRboretum()

# Read in multiPhylo
myTrees <- readRooted(rb_all_unrooted,root_taxa = c('Species_C','Species_H'))

names(myTrees) # Names are based on filenames by default

# Remove names from trees
names(myTrees) <- NULL
names(myTrees)

# Add dummy names
named_multiPhlyo <- treeNamer(myTrees)
names(named_multiPhlyo)
