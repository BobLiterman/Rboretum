library(Rboretum)
sourceRboretum()

# Read in multiPhylo
noname_multiPhylo <- readRooted(rb_all_unrooted,root_taxa = c('Species_C','Species_H'))

names(noname_multiPhylo) # Names are based on filenames by default

# Add dummy names
named_multiPhlyo <- treeNamer(noname_multiPhylo)
names(named_multiPhlyo)
