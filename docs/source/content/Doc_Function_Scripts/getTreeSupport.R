library(Rboretum)
sourceRboretum()

# Read in a single phylogeny
myTree <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))

# Extract the signal from an alignment file
mySignal <- getAlignmentSignal(alignment_path = rb_align1_path,species_info = myTree)

# Map the signal from a single alignment onto a single tree topology
myTreeSupport <- getTreeSupport(tree = myTree,signal = mySignal,include_root = TRUE)

myTreeSupport

# Map the signal from a multiple alignments onto a single tree topology
mySignals <- getAlignmentSignal(alignment_path = rb_all_align,species_info = myTree)

# Get tree support, separated by alignment
myTreeSupport_Sep <- getTreeSupport(tree = myTree,signal = mySignals,include_root = TRUE)

# Get tree support, summed across alignments
myTreeSupport_Sum <- getTreeSupport(tree = myTree,signal = mySignals,include_root = TRUE,separate_signal = FALSE)

myTreeSupport_Sep
myTreeSupport_Sum

# Get support for a specific clade not found in any trees
myTreeSupport_Clade <- getTreeSupport(clade='Species_C;Species_F',signal = mySignals)
myTreeSupport_Clade
