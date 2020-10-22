library(Rboretum)
sourceRboretum()

# Read in a single phylogeny
myTree <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))

# Extract the signal from an alignment file
mySignal <- getAlignmentSignal(alignment_path = rb_align1_path,species_info = myTree)

# Map the signal from a single alignment onto a single tree topology
myTreeSupport <- getAlignmentSupport(tree = myTree,signal = mySignal,include_root = TRUE)

# Note the column name specifies how many taxa were allowed to have missing data (default: Total - 3)
myTreeSupport

# Map the signal from a single alignment onto a single tree topology, and only consider sites where one taxon has missing data. Add this to previous data.
myTreeSupport_m1 <- getAlignmentSupport(tree = myTree,signal = mySignal,include_root = TRUE,max_missing = 1,existing_support = myTreeSupport)
myTreeSupport_m1

# Map the signal from a multiple alignments onto a single tree topology
mySignals <- getAlignmentSignal(alignment_path = rb_all_align,species_info = myTree)

# Get tree support, separated by alignment. Treat gaps as missing data. Only include parsimony-informative sites
myTreeSupport_Sep <- getAlignmentSupport(tree = myTree,signal = mySignals,include_root = TRUE,include_gap = FALSE,only_parsinf = TRUE)

# Get same tree support, summed across alignments
myTreeSupport_Sum <- getAlignmentSupport(tree = myTree,signal = mySignals,include_root = TRUE,separate_signal = FALSE,include_gap = FALSE,only_parsinf = TRUE)

myTreeSupport_Sep
myTreeSupport_Sum

# Get support for a specific
myTreeSupport_Clade <- getAlignmentSupport(clade='Species_C;Species_H',signal = mySignals)
myTreeSupport_Clade

# Get support for a specific clade not found in any trees
myTreeSupport_Clade <- getAlignmentSupport(clade='Species_C;Species_F',signal = mySignals)
myTreeSupport_Clade
