library(Rboretum)
sourceRboretum()

# Date a phylo object where branch lengths are substitution rates
myTree <- readRooted(rb_tree1_path,c('Species_C','Species_H'))

# Estimate node ages by calibrating the root node to between 100MY and 120MY, iterating estimates 100 times
myDatedTree <- treeDater(tree = myTree, taxa='Species_C;Species_M',min_max = c(100,120),iterations = 100)

is.ultrametric(myDatedTree)
extractNodeAges(myDatedTree)

# Estimate node ages by calibrating at two nodes

myCalibration <- tibble(Taxon_A = c('Species_C','Species_A'),
                        Taxon_B = c('Species_M','Species_F'),
                        Min = c(100,15),
                        Max = c(120,17))
myCalibration

myRedatedTree <- treeDater(tree = myTree, calibration_df = myCalibration,iterations = 100)

is.ultrametric(myRedatedTree)
extractNodeAges(myRedatedTree)

# Date a multiPhylo object
myTrees <- readRooted(c(rb_tree1_path,rb_tree2_path,rb_timeTree3_path),c('Species_C','Species_H'))

myDatedTrees <- treeDater(tree = myTrees, taxa='Species_C;Species_M',min_max = c(100,120),iterations = 100)

all(is.ultrametric(myDatedTrees))
extractNodeAges(myDatedTrees,return_summary = 'both')
