library(Rboretum)
sourceRboretum()

# Read in ultrametric timetrees (branch lengths ~ time)
rb_timeTree1_path
rb_timeTree2_path
rb_timeTree3_path

timeTree_1 <- readRooted(to_root = rb_timeTree1_path,root_taxa = c('Species_C','Species_H'))
timeTree_2 <- readRooted(to_root = rb_timeTree2_path,root_taxa = c('Species_C','Species_H'))
timeTree_3 <- readRooted(to_root = rb_timeTree3_path,root_taxa = c('Species_C','Species_H'))

# Create multiPhylo of trees
timeTrees <- c(timeTree_1,timeTree_2,timeTree_3)

extractNodeAges(timeTree_1)
extractNodeAges(timeTree_2)
extractNodeAges(timeTree_3)

print(extractNodeAges(timeTrees),n = 42)
  
extractNodeAges(timeTrees,return_summary = 'mean')
extractNodeAges(timeTrees,return_summary = 'median')
extractNodeAges(timeTrees,return_summary = 'both')
