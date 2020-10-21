library(Rboretum)
sourceRboretum()

# Read in ultrametric timetrees (branch lengths ~ time)
timeTree_1 <- readRooted(to_root = rb_timeTree1_path,root_taxa = c('Species_C','Species_H'))

extractNodeAges(timeTree_1)

# Create multiPhylo of trees that share a topology
timeTrees <- c(timeTree_1,timeTree_2,timeTree_3) %>% treeNamer()

extractNodeAges(timeTrees) %>% arrange(Clade)
extractNodeAges(timeTrees,return_summary = 'mean')
extractNodeAges(timeTrees,return_summary = 'median')
extractNodeAges(timeTrees,return_summary = 'both')
