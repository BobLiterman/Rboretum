library(Rboretum)
sourceRboretum()

# Read in ultrametric timetrees (branch lengths ~ time)
timeTree_1 <- readRooted(to_root = rb_timeTree1_path,root_taxa = c('Species_C','Species_H'))
timeTree_2 <- readRooted(to_root = rb_timeTree2_path,root_taxa = c('Species_C','Species_H'))
timeTree_3 <- readRooted(to_root = rb_timeTree3_path,root_taxa = c('Species_C','Species_H'))
timeTree_4 <- readRooted(to_root = rb_timeTree4_path,root_taxa = c('Species_C','Species_H'))
timeTree_5 <- readRooted(to_root = rb_timeTree5_path,root_taxa = c('Species_C','Species_H'))

extractNodeAges(timeTree_1)
extractNodeAges(timeTree_2)
extractNodeAges(timeTree_3)
extractNodeAges(timeTree_4)
extractNodeAges(timeTree_5)

# Create multiPhylo of trees that share a topology
timeTrees <- c(timeTree_1,timeTree_2,timeTree_3) %>% treeNamer()

extractNodeAges(timeTrees) %>% arrange(Clade)
extractNodeAges(timeTrees,return_summary = 'mean')
extractNodeAges(timeTrees,return_summary = 'median')
extractNodeAges(timeTrees,return_summary = 'both')
