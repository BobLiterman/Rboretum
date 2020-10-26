library(Rboretum)
sourceRboretum()

# Read in ultrametric timetrees (branch lengths ~ time)
timeTrees <- readRooted(c(rb_timeTree1_path,rb_timeTree2_path,rb_timeTree3_path),root_taxa = c('Species_C','Species_H')) %>% treeNamer()

# Extract dates from 3 timetrees and get the mean node ages
tree_dates <- extractNodeAges(timeTrees,return_summary = 'mean')

# Read in signal from 5 multiple sequence alignments
mySignals <- getAlignmentSignal(rb_alignment_dir,species_info = timeTrees,suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))
mySupports <- getAlignmentSupport(mySignals,timeTrees,include_root = TRUE,dataset_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Plot percent breakdown of split support by alignment over time, looking for time effect at alpha = 0.05 (pre-Bonferroni correction)
# Ie. For each node in the tree, what percentage of support comes from each gene/subset?
timePlotter(node_age_df = tree_dates,tree_support = mySupports,lm_alpha = 0.05,return_stats = TRUE)
