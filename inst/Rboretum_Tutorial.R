# This is a brief walkthrough of some of the function of Rboretum 
# Example Data: 5 genes were sequenced for 20 species. Maxmimum likelihood trees were inferred from the alignments, and results are analyzed here.

library(Rboretum)

# In order to use alignment features, the user must source the scripts
source_python(system.file("", "Split_Processor.py", package = "Rboretum"))

# Rboretum can read trees in as rooted phylo objects if the outgroup can be specified
Gene1_file <- system.file("extdata", "Gene_1.nwk", package = "Rboretum")
geneTree1 <- readRooted(Gene1_file,root_taxa = c('Species_C','Species_H'))

# Extract monophyletic clades from trees
getTreeClades(geneTree1)

# Get tree spilts
getTreeSplits(geneTree1)

# We can quickly plot the tree using treePlotter
treePlotter(geneTree1)
treePlotter(geneTree1,xmax = 7) # Add buffer on right side of plot
treePlotter(geneTree1,branch_length = TRUE,xmax=3) # Plot tree with branch lengths
treePlotter(geneTree1,xmax = 7,node_label = "node") # Plot ggtree node IDs
treePlotter(geneTree1,xmax = 7,node_label = "node",node_label_box = TRUE) # Plot ggtree node IDs in a label box

# Highlight taxa
treePlotter(geneTree1,xmax = 7,to_color = c('Species_A','Species_B')) # Highlight taxa
treePlotter(geneTree1,xmax = 7,to_color = c('Species_A','Species_B'),colors = "blue") # Highlight taxa in  a different color

# Highlight groups of taxa
cladesOfInterest = list('Group 1' = c('Species_A','Species_F'),'Group 2'=c('Species_B','Species_O')) # Create a list of taxa to highlight
treePlotter(geneTree1,xmax = 7,to_color = cladesOfInterest) # Highlight taxa by group
treePlotter(geneTree1,xmax = 7,to_color = cladesOfInterest,highlight_legend = TRUE,color_branches = TRUE) # Highlight taxa by group

# If (1) there are multiple trees to read in, and (2) all trees share the same outgroup, they can be read in together in a few different ways
file_names <- c('Gene_1.nwk','Gene_2.nwk','Gene_3.nwk','Gene_4.nwk','Gene_5.nwk')
tree_paths <- paste(system.file("extdata",
                     file_names,
                     package = "Rboretum"), sep = ",")

allTrees <- readRooted(to_root = tree_paths,root_taxa = c('Species_C','Species_H')) # Read in multiple trees by passing multiple paths
data_dir <- system.file("extdata", package = "Rboretum")
allTrees <- readRooted(to_root = data_dir,suffix = 'nwk',root_taxa = c('Species_C','Species_H'),tree_names = c('Gene_1','Gene_2','Gene_3','Gene_4','Gene_5')) # Read in trees from <dir> that end with <suffix>, and assign new names

# Summarize a multiPhlyo
summarizeMultiPhylo(allTrees)
getTreeClades(allTrees,print_counts = TRUE)

# Get unique topologies from a multiPhylo
uniqueTrees <- getUniqueTopologies(allTrees,print_table = TRUE)

# Trim a tree to a set of taxa
trimmedTree <- treeTrimmer(geneTree1,taxa = c('Species_A','Species_B','Species_C'))
treePlotter(trimmedTree)

noSpeciesETree <- treeTrimmer(geneTree1,taxa = 'Species_E',remove = TRUE) # Remove taxa rather than retain
treePlotter(noSpeciesETree)

# Check if tips are in a tree
checkTips(trimmedTree,'Species_A')
checkTips(trimmedTree,'Species_F')

checkTips(allTrees,c('Species_A','Species_F')) # Check all trees in a multiPhylo
checkTips(allTrees,c('Species_A','Species_F'),check_mono = TRUE) # Also check if species are monophyletic
checkTips(allTrees,c('Species_A','Species_F'),check_mono = TRUE,check_root = TRUE) # Check if species are monophyletic and root
checkTips(allTrees,c('Species_C','Species_H'),check_mono = TRUE,check_root = TRUE) # Check if species are monophyletic and root

# Convert labels on tree
name_df <- read_tsv('Name_Conversion_Table.tsv')
renamed_tree <- convertLabels(geneTree1,name_df)
treePlotter(renamed_tree,xmax = 7,node_label_nudge = 0.2)

# Reading in mulitple sequence alignments to parse signal
alignmentSignal <- getAlignmentSignal(data_dir,suffix = 'phylip',use_gaps = FALSE,species_info = allTrees) # Read in alignments from <dir> that end with 'phylip', and process singal for all taxa in 'species_info'

# Get support for clades in allTrees
treeSupport <- getTreeSupport(alignmentSignal,allTrees,dataset_name = c('Gene_1','Gene_2','Gene_3','Gene_4','Gene_5'))
cladeSupport <- getTreeClades(allTrees,return_counts = TRUE)

# Plot support
treePlotter(geneTree1,tree_support = treeSupport,geom_size = c(1,3),xmax=8,node_label_nudge = 0.5,geom_alpha=0.5,node_label_fontface = 'bold',use_pies = TRUE,pie_legend_position = c(8,8,8,8))
treePlotter(allTrees,clade_support = cladeSupport,tree_support = treeSupport,geom_size = c(3,20),node_label = 'support',xmax=8,node_label_nudge = 0.5,geom_alpha=0.5,node_label_fontface = 'bold') %>% tandemPlotter()
treePlotter(allTrees,tree_support = treeSupport,geom_size = c(1,3),node_label = 'clade',xmax=8,node_label_nudge = 0.5,geom_alpha=0.5,node_label_fontface = 'bold',use_pies = TRUE,pie_legend_position = c(-1,-1,-1,-1)) %>% tandemPlotter()