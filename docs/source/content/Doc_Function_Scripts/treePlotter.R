# remove.packages('Rboretum')
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
# devtools::install_github('BobLiterman/Rboretum',upgrade=FALSE)

library(Rboretum)
sourceRboretum()

# Read in a single tree
myTree <- readRooted(rb_tree1_path,'Species_C;Species_H')

# Plot tree
treePlotter(myTree)

# Adjust xmin and xmax to fix cut-off labels
treePlotter(myTree,xmin=-1,xmax=8)

# Increase the branch weight
treePlotter(myTree,xmin=-1,xmax=8,
            branch_weight=3)

# Adjust tip labels
treePlotter(myTree,xmin=-1,xmax=8,branch_weight=3,
            taxa_font_size=6,taxa_fontface = 'bold',taxa_offset = 0.2)

# Plot tree with branch lengths
treePlotter(myTree,
            branch_length=TRUE,xmax = 0.09)

# Adjust node labels
treePlotter(myTree,xmax=8,
            node_label = 'node',node_label_nudge = 0.1,node_label_box = FALSE,node_label_color = 'red')

# Highlighting specific taxa

myFavoriteSpecies <- c('Species_N','Species_G','Species_I')
myWorstSpecies <- c('Species_A','Species_B')

species_list <- list("Good" = myFavoriteSpecies, "Bad" = myWorstSpecies)

treePlotter(myTree,xmin=-1,xmax=8,branch_weight=3,taxa_font_size=6,taxa_fontface = 'bold',taxa_offset = 0.2,
            to_color = myFavoriteSpecies,colors = 'red',color_branches = TRUE)

treePlotter(myTree,xmin=-1,xmax=8,branch_weight=3,taxa_font_size=6,taxa_fontface = 'bold',taxa_offset = 0.2,
            to_color = myWorstSpecies,colors = 'blue')

treePlotter(myTree,xmin=-1,xmax=8,branch_weight=3,taxa_font_size=6,taxa_fontface = 'bold',taxa_offset = 0.2,
            to_color = species_list)

treePlotter(myTree,xmin=-1,xmax=8,branch_weight=3,taxa_font_size=6,taxa_fontface = 'bold',taxa_offset = 0.2,
            to_color = species_list,colors = c('green3','pink3'),highlight_legend = TRUE)

# Reversing the X-axis
plot1 <- treePlotter(myTree,xmin=-1,xmax=8,branch_weight=3,taxa_font_size=6,taxa_fontface = 'bold',taxa_offset = 0.2,to_color = myFavoriteSpecies,colors="red")

plot2 <- treePlotter(myTree,xmin=-1,xmax=8,branch_weight=3,taxa_font_size=6,taxa_fontface = 'bold',taxa_offset = 0.2,to_color = myWorstSpecies,colors="blue",
                     reverse_x=TRUE)

tandemPlotter(plot1,plot2)

# Read in multiPhylo of trees
myTrees <- readRooted(rb_unroot_dir,'Species_C;Species_H')

# Create a tree without Species A or B
myTrimmedTree <- treeTrimmer(myTree,myWorstSpecies,remove = TRUE)

# Add trimmed tree to mutliPhylo and rename trees
myNewTrees <- c(myTrees,myTrimmedTree) %>% treeNamer()

treePlotter(myTrees,xmin=-1,xmax=8,node_label_font_size = 3,taxa_offset = 0.2,taxa_fontface = "bold",taxa_font_size = 3,to_color=species_list)
treePlotter(myTrees,xmin=-1,xmax=8,node_label_font_size = 3,taxa_offset = 0.2,taxa_fontface = "bold",taxa_font_size = 3,to_color=species_list,basic_plot = TRUE)

treePlotter(myNewTrees,xmin=-1,xmax=8,node_label_font_size = 3,taxa_offset = 0.2,taxa_fontface = "bold",taxa_font_size = 3,to_color=species_list)
treePlotter(myNewTrees,xmin=-1,xmax=8,node_label_font_size = 3,taxa_offset = 0.2,taxa_fontface = "bold",taxa_font_size = 3,to_color=species_list,basic_plot = TRUE)

# Visualize node support
tree_clades <- getTreeClades(myTrees,return_counts = TRUE)
treePlotter(myTrees,xmin=-1,xmax=8,taxa_offset = 0.2,taxa_fontface = "bold",to_color=species_list,
            clade_support = tree_clades,node_label = 'none',geom_alpha = 1,geom_size = 10,legend_shape_size = 10,legend_font_size = 15)

# Plot support from a multiple sequence alignment onto the phylogeny

# Get signal from an alignment file, as it relates to species from myTree
mySignal <- getAlignmentSignal(rb_align1_path,species_info = myTree,alignment_name = "Gene A")

# Break down support from mySignal as it relates to splits from myTree
mySupport <- getAlignmentSupport(mySignal,myTree,include_root = TRUE,dataset_name = "Gene A")

# Plot tree with node labels corresponding to the number of sites supporting each split
treePlotter(myTree,xmin=-1,xmax=8,branch_weight=3,taxa_font_size=6,taxa_fontface = 'bold',taxa_offset = 0.2,
            tree_support = mySupport,node_label='support',plot_root_support = TRUE)

# Plot tree with node geoms, with sizes scaled to the number of sites supporting each split (rescaled from 5 - 30)
treePlotter(myTree,xmin=-1,xmax=8,branch_weight=3,taxa_font_size=6,taxa_fontface = 'bold',taxa_offset = 0.2,
            tree_support = mySupport,node_label='support',node_label_nudge = 0.1,plot_root_support = TRUE,geom_size = c(5,30),geom_alpha = 1)

# Plot support from 2+ multiple sequence alignments onto the phylogeny

# Get signal from an alignment files, as it relates to species from myTrees
mySignals <- getAlignmentSignal(rb_alignment_dir,species_info = myTrees,suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Break down support from mySignals as it relates to splits from myTrees
mySupports <- getAlignmentSupport(mySignals,myTrees,include_root = TRUE,dataset_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Plot support values as geoms scaled to the total support summed over all alignments
treePlotter(myTrees,xmin=-1,xmax=8,branch_weight=3,taxa_font_size=4,taxa_fontface = 'bold',taxa_offset = 0.1,
            tree_support = mySupports,node_label='support',node_label_nudge = 0.1,plot_root_support = TRUE,geom_size = c(5,30),geom_alpha = 1)

# Plot support values as pies scaled to the total support summed over all alignments, and color coded by alignment
# Note when using pies that the position of the legend needs to be set manually
treePlotter(myTrees,xmin=-1,xmax=8,branch_weight=3,taxa_font_size=4,taxa_fontface = 'bold',taxa_offset = 0.1,
            tree_support = mySupports,use_pies = TRUE,node_label='support',node_label_fontface = "bold",geom_size = 4,node_label_color = "red",node_label_font_size = 3,node_label_nudge = 0.4,plot_root_support = TRUE,geom_alpha = 1,pie_scaler=0.45,pie_colors = c('red','blue','green','orange','purple'),pie_legend_position = c(1,1,14,14),scale_range = c(75,250))
