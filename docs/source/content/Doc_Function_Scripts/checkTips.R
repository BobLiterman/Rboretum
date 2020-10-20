library(Rboretum)
sourceRboretum()

# Read in two test trees with different topologies
myTree_1 <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))
myTree_2 <- readRooted(rb_tree4_path,root_taxa = c('Species_C','Species_H'))

# Visualize trees

# Set some clades of interest
coi_1 <- c('Species_J','Species_M')
coi_2 <- c('Species_C','Species_H')
coi_3 <- c('Species_D','Species_K','Species_L')
coi_list <- list("One"=coi_1,"Two"=coi_2,"Three"=coi_3)

# Make plots
plot_1 <- treePlotter(myTree_1,branch_weight = 1.5,xmax=8,to_color = coi_list,taxa_fontface = "bold",taxa_offset = 0.2,xmin=-0.2)
plot_2 <- treePlotter(myTree_2,branch_weight = 1.5, xmax=8,to_color = coi_list,taxa_fontface = "bold",taxa_offset = 0.2,reverse_x=TRUE,xmin=-0.2)
tandemPlotter(plot_1,plot_2)

# Check if taxa from Clade of Interest 1 (Species_J/Species_M) are in the trees
checkTips(tree=myTree_1,taxa = coi_1)
checkTips(tree=myTree_2,taxa = coi_1)

# Chec if nonsense taxa are in the trees
checkTips(tree=myTree_1,taxa = "Not_A_Species")
checkTips(tree=myTree_2,taxa = "Also_Not_A_Species")

# Check if taxa from Clade of Interest 1 (Species_J/Species_M) form a monophyletic group
checkTips(tree=myTree_1,taxa = coi_1,check_mono = TRUE)
checkTips(tree=myTree_2,taxa = coi_1,check_mono = TRUE)

# Check if taxa from Clade of Interest 1 (Species_J/Species_M) form a monophyletic group at the root of the tree
checkTips(tree=myTree_1,taxa = coi_1,check_mono = TRUE, check_root=TRUE)
checkTips(tree=myTree_2,taxa = coi_1,check_mono = TRUE, check_root=TRUE)

# Check if taxa from Clade of Interest 2 (Species_C/Species_H) form a monophyletic group at the root of the tree
checkTips(tree=myTree_1,taxa = coi_2,check_mono = TRUE, check_root=TRUE)
checkTips(tree=myTree_2,taxa = coi_2,check_mono = TRUE, check_root=TRUE)

# Check if taxa from Clade of Interest 3 (Species_D/Species_K/Species_L)form a monophyletic group
checkTips(tree=myTree_1,taxa = coi_3,check_mono = TRUE)
checkTips(tree=myTree_2,taxa = coi_3,check_mono = TRUE)
