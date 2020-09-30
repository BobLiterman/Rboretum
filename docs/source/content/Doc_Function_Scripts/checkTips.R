library(Rboretum)
sourceRboretum()

myTree_1 <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))
myTree_2 <- readRooted(rb_tree4_path,root_taxa = c('Species_C','Species_H'))

# Set some clades of interest
coi_1 <- c('Species_J','Species_M')
coi_2 <- c('Species_C','Species_H')
coi_3 <- c('Species_D','Species_K','Species_L')

# Visualize trees
coi_list <- list("One"=coi_1,"Two"=coi_2,"Three"=coi_3)
treePlotter(myTree_1,branch_weight = 1.5,xmax=10,to_color = coi_list,taxa_fontface = "bold")
treePlotter(myTree_2,branch_weight = 1.5, xmax=10,to_color = coi_list,taxa_fontface = "bold")

# Check if taxa from Clade of Interest 1 are in the trees
checkTips(tree=myTree_1,taxa = coi_1)
checkTips(tree=myTree_2,taxa = coi_1)

# Check if taxa from Clade of Interest 1 form a monophyletic group
checkTips(tree=myTree_1,taxa = coi_1,check_mono = TRUE)
checkTips(tree=myTree_2,taxa = coi_1,check_mono = TRUE)

# Check if taxa from Clade of Interest 1 form a monophyletic group at the root of the tree
checkTips(tree=myTree_1,taxa = coi_1,check_mono = TRUE, check_root=TRUE)
checkTips(tree=myTree_2,taxa = coi_1,check_mono = TRUE, check_root=TRUE)

# Check if taxa from Clade of Interest 2 form a monophyletic group at the root of the tree
checkTips(tree=myTree_1,taxa = coi_2,check_mono = TRUE, check_root=TRUE)
checkTips(tree=myTree_2,taxa = coi_2,check_mono = TRUE, check_root=TRUE)

# Check if taxa from Clade of Interest 3 form a monophyletic group
checkTips(tree=myTree_1,taxa = coi_3,check_mono = TRUE)
checkTips(tree=myTree_2,taxa = coi_3,check_mono = TRUE)
