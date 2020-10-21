library(Rboretum)
sourceRboretum()

# Read in rooted tree
myTree <- readRooted(rb_tree1_path,'Species_C;Species_H')
getTreeRoot(myTree)

myTree <- readRooted(rb_tree1_path,'Species_B;Species_O')
getTreeRoot(myTree)

myTree <- readRooted(rb_tree1_path,'Species_F')
getTreeRoot(myTree)
