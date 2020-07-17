context("getTreeSplits")

# Read in phylo
tree_1 <- readRooted("inst/extdata/unrootedTrees/Gene_1.nwk",root_taxa = c('Species_C','Species_H'))
tree_1_species <- tree_1$tip.label
tree_1_subtree <- subtrees(tree_1)
tree_1_subtree_clades <- purrr::map(.x=tree_1_subtree,.f=function(x){semiSorter(x$tip.label)}) %>% unlist()

tree_2 <- readRooted("inst/extdata/unrootedTrees/Gene_1.nwk",root_taxa = 'Species_B')
tree_2_species <- tree_2$tip.label
tree_2_subtree <- subtrees(tree_2)
tree_2_subtree_clades <- purrr::map(.x=tree_2_subtree,.f=function(x){semiSorter(x$tip.label)}) %>% unlist()

tree_1_splits <- getTreeSplits(tree_1)
tree_2_splits <- getTreeSplits(tree_2)

test_that('getTreeSplits returning correct data structures for phylo...', {
  expect_s3_class(tree_1_splits,"data.frame")
  expect_s3_class(tree_2_splits,"data.frame") 
  expect_type(tree_1_splits, "list")
  expect_type(tree_2_splits, "list")

  expect_type(tree_1_splits$Clade,"character")
  expect_type(tree_1_splits$Mirror_Clade,"character")
  expect_type(tree_1_splits$Split_Node,"integer")
  
  expect_type(tree_2_splits$Clade,"character")
  expect_type(tree_2_splits$Mirror_Clade,"character")
  expect_type(tree_2_splits$Split_Node,"integer")
})

test_that('getTreeSplits returning correct splits for phylo...', {
  expect_true(all(tree_1_splits$Clade %in% tree_1_subtree_clades))
  expect_true(all(tree_2_splits$Clade %in% tree_2_subtree_clades))
  expect_true(nrow(tree_1_splits) - nrow(tree_2_splits) = -1)
})

# Process as multiPhylo
trees <- c(tree_1,tree_2)
names(trees) <- c('Two_Root','One_Root')
tree_splits <- getTreeSplits(trees)

test_that('getTreeSplits returning correct data structures for multiPhylo...', {
  expect_s3_class(tree_splits,"list")
  expect_true(length(tree_splits) = length(trees))
  expect_s3_class(tree_splits[[1]],"data.frame")
  expect_s3_class(tree_splits[[2]],"data.frame")

  expect_type(tree_splits[[1]]$Clade,"character")
  expect_type(tree_splits[[1]]$Mirror_Clade,"character")
  expect_type(tree_splits[[1]]$Split_Node,"integer")
  
  expect_type(tree_splits[[2]]$Clade,"character")
  expect_type(tree_splits[[2]]$Mirror_Clade,"character")
  expect_type(tree_splits[[2]]$Split_Node,"integer")
})

test_that('getTreeSplits returning correct splits for multiPhylo...', {
  expect_true(all(tree_1_splits$Clade == tree_splits[[1]]$Clade))
  expect_true(all(tree_1_splits$Mirror_Clade == tree_splits[[1]]$Mirror_Clade))
  expect_true(all(tree_2_splits$Clade == tree_splits[[2]]$Clade))
  expect_true(all(tree_2_splits$Mirror_Clade == tree_splits[[2]]$Mirror_Clade))
})
