context("checkTips")

# Read in test tree
test_tree <- Rboretum::readRooted("inst/extdata/unrootedTrees/Gene_1.nwk",root_taxa = c('Species_C','Species_H'))

# Remove tip 'Species_A'
trimmed_tree <- Rboretum::treeTrimmer(test_tree,taxa = 'Species_A',remove=TRUE)

test_that("taxa present as expected", {
  expect_true(checkTips(test_tree,'Species_A'))
  expect_true(checkTips(test_tree,'Species_C',check_root=TRUE))
  expect_true(checkTips(test_tree,'Species_H',check_root=TRUE))
  expect_true(checkTips(test_tree,c('Species_C','Species_H'),check_root=TRUE))
  expect_false(checkTips(test_tree,'Species_Q'))
  expect_false(checkTips(trimmed_tree,'Species_A'))
})

# Read in  MultiPhylo
file_names <- c('inst/extdata/unrootedTrees/Gene_1.nwk','inst/extdata/unrootedTrees/Gene_2.nwk','inst/extdata/unrootedTrees/Gene_3.nwk')
three_trees <- readRooted(to_root = file_names,root_taxa = c('Species_C','Species_H'))

# Set tree 1 to trimmed tree without 'Species_A'
three_trees_b <- three_trees
three_trees_b[[1]] <- trimmed_tree

test_that("taxa present as expected multiphylo", {
  expect_true(checkTips(three_trees,'Species_A'))
  expect_false(checkTips(three_trees,'Species_Q'))
  expect_true(checkTips(three_trees,'Species_C',check_root=TRUE))
  expect_true(checkTips(three_trees,'Species_H',check_root=TRUE))
  expect_true(checkTips(three_trees,c('Species_C','Species_H'),check_root=TRUE))
  expect_false(checkTips(three_trees_b,'Species_A')) #one tree is missing the species
  expect_true(checkTips(three_trees_b,'Species_C',check_root=TRUE))
  expect_true(checkTips(three_trees_b,'Species_H',check_root=TRUE))
})
