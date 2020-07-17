context("treeNamer")

# Read in multiPhylo
file_names <- c('inst/extdata/unrootedTrees/Gene_1.nwk','inst/extdata/unrootedTrees/Gene_2.nwk','inst/extdata/unrootedTrees/Gene_3.nwk')
three_trees <- readRooted(to_root = file_names,root_taxa = c('Species_C','Species_H'))
named_trees <- Rboretum::treeNamer(three_trees)

test_that("treeNamer returns a multiPhylo...", {
  expect_s3_class(named_trees, "multiPhylo")
})

test_that("Trees named correctly...", {
  expect_true(all(names(three_trees)==c('Gene_1.nwk','Gene_2.nwk','Gene_3.nwk')))
  expect_true(all(names(named_trees)==c('Tree_1','Tree_2','Tree_3')))
})