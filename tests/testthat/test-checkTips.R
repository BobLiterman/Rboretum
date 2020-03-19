context("checkTips")

#read tree
tree2 <- readRooted("inst/extdata/Gene_1.nwk",root_taxa = c('Species_C','Species_H'))
trimmedTree <- treeTrimmer(tree2,taxa = paste0("Species_", LETTERS[2:14])) #base

test_that("taxa present as expected", {
  expect_true(checkTips(tree2,'Species_A'))
  expect_false(checkTips(tree2,'Species_Q'))
  expect_false(checkTips(trimmedTree,'Species_A'))
})

#read multiphylo
file_names <- c('inst/extdata/Gene_1.nwk','inst/extdata/Gene_2.nwk','inst/extdata/Gene_3.nwk')
allTrees <- readRooted(to_root = file_names,root_taxa = c('Species_C','Species_H'))
allTrees_difspp <- allTrees
allTrees_difspp[[1]] <- trimmedTree

test_that("taxa present as expected multiphylo", {
  expect_true(checkTips(allTrees,'Species_A'))
  expect_false(checkTips(allTrees,'Species_Q'))
  expect_false(checkTips(allTrees_difspp,'Species_A')) #one tree is missing the species
})
