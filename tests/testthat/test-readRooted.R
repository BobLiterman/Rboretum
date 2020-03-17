context("readRooted")

#read tree
tree <- readRooted("inst/extdata/Gene_1.nwk",root_taxa = c('Species_C','Species_H'))

test_that("reads phylo", {
  expect_s3_class(tree, "phylo")
})

test_that("test single root", {
  

})

test_that("multi root", {
  
  
  
})

test_that("warning when there is no data", {
  expect_warning(, "  ")
})

#read multiphylo
file_names <- c('inst/extdata/Gene_1.nwk','inst/extdata/Gene_2.nwk','inst/extdata/Gene_3.nwk')
allTrees <- readRooted(to_root = file_names,root_taxa = c('Species_C','Species_H'))

test_that("reads multiphylo", {
  expect_s3_class(allTrees, "multiPhylo")
})