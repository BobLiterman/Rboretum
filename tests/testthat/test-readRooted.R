context("readRooted")

#read tree
tree1 <- readRooted("inst/extdata/Gene_1.nwk",root_taxa = 'Species_H')
tree2 <- readRooted("inst/extdata/Gene_1.nwk",root_taxa = c('Species_C','Species_H'))

test_that("tree reads as phylo", {
  expect_s3_class(tree1, "phylo")
  expect_s3_class(tree2, "phylo")
})


test_that("tree is rooted", {
  expect_true(is.rooted(tree1))
  expect_true(is.rooted(tree2))
})

test_that("test rooting produces expected tree for 1 or 2 root species", {
  expect_equal(sort(tree2$tip.label), paste0("Species_", LETTERS[1:15])) #all species found

})

#read multiphylo
file_names <- c('inst/extdata/Gene_1.nwk','inst/extdata/Gene_2.nwk','inst/extdata/Gene_3.nwk')
allTrees <- readRooted(to_root = file_names,root_taxa = c('Species_C','Species_H'))

test_that("reads multiphylo", {
  expect_s3_class(allTrees, "multiPhylo")
})

test_that("trees are rooted", {
  expect_true(all(is.rooted(allTrees)))
})