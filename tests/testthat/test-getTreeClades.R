context("getTreeClades")

#read tree
tree2 <- readRooted("inst/extdata/Gene_1.nwk",root_taxa = c('Species_C','Species_H'))
clades <- getTreeClades(tree2)

test_that('data types correct', {
  expect_vector(clades)
  expect_type(clades, 'character')
})

test_that('result is correct', {

})

#repeat for return_counts (needs multiPhylo)
#read multiphylo
file_names <- c('inst/extdata/Gene_1.nwk','inst/extdata/Gene_2.nwk','inst/extdata/Gene_3.nwk')
allTrees <- readRooted(to_root = file_names,root_taxa = c('Species_C','Species_H'))
cladeSupport <- getTreeClades(allTrees,return_counts = TRUE)
#what should you get w only a phylo????

test_that('return counts returns df', {
  expect_type(cladeSupport,'list')  #why is it a list and not a df??
  expect_type(cladeSupport$Clade, 'character')
  expect_type(cladeSupport$Count, 'integer')
  expect_type(cladeSupport$Trees, 'character') #wondering if the Trees owuld be better as a vector of chr
})

#repeat for multiphylo
test_that('result is correct multiphylo', {
  
})