context("getTreeSplits")

#read tree
tree2 <- readRooted("inst/extdata/Gene_1.nwk",root_taxa = c('Species_C','Species_H'))
#generate data frame
splits <- getTreeSplits(tree2)

#test ideas from https://katherinemwood.github.io/post/testthat/

test_that('data types correct', {
  expect_s3_class(splits,"data.frame") 
  expect_type(splits, "list") # A data frame is built from a list
  expect_s3_class(splits$Clade, 'factor')
  #expect_type(splits$Clade,"character")  #wondering if we want chr not factor (prior line)
  expect_s3_class(splits$Mirror_Clade, 'factor')
  #expect_type(splits$Mirror_Clade,"character")
  expect_type(splits$Split_Node, 'integer')
})

test_that('data dimensions correct', {
  expect_equal(ncol(splits), 3) 
  expect_equal(nrow(splits), 12) 
})

#read multiphylo
file_names <- c('inst/extdata/Gene_1.nwk','inst/extdata/Gene_2.nwk','inst/extdata/Gene_3.nwk',
                'inst/extdata/Gene_4.nwk','inst/extdata/Gene_5.nwk')
allTrees <- readRooted(to_root = file_names,root_taxa = c('Species_C','Species_H'))
splits <- getTreeSplits(allTrees)

test_that('data types correct', {
  expect_type(splits, "list")
  map(splits, ~expect_s3_class(.,"data.frame")) 
  map(splits, "Clade") %>% map(expect_s3_class, 'factor')
  map(splits, "Mirror_Clade") %>% map(expect_s3_class, 'factor')
  map(splits, "Split_Node") %>% map(expect_type, 'integer')
})

test_that('data dimensions correct', {
  expect_equal(length(splits),3)
  map(splits,ncol) %>% map(expect_equal,3) 
  map(splits,nrow) %>% map(expect_equal,12) 
})