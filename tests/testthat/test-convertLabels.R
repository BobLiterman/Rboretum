context("convertLabels")

#read phylo

#run

test_that("test multiphylo", {
  
  expect_s3_class(tree, "phylo")
  
})

test_that("names substituted correctly", {
  

  
})


#repeat for multiphylo

test_that("test multiphylo", {
  
  expect_s3_class(tree, "multiPhylo")
  
})