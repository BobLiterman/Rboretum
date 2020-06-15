context("getUniqueTopologies")

#read multiphylo

#run

test_that("test multiphylo", {
  
  expect_s3_class(tree, "multiPhylo")
  
})

#return table option
test_that('data types correct', {
  expect_is(testing_data,'table')
})



