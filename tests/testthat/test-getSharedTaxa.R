context("getSharedTaxa")

#read multiphylo

#run

test_that('data types correct', {
  expect_is(testing_data,'vector')
  expect_is(testing_data, 'character')
})

test_that("return correct subset", {
  
  expect_equal(testingdata, paste0("t", 4:6)) #include expected tips
  
})