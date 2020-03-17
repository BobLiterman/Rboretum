context("extractNodeAges")

#read phylo

#run

test_that('data types correct', {
  expect_is(testing_data,'data.frame')
  expect_is(testing_data$numbers, 'double') #fix
})



#repeat for multiphylo

test_that('node summary summarizes', {
  expect_is(testing_data,'data.frame')
  expect_is(testing_data$numbers, 'double') #fix
})