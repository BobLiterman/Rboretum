context("getTreeSupport")

#read phylo

#run

test_that('data types correct', {
  expect_is(testing_data,'data.frame')
  expect_is(testing_data$numbers, 'double') #fix
  expect_is(testing_data,'vector')  #given return_integer used
  expect_is(testing_data$numbers, 'double') #given return_integer used
})

test_that('right answers (singlton, bi, tri, quad, pent)  ', {

})


#repeat for multiphylo