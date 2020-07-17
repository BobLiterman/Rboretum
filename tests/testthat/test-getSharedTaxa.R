context("getSharedTaxa")

# Generate phylos
phylo_1 <- ape::rtree(10)
phylo_2 <- ape::rtree(20)
phylo_3 <- ape::rtree(30)

ten_taxa <- phylo_1$tip.label %>% naturalsort()
twenty_taxa <- phylo_2$tip.label %>% naturalsort()
thirty_taxa <- phylo_3$tip.label %>% naturalsort()

multi_12 <- c(phylo_1,phylo_2)
multi_13 <- c(phylo_1,phylo_3)
multi_23 <- c(phylo_2,phylo_3)
multi_33 <- c(phylo_3,phylo_3)

shared_12 <- getSharedTaxa(multi_12)
shared_13 <- getSharedTaxa(multi_13)
shared_23 <- getSharedTaxa(multi_23)
shared_33 <- getSharedTaxa(multi_33)

test_that('getSharedTaxa retuns a character vector...', {
  expect_is(shared_12,'vector')
  expect_type(shared_12, 'character')
  expect_is(shared_13,'vector')
  expect_type(shared_13, 'character')
  expect_is(shared_23,'vector')
  expect_type(shared_23, 'character')
  expect_is(shared_33,'vector')
  expect_type(shared_33, 'character')
})

test_that("getSharedTaxa retuns the correct taxa...", {
  expect_true(all(shared_12 == ten_taxa))
  expect_true(all(ten_taxa == shared_12))
  
  expect_true(all(shared_13 == ten_taxa))
  expect_true(all(ten_taxa == shared_13))
  
  expect_true(all(shared_23 == twenty_taxa))
  expect_true(all(twenty_taxa == shared_23))
  
  expect_true(all(shared_33 == thirty_taxa))
  expect_true(all(thirty_taxa == shared_33))
})