library(Rboretum)

myVector <- c('a','b','b','c','c','c')
myTable <- table(myVector)

myTable[['a']]
tableCount(search_table = myTable,name = 'a')

myTable[['b']]
tableCount(myTable,'b')

myTable[['c']]
tableCount(myTable,'c')

myTable[['d']]
tableCount(myTable,'d')
