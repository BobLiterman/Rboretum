library(Rboretum)

# Set test data directory
sourceRboretum()

semicolon_single_char <- 'a;b;c;d;e;f'
semiVector(semicolon_single_char)

semicolon_multi_char <- c('a;b;c','d;e;f')
semiVector(semicolon_multi_char)
