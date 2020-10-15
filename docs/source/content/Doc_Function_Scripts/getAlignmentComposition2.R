devtools::install_github('BobLiterman/Rboretum',upgrade=FALSE)
library(Rboretum)

# Set path to alignment data
myManualAlignment <- 'inst/extdata/alignments/Gene_1.phylip'
mySpecies <- getAlignmentSpecies(myManualAlignment)
myName <- 'Gene_1'

rb_run_align_comp('inst/',myManualAlignment,mySpecies,myName)
rb_run_species_comp('inst/',myManualAlignment,mySpecies,myName)
