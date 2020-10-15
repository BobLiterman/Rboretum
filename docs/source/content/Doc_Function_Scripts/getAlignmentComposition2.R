devtools::install_github('BobLiterman/Rboretum',upgrade=FALSE)
library(Rboretum)
sourceRboretum()

# Set path to alignment data
myManualAlignment <- 'inst/extdata/alignments/Gene_1.phylip'
mySpecies <- getAlignmentSpecies(myManualAlignment)
myName <- 'Gene_1'

getAlignmentComposition(myManualAlignment,mySpecies,myName)
getSpeciesComposition(myManualAlignment,mySpecies,myName)
getAlignmentSignal(myManualAlignment,mySpecies,TRUE,myName)
getAlignmentPatterns(myManualAlignment,mySpecies,TRUE,myName)