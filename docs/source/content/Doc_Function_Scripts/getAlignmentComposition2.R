devtools::install_github('BobLiterman/Rboretum',upgrade=FALSE)
library(Rboretum)
sourceRboretum()

# Set path to alignment data
myAlignmentFile <- rb_align1_path
myAlignmentDir <- rb_alignment_dir

# Get alignment composition information for a single alignment
getAlignmentComposition(alignment_path = myAlignmentFile)
getSpeciesComposition(alignment_path = myAlignmentFile)

# Get alignment composition information from all .phylip files in a directory, providing new names
getAlignmentComposition(alignment_path = myAlignmentDir,suffix = ".phylip",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Get alignment composition from dummy alignment
getAlignmentComposition(alignment_path = rb_dummy_align_path)

