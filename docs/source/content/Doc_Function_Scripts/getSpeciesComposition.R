library(Rboretum)
sourceRboretum()

# Set path to alignment data
myAlignmentFile <- rb_align1_path
myAlignmentDir <- rb_alignment_dir

# Get species composition information for a single alignment
getSpeciesComposition(alignment_path = myAlignmentFile)

# Get species composition information from all .phylip files in a directory, providing new names
getSpeciesComposition(alignment_path = myAlignmentDir,suffix = ".phylip",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))