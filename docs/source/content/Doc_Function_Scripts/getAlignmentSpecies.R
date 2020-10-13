library(Rboretum)
sourceRboretum()

# Set alignment path
myAlignmentFile <- rb_align1_path

# Get sample IDs from alignment
getAlignmentSpecies(alignment_path = myAlignmentFile)
