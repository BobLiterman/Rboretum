library(Rboretum)
sourceRboretum()

# Set path to alignment data
myAlignmentFile <- rb_align1_path
myAlignmentDir <- rb_alignment_dir

# Get alignment composition information for a single alignment
getAlignmentComposition(alignment_path = myAlignmentFile)

# Get alignment composition information from all .phy files in a directory, providing new names
getAlignmentComposition(alignment_path = myAlignmentDir,suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Get alignment composition information from all .phy files in a directory, providing new names, considering only Species A - E
getAlignmentComposition(alignment_path = myAlignmentDir,species_info = 'Species_A;Species_B;Species_C;Species_D;Species_E',suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Get alignment composition from dummy alignment
getAlignmentComposition(alignment_path = rb_dummy_align_path)
