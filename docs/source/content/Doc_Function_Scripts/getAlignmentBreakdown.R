library(Rboretum)
sourceRboretum()

# Get alignment Breakdown information for a single alignment
getAlignmentBreakdown(alignment_path = rb_align1_path)

# Get alignment Breakdown information from all .phy files in a directory, providing new names, consider gaps as missing data
getAlignmentBreakdown(alignment_path = rb_alignment_dir,use_gaps = FALSE,suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Get alignment Breakdown information from all .phy files in a directory, providing new names, consider gaps as indel data, and only analyze species A-E
getAlignmentBreakdown(alignment_path = rb_alignment_dir,species_info = c('Species_A','Species_B','Species_C','Species_D','Species_E'),use_gaps = TRUE,suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Get alignment Breakdown from dummy alignment, with and without gap support
getAlignmentBreakdown(alignment_path = rb_dummy_align_path)

getAlignmentBreakdown(alignment_path = rb_dummy_align_path,use_gaps = FALSE)
