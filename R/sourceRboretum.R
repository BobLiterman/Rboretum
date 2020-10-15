#' Rboretum Python Script Sourcer
#'
#' This function sources all Rborteum Python scripts for 'reticulate', and globally assigns the test data directories
#' @export
#' 
sourceRboretum <- function(){
  
  # Source python scripts
  source_python(system.file("", "Alignment_Species.py", package = "Rboretum"),envir=globalenv())
  source_python(system.file("", "Rboretum_pyRun.py", package = "Rboretum"),envir=globalenv())

  
  #source_python(system.file("", "Alignment_Splits.py", package = "Rboretum"),envir=globalenv())
  #source_python(system.file("", "Alignment_Patterns.py", package = "Rboretum"),envir=globalenv())
  #source_python(system.file("", "Alignment_Composition.py", package = "Rboretum"),envir=globalenv())
  #source_python(system.file("", "Species_Composition.py", package = "Rboretum"),envir=globalenv())
  #source_python(system.file("", "readAlignment.py", package = "Rboretum"),envir=globalenv())
  #source_python(system.file("", "siteCounter.py", package = "Rboretum"),envir=globalenv())
  
  rboretum_example_data_dir <<- system.file("extdata",package = 'Rboretum')
  rb_unroot_dir <<- paste0(rboretum_example_data_dir,'/unrootedTrees')
  rb_alignment_dir <<- paste0(rboretum_example_data_dir,'/alignments')
  rb_timeTree_dir <<- paste0(rboretum_example_data_dir,'/timeTrees')
  rb_name_file <<- paste0(rboretum_example_data_dir,'/Name_Conversion_Table.tsv')
  
  rb_tree1_path <<- paste0(rb_unroot_dir,'/Gene_1.nwk')
  rb_tree2_path <<- paste0(rb_unroot_dir,'/Gene_2.nwk')
  rb_tree3_path <<- paste0(rb_unroot_dir,'/Gene_3.nwk')
  rb_tree4_path <<- paste0(rb_unroot_dir,'/Gene_4.nwk')
  rb_tree5_path <<- paste0(rb_unroot_dir,'/Gene_5.nwk')
  
  rb_all_unrooted <<- c(rb_tree1_path,rb_tree2_path,rb_tree3_path,rb_tree4_path,rb_tree5_path)
  
  rb_align1_path <<- paste0(rb_alignment_dir,'/Gene_1.phylip')
  rb_align2_path <<- paste0(rb_alignment_dir,'/Gene_2.phylip')
  rb_align3_path <<- paste0(rb_alignment_dir,'/Gene_3.phylip')
  rb_align4_path <<- paste0(rb_alignment_dir,'/Gene_4.phylip')
  rb_align5_path <<- paste0(rb_alignment_dir,'/Gene_5.phylip')

  rb_all_align <<- c(rb_align1_path,rb_align2_path,rb_align3_path,rb_align4_path,rb_align5_path)
  
  rb_dummy_align_path <<- paste0(rb_alignment_dir,'/Gap_GC_N.fa')
  
  rb_timeTree1_path <<- paste0(rb_timeTree_dir,"/Chronogram_1.nwk")
  rb_timeTree2_path <<- paste0(rb_timeTree_dir,"/Chronogram_2.nwk")
  rb_timeTree3_path <<- paste0(rb_timeTree_dir,"/Chronogram_3.nwk")
  rb_timeTreeM_path <<- paste0(rb_timeTree_dir,"/Chronogram_MultiPhylo.nwk")
  
  rb_all_timeTree <<- c(rb_timeTree1_path,rb_timeTree2_path,rb_timeTree3_path)
}