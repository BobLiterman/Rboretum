#' Rboretum Python Script Sourcer
#'
#' This function sources all Rborteum Python scripts for 'reticulate', and globally assigns the test data directories
#' @export
#' 
sourceRboretum <- function(){
  
  # Set python and local package directory
  user_python <<- py_config()$python
  rb_lib_dir <<- system.file(package = "Rboretum")
  
  # Source python scripts
  source_python(system.file("", "Alignment_Species.py", package = "Rboretum"),envir=globalenv())
  source_python(system.file("", "Rboretum_pyRun.py", package = "Rboretum"),envir=globalenv())
  
  # Set data directory
  rboretum_example_data_dir <<- system.file("extdata",package = 'Rboretum')
  
  # Name Conversion Table
  rb_name_file <<- paste0(rboretum_example_data_dir,'/Name_Conversion_Table.tsv')
  
  # Trees
  rb_unroot_dir <<- paste0(rboretum_example_data_dir,'/unrootedTrees')
  
  rb_tree1_path <<- paste0(rb_unroot_dir,'/Gene_1.nwk')
  rb_tree2_path <<- paste0(rb_unroot_dir,'/Gene_2.nwk')
  rb_tree3_path <<- paste0(rb_unroot_dir,'/Gene_3.nwk')
  rb_tree4_path <<- paste0(rb_unroot_dir,'/Gene_4.nwk')
  rb_tree5_path <<- paste0(rb_unroot_dir,'/Gene_5.nwk')
  
  rb_all_unrooted <<- c(rb_tree1_path,rb_tree2_path,rb_tree3_path,rb_tree4_path,rb_tree5_path)
  
  # Alignments
  rb_alignment_dir <<- paste0(rboretum_example_data_dir,'/alignments')
  
  rb_align1_path <<- paste0(rb_alignment_dir,'/Gene_1.phy')
  rb_align2_path <<- paste0(rb_alignment_dir,'/Gene_2.phy')
  rb_align3_path <<- paste0(rb_alignment_dir,'/Gene_3.phy')
  rb_align4_path <<- paste0(rb_alignment_dir,'/Gene_4.phy')
  rb_align5_path <<- paste0(rb_alignment_dir,'/Gene_5.phy')

  rb_all_align <<- c(rb_align1_path,rb_align2_path,rb_align3_path,rb_align4_path,rb_align5_path)
  
  rb_dummy_align_path <<- paste0(rb_alignment_dir,'/Dummy_Alignment.fa')
  
  # TimeTrees
  rb_timeTree_dir <<- paste0(rboretum_example_data_dir,'/timeTrees')
  
  rb_timeTree1_path <<- paste0(rb_timeTree_dir,"/Gene_1_TimeTree.nwk")
  rb_timeTree2_path <<- paste0(rb_timeTree_dir,"/Gene_2_TimeTree.nwk")
  rb_timeTree3_path <<- paste0(rb_timeTree_dir,"/Gene_3_TimeTree.nwk")
  rb_timeTree4_path <<- paste0(rb_timeTree_dir,"/Gene_4_TimeTree.nwk")
  rb_timeTree5_path <<- paste0(rb_timeTree_dir,"/Gene_5_TimeTree.nwk")
}