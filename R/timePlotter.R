#' Rboretum Time-Dependent Signal Plotter
#'
#' This function takes information about node ages and node support, and plots relative changes in node support over time among datasets. 
#' @param node_age_df Output from extractNodeAges(phylo) or extractNodeAges(multiPhylo,return_summary = 'mean' or 'median') [Or a two-column dataframe containing (1) semicolon separated clades and (2) estimated node ages]
#' @param tree_support Output from getAlignmentSupport with information about all clades in node_age_df and 2+ alignment support columns
#' @param plot_datasets OPTIONAL: Character vector specifying the datasets desired for the plot (Must be 2+ and present in 'tree_support')
#' @param all_sites_col OPTIONAL (Very special case): If the alignments from tree_support contain overlapping sites AND there is is a column in the data corresponding to all non-overlapping sites, specify the column name for all nonoverlapping sites to be used when calculating proportions
#' @param lm_alpha OPTIONAL: Run linear models at this alpha level (must be > 0 and < 1; Alpha will be automatically Bonferroni corrected based on dataset count) [Default: Do not run linear models]
#' @param wrap OPTIONAL (TRUE/FALSE); Plot datasets in separate facets [Default: FALSE; do not wrap and produce 1 plot]
#' @param wrap_scales OPTIONAL ('free'/'fixed'): If wrapping, allow free or fixed scales on each facet [Default: 'fixed']
#' @param return_stats OPTIONAL (TRUE/FALSE): In addition to returning the plot, also return the statistical output from linear models in a list
#' @return Plot or list (depending on return_stats)
#' @export
#' 
timePlotter <- function(node_age_df,tree_support,plot_datasets,all_sites_col,lm_alpha,wrap,wrap_scales,return_stats){
  
  # Check node age df (Should have 2 columns with prescribed names)
  if(missing(node_age_df)){
    stop("'node_age_df' [output from extractNodeAges() or extractNodeAges(return_summary=mean/median)] is required")
  } else if(!is.data.frame(node_age_df)){
    stop("'node_age_df' should be the output from extractNodeAges() or extractNodeAges(return_summary=mean/median)")
  } 
  
  # Remove variation column
  if(ncol(node_age_df)==3){
    node_age_df <- node_age_df %>% select(-3)
  }
  
  if(!ncol(node_age_df)==2){
    stop("'node_age_df' should be the output from extractNodeAges() or extractNodeAges(return_summary=mean/median)")
  } else{
    
    # Check column names
    if(names(node_age_df)[1] == 'Clade' && names(node_age_df)[2] %in% c('Node_Age','Mean_Node_Age','Median_Node_Age')){
      node_age_df <- node_age_df %>% `names<-`(c('Clade','Node_Age'))
    } else{
      stop("'node_age_df' should be the output from extractNodeAges() or extractNodeAges(return_summary=mean/median)")
    }
  }

  # Ensure tree_support contains at least two columns of data
  if(missing(tree_support)){
    stop("'tree_support' (with 2+ data columns) is required.")
  } else if(!ncol(tree_support) > 2){
    stop("'tree_support' needs to contain 2+ columns of support data for timePlotter().")    
  } else{
    
    # Ensure node_age_df contains information about all splits in tree_support
    node_age_clades <- node_age_df$Clade
    
    if(!Rboretum::isAlignmentSupport(test_object = tree_support,test_clade = node_age_clades,partial = TRUE)){
      stop("Clades from 'node_age_df' are missing from 'tree_support")
    } else{
      
      tree_support <- tree_support %>% filter(Clade %in% node_age_clades)
      tree_support_clades <- tree_support %>% pull(Clade)
      
      # Check for an 'All Sites' column
      if(missing(all_sites_col)){
        tree_support <- tree_support
        all_sites_filter <- FALSE
      } else if(!is.character(all_sites_col)){
        stop("'all_sites_col' should be a character identifying a column in the 'tree_support' dataframe")
      } else if(length(all_sites_col)!=1){
        stop("'all_sites_col' should be a character identifying a column in the 'tree_support' dataframe")
      } else if(!all_sites_col %in% names(tree_support)){
        stop("'all_sites_col' should be a character identifying a column in the 'tree_support' dataframe")
      } else{
        all_sites_support <- tree_support %>% pull(all_sites_col)
        tree_support <- tree_support %>% select(-all_sites_col)
        all_sites_filter <- TRUE
      }
      
      support_col_nums <- 2:ncol(tree_support)
      support_col_names <- names(tree_support)[support_col_nums]
      support_count <- length(support_col_nums)
      dataset_colors <- viridisLite::viridis(support_count)
      
      node_age_vector <- node_age_df$Node_Age
      names(node_age_vector) <- node_age_df$Clade
      tree_support_ages <- purrr::map(.x=tree_support_clades,.f=function(x){node_age_vector[x]}) %>% unlist()
    }
  }
  
  # Get wrap/scale information
  if(missing(wrap)){
    wrap <- FALSE
    wrap_scales <- NULL
  } else if(!is.logical(wrap)){
    warning("'wrap' should be TRUE or FALSE. Plotting without facet wrapping...")
    wrap <- FALSE
  } else if(wrap){
    if(missing(wrap_scales)){
      wrap_scales <- 'free'
    } else if(!is.character(wrap_scales)){
      wrap_scales <- 'free'
      warning("'wrap_scales' should be 'free' [default] or 'fixed'. Plotting with free scales...")
    } else if(length(wrap_scales) > 1){
      wrap_scales <- 'free'
      warning("'wrap_scales' should be 'free' [default] or 'fixed'. Plotting with free scales...")
    } else if(!(wrap_scales %in% c('free','fixed'))){
      wrap_scales <- 'free'
      warning("'wrap_scales' should be 'free' [default] or 'fixed'. Plotting with free scales...")
    }
  }
  
  plot_df <- tibble(Clade=character(),Node_Age=numeric(),Dataset=character(),Total_Support=numeric(),Percent_Support=numeric(),Plot_Color=character())
  
  if(!all_sites_filter){
    # 'Percent_Support' = Dataset support / Sum of Node support across datasets
    
    # Get total support for each clade across datasets
    total_clade_support <- rowSums(tree_support[support_col_nums])
    
    for(i in 1:support_count){
      dataset_name <- support_col_names[i]
      dataset_counts <- pull(tree_support,support_col_nums[i])
      dataset_percents = 100*(dataset_counts/total_clade_support)
      temp_df <- tibble(Clade=tree_support_clades,Node_Age=tree_support_ages,Dataset=dataset_name,Total_Support=total_clade_support,Percent_Support=dataset_percents,Plot_Color=dataset_colors[[i]])
      plot_df <- rbind(plot_df,temp_df)
    }
  } else{
    # 'Percent_Support' = Dataset support / Support from 'all sites' 
    
    for(i in 1:support_count){
      dataset_name <- support_col_names[i]
      dataset_counts <- pull(tree_support,support_col_nums[i])
      dataset_percents = 100*(dataset_counts/all_sites_support)
      temp_df <- tibble(Clade=tree_support_clades,Node_Age=tree_support_ages,Dataset=dataset_name,Total_Support=all_sites_support,Percent_Support=dataset_percents,Plot_Color=dataset_colors[[i]])
      plot_df <- rbind(plot_df,temp_df)
    }
  }
  
  # Run linear models if 'lm_alpha' arguement is passed
  if(missing(lm_alpha)){
    run_lm <- FALSE
  } else if(!is.numeric(lm_alpha)){
    warning("'lm_alpha' should be a critical p-value threshold...performing statstical assessment at default alpha of 0.05...")
    lm_alpha <- 0.05
    run_lm <- TRUE
  } else if(lm_alpha >= 1 | lm_alpha <= 0){
    warning("'lm_alpha' should be a critical p-value threshold between 0 and 1...performing statstical assessment at default alpha of 0.05...")
    lm_alpha <- 0.05
    run_lm <- TRUE
  } else{
    run_lm <- TRUE
  }
  
  if(run_lm){
    lm_df <- time_lm(plot_df %>% select(Node_Age,Dataset,Percent_Support),lm_alpha)
    
    plot_df <- plot_df %>% 
      left_join(lm_df,by='Dataset') %>%
      mutate(Slope = 100*Slope) %>%
      mutate(Plot_Fill = ifelse(BF_Sig =="N","white",Plot_Color)) %>%
      mutate(Plot_Shape = ifelse(BF_Sig =="N",1,23))
  } else{
    plot_df <- plot_df %>% mutate(Plot_Fill = Plot_Color) %>%
      mutate(Plot_Shape = 21)
  }
  
  # Pull out color information for legend
  color_df <- plot_df[!duplicated(plot_df$Dataset),] %>% select(Dataset,Plot_Color)
  plot_colors <- pull(color_df,Plot_Color)
  names(plot_colors) <- pull(color_df,Dataset)
  
  # Pull out desired datasets if requested
  if(missing(plot_datasets)){
    plot_df <- plot_df
  } else if(!is.character(plot_datasets)){
    warning("'plot_datasets' should be a character vector of desired datasets to plot. Plotting all datasets...")
  } else if(length(plot_datasets)<2){
    warning("'plot_datasets' should be a character vector of two or more desired datasets to plot. Plotting all datasets...")
  } else if(!all(plot_datasets %in% support_col_names)){
    warning("Column names from 'plot_datasets' do not all occur in 'tree_support'. Plotting all datasets...")
  } else{
    plot_df <- plot_df %>% filter(Dataset %in% plot_datasets)
  }
  
  base_plot <- ggplot(plot_df,aes(x=Node_Age,y=Percent_Support,color=Dataset)) +
    geom_jitter(size=4,aes(fill=Plot_Fill,shape=Plot_Shape)) +
    geom_smooth(size=1.5,method=lm,se=FALSE,show.legend = FALSE,aes(fill=Plot_Fill)) +
    scale_fill_identity() +
    scale_color_manual(values = plot_colors,guide="legend") +
    scale_shape_identity() +
    theme_bw() +
    xlab("\nTime Since Split\n") +
    ylab("Proportion of Phylogenetic Signal (%)\n") +
    scale_x_reverse() +
    theme(axis.title=element_text(size=14,face="bold"),
          axis.text.x=element_text(size=13),
          axis.text.y = element_text(size=13,face="bold"))
  
  if(!wrap){ # One plot
    return_plot <- base_plot +
      guides(shape = guide_legend(override.aes = list(size = 4))) +
      theme(legend.title = element_text(hjust = 0.5, size=14, face="bold"),
            legend.text = element_text(size=12))
    
  } else{ # Wrapped plot
    return_plot <- base_plot +
      facet_wrap(~Dataset,scales=wrap_scales) +
      theme(panel.spacing = unit(1.5, "lines"),
            legend.position = "none",
            strip.text.x = element_text(size = 12,face="bold")) 
  }
  
  if(missing(return_stats)){
    return(return_plot)
  } else if(!is.logical(return_stats)){
    warning("'return_stats' should be TRUE or FALSE. Not returning stats table...")
    return(return_plot)
  } else if(run_lm & return_stats){
    return(list(Plot=return_plot,Stats=lm_df))
  } else{
    return(return_plot)
  }
}