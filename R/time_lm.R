#' Rboretum Time Linear Modeler
#'
#' This function takes information about node ages and node support, and uses linear models and Bonferroni correction to assess changes in relative phylogenetic utility over time among datasets
#' @param age_support_df 3-column dataframe including (1) Node Age, (2) Dataset Name, and (3) Percent Support
#' @param lm_alpha Run linear models at this alpha level [Default: 0.05]
#' @return 
#' @export
#' 
time_lm <- function(age_support_df,lm_alpha){
  
  # Check age_support_df
  if(!is.data.frame(age_support_df)){
    stop("'age_support_df' [usually passed from timePlotter()] does not appear to be a dataframe.")
  } else if(ncol(age_support_df)!=3){
    stop("'age_support_df' [usually passed from timePlotter()] does not contain 3 columns.")
  } else if(!all(names(age_support_df)==c('Node_Age','Dataset','Percent_Support'))){
    stop("'age_support_df' [usually passed from timePlotter()] does not contain columns named ['Node_Age','Dataset','Percent_Support'] .")
  }
  
  if(missing(lm_alpha)){
    lm_alpha <- 0.05
  }
  
  # Get dataset counts and set signficance based on Bonferroni
  datasets <- unique(age_support_df$Dataset)
  critical_value <- lm_alpha/length(datasets)
  
  stat_list <- c()
  r_list  <- c()
  intercept_list <- c()
  
  for(dataset in datasets){
    
    focal_lm <- age_support_df %>% filter(Dataset==dataset) %>% lm(Percent_Support ~ Node_Age, data=.)
    
    stat_list[[length(stat_list)+1]] <-t(data.frame(summary(focal_lm)$coefficients[2,]))
    r_list[[length(r_list)+1]] <-summary(focal_lm)$adj.r.squared
    intercept_list[[length(intercept_list)+1]] <-focal_lm$coefficients[[1]]
  }
  
  return(data.frame(matrix(unlist(stat_list), nrow=length(stat_list), byrow=T)) %>% 
           rename(Slope='X1',StdErr='X2',t_value='X3',p_value = 'X4') %>% 
           mutate(Dataset=datasets) %>%
           mutate(BF_Sig = ifelse(p_value <= critical_value,"Y","N")) %>%
           mutate(Adj_R_Sq=r_list) %>%
           mutate(Intercept=intercept_list))
}