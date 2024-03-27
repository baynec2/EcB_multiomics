## Function to assess the metadata for potential statistically significant confounding variables. ##
metadata_assesment = function(data,
                              sample_id = sample_id,
                              metric = row_id,
                              features,
                              outcome = norm_value,
                              outcome_c = "norm_value"){
  
  
  # Splitting the inital data frame into the data columns and the 
  input = data %>% 
    dplyr::select({{sample_id}},{{metric}},{{outcome}})
  
  md_in_all = data %>% 
    dplyr::select({{sample_id}},{{metric}})
  
  md = data %>% 
    dplyr::select(-{{sample_id}},-{{metric}},-{{outcome}})
  
  # What columns are categorical
  categorical_cols = purrr::map_lgl(md,~!is.numeric(.))
  
  categorical_md = md %>% 
    dplyr::select(which(categorical_cols))
  
  # number of distinct things in each categorical category
  col_numbers_in_cat = categorical_md %>% 
    purrr::map_dbl(~length(unique(na.omit(.))))
  
  ##Of the categories that are categorical, which have 2 categories
  # Sample ID should be in the first column. Adding it here. 
  cols_2 = categorical_md[,which(col_numbers_in_cat == 2)]
  
  ##Of the categories that are categorical, which have more than 2 categories. 
  cols_g_2_cat = categorical_md[,which(col_numbers_in_cat > 2)]
  
  # What columns are continuous
  continuous_cols = purrr::map_lgl(md,is.numeric)
  
  continuous_md = md[,which(continuous_cols)]
  
 
  input_md_continuous = dplyr::bind_cols(input,continuous_md)
  input_md_2_col = dplyr::bind_cols(input,cols_2)
  input_md_g_2_col = dplyr::bind_cols(input,cols_g_2_cat)
 
  
  # Lets do the continuous stats first # 
  stats = data.frame()
 
  # Outer loop with the different variables to test 
for(x in names(which(continuous_cols))){
  
  for(i in features){
    loop = dplyr::filter(input_md_continuous,{{metric}} == i)
    #Let's look at the continuous stuff first 
    c_mod = loop %>% 
      rstatix::cor_test({{outcome}},x,method = "pearson",use="complete.obs") %>% 
      dplyr::select(variable = var2,p) %>% 
      dplyr::mutate({{metric}} := i)
    
    stats = bind_rows(stats,c_mod)
  }
}
  
  # Outer loop 
  for(x in names(cols_2)){
    # Categorical two categories
    for(i in features){
      loop = dplyr::filter(input_md_2_col,{{metric}} == i)
      
      #Let's do the stats for the things with only two values first. 
      inner_loop = rstatix::wilcox_test(data = loop,as.formula(paste0(outcome_c," ~ ",x))) %>% 
       dplyr::mutate(variable = x,
               {{metric}} := i) %>% 
        select(variable,p,{{metric}})
    
      stats = bind_rows(stats,inner_loop)
      
    }
    
}
  
  # Outer loop 
  for(x in names(cols_g_2_cat)){
    # Categorical two categories
    for(i in features){
      loop = dplyr::filter(input_md_g_2_col,{{metric}} == i)
      
      #Let's do the stats for the things with more than two values now. 
      inner_loop = rstatix::kruskal_test(data = loop,as.formula(paste0(outcome_c," ~ ",x))) %>% 
        dplyr::mutate(variable = x,{{metric}} := i) %>% 
        dplyr::select(variable,p,{{metric}})
      
      
      stats = dplyr::bind_rows(stats,inner_loop)
      
    }
    
  }
  return(stats)
}


## Function for making ROC plots using proteomics data ##

make_roc_proteomics = function(data,topn,ncol = 5){
  library(gridExtra)
  
  list = list()
  
  for(i in topn){
    
    loop = filter(data, Gene == i)
    
    mod = glm(as.factor(Condition) ~ Abundance, 
              data = data,
              family = "binomial") 
    
    predicted <- predict(mod ,family = "binomial", loop, type="response")
    
    #define object to plot
    rocobj <- pROC::roc(loop$Condition, predicted)
    
    auc = round(as.numeric(gsub(".*: ", "",rocobj$auc)),2)
    #create ROC plot
    p1 = pROC::ggroc(rocobj)+
      ggtitle(paste0(i," AUC: ",auc))
    
    
    list[[i]] = p1
  }
  plot = do.call("grid.arrange", c(list, ncol = ncol))
  
  return(plot)
}



make_violin_plot_proteomics = function(data,topn,ncol=5){
  library(gridExtra)
  # initializing list
  list = list()
  #looping through each biomarker
  for(i in topn){
    
    loop = filter(data, Gene == i)
    stats = rstatix::t_test(loop, Abundance ~ Condition) %>% 
    rstatix::add_y_position(fun = "median")
  
  p1 = loop %>% 
    ggplot(aes(Condition,Abundance))+
    geom_violin()+
    geom_point()+
    theme(axis.text.x = element_text(angle = 90))+
    ggpubr::stat_pvalue_manual(stats, label = "p.adj.signif", tip.length = 0.001)+
    ggtitle(i)+
    theme(axis.title.x = element_blank())
  
  list[[i]] = p1
  }
  
  plot = do.call("grid.arrange", c(list, ncol = ncol))
  
  return(plot)

}



make_violin_plot_proteomics = function(data,topn,ncol=5,pname="p.adj.signif"){
  library(gridExtra)
  # initializing list
  list = list()
  #looping through each biomarker
  for(i in topn){
    
    loop = filter(data, Gene == i)
    stats = rstatix::t_test(loop, Abundance ~ Condition) %>% 
    rstatix::add_y_position() %>% 
    rstatix::add_significance()
  
  p1 = loop %>% 
    ggplot(aes(Condition,Abundance))+
    geom_violin()+
    geom_point()+
    theme(axis.text.x = element_text(angle = 90))+
    ggpubr::stat_pvalue_manual(data = stats, label = pname, tip.length = 0.005)+
    ggtitle(i)+
    theme(axis.title.x = element_blank())
  
  list[[i]] = p1
  }
  
  plot = do.call("grid.arrange", c(list, ncol = ncol))
  
  return(plot)

}


## Function to make ROC curves from metabolomics data ##
## Defining a function to make the ROC curves ##
make_roc_metabolomics = function(data,ncol=5){
  
  list = list()
  
  for(i in unique(data$compound_name_cleaned)){
    
    loop = filter(data, compound_name_cleaned == i)
    
    mod = glm(as.factor(condition) ~ norm_value, 
              data = loop,
              family = "binomial") 
    
    predicted <- predict(mod ,family = "binomial", loop, type="response")
    
    #define object to plot
    rocobj <- pROC::roc(loop$condition, predicted)
    
    auc = round(as.numeric(gsub(".*: ", "",rocobj$auc)),2)
    
    #Dealing with long names 
    wrapper <- function(x, ...) 
    {
      paste(strwrap(x, ...), collapse = "\n")
    }
    
    #create ROC plot
    p1 = pROC::ggroc(rocobj)+
      ggtitle(paste0(i," AUC: ",auc))
    
    
    list[[i]] = p1
  }
  plot = do.call("grid.arrange",c(list, ncol = ncol))
  
  return(plot) 
}

make_violin_plot_metabolomics = function(data,topn,ncol=5,pname="p.adj.signif"){
  library(gridExtra)
  # initializing list
  list = list()
  #looping through each biomarker
  for(i in topn){
    
    loop = filter(data, row_id == i)
    stats = rstatix::t_test(loop, norm_value ~ condition) %>% 
      rstatix::add_y_position()%>% 
      rstatix::add_significance()
    
    cn = loop %>% 
      pull(compound_name_cleaned) %>% 
      unique()
    
    p1 = loop %>% 
      ggplot(aes(condition,norm_value))+
      geom_violin()+
      geom_point()+
      theme(axis.text.x = element_text(angle = 90))+
      ggpubr::stat_pvalue_manual(stats, label = pname, tip.length = 0.005)+
      ggtitle(cn)+
      theme(axis.title.x = element_blank())
    
    list[[i]] = p1
  }
  
  plot = do.call("grid.arrange", c(list, ncol = ncol))
  
  return(plot)
  
}





