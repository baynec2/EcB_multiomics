## Function for making ROC plots using proteomics data ##

make_roc = function(data,topn){
  library(gridExtra)
  
  df = data.frame()
  
  for(i in topn){
    
    loop = filter(data, Feature == i)
    
    mod = glm(as.factor(condition) ~ Abundance, 
              data = data,
              family = "binomial") 
    
    predicted <- predict(mod ,family = "binomial", loop, type="response")
    
    #define object to plot
    rocobj <- pROC::roc(loop$condition, predicted)
    
    auc = round(as.numeric(gsub(".*: ", "",rocobj$auc)),2)
    #create ROC plot
    p1 = pROC::ggroc(rocobj)+
      ggtitle(paste0(i," AUC: ",auc))
    
    
   p1 =  tibble(Feature = i, sensitivity = p1$data$sensitivity, specificity = p1$data$specificity)
    
    
    df = rbind(df, p1)
  }
  return(df)
}


make_roc_m = function(data,topn){
  library(gridExtra)
  
  df = data.frame()
  
  for(i in topn){
    
    loop = filter(data, Feature == i)
    
    cn = unique(loop$compound_name_cleaned)
    
    mod = glm(as.factor(condition) ~ Abundance, 
              data = data,
              family = "binomial") 
    
    predicted <- predict(mod ,family = "binomial", loop, type="response")
    
    #define object to plot
    rocobj <- pROC::roc(loop$condition, predicted)
    
    auc = round(as.numeric(gsub(".*: ", "",rocobj$auc)),2)
    #create ROC plot
    p1 = pROC::ggroc(rocobj)+
      ggtitle(paste0(cn," AUC: ",auc))
    
    
    p1 =  tibble(Feature = i,compound_name = cn, sensitivity = p1$data$sensitivity, specificity = p1$data$specificity)
    
    
    df = rbind(df, p1)
  }
  return(df)
}




