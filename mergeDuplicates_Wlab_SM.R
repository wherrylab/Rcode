mergeDuplicates <- function(data, mergeMethod){
  #Function to merge duplicate gene level expression data using one of the mean,median,max methods,the input data should
  # have first column as the gene and rest of the data numeric expression data.
  all_nms <- data[,1]
  colnms <- colnames(data);
  l<-length(colnames(data))
  dup_inx <- duplicated(all_nms);
  
  if(sum(dup_inx) > 0){
    uniq_nms <- all_nms[!dup_inx];
    uniq_data <- data[!dup_inx,,drop=F];
    
    dup_nms <- all_nms[dup_inx];
    uniq_dupnms <- unique(dup_nms);
    uniq_duplen <- length(uniq_dupnms);
    
    for(i in 1:uniq_duplen){
      nm <- uniq_dupnms[i];
      hit_inx_all <- which(all_nms == nm);
      hit_inx_uniq <- which(uniq_nms == nm);
      
      # average the whole sub matrix 
      if(mergeMethod == "mean"){
        uniq_data[hit_inx_uniq,2:l ]<- apply(data[hit_inx_all,2:l,drop=F], 2, mean, na.rm=T);
      }else if(mergeMethod == "median"){
        uniq_data[hit_inx_uniq,2:l ]<- apply(data[hit_inx_all,2:l,drop=F], 2, median, na.rm=T);
      }else if(mergeMethod == "max"){
        uniq_data[hit_inx_uniq,2:l ]<- apply(data[hit_inx_all,2:l,drop=F], 2, max, na.rm=T);
      }else{ # sum
        uniq_data[hit_inx_uniq,2:l ]<- apply(data[hit_inx_all,2:l,drop=F], 2, sum, na.rm=T);
      }
    }
  }
  else{
    uniq_data<-data
  }
  return(uniq_data)
}