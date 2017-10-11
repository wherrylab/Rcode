#Script to filter genes with mean counts below first qu from the PORT output file

#Example: Rscript Filter_PORT_output.R Port_ouput.txt

library(dplyr)
library(tibble)

args = commandArgs(TRUE)
port_input_file = args[1]


read_count =read.delim(file = port_input_file, row.names = NULL) %>% 
  #Strip off "gene:" prefix for gene ID column (if present)
  mutate(id = gsub("gene:", "", id))

data_df = tbl_df(read_count)  %>% 
  column_to_rownames(var = "id") 

l<-ncol(data_df)

genes_mean<-apply(data_df[,-c(l,(l-1))],MARGIN = 1,mean)

genes_mean_summ<-summary(genes_mean)

filtered_ids<-which(genes_mean>genes_mean_summ[2])

data_matrix<-data.frame(data_df)[filtered_ids,]

write.table(data_matrix,gsub(".txt","_filtered.txt",port_input_file),sep="\t",quote = F,row.names = F)
