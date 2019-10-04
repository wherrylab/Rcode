##small Functions for expression data and gene lists
write.gct <- function(dataset, filename) {
  EOL = "\n"
  h1 = c("#1.2", EOL)
  h2 = c(nrow(dataset), "\t", ncol(dataset)-2, EOL)
  h3 = c(paste(colnames(dataset), collapse="\t"), EOL)
  cat(h1, file=filename, append=F, sep="")
  cat(h2, file=filename, append=T, sep="")
  cat(h3, file=filename, append=T, sep="")
  write.table(dataset, file=filename, sep="\t", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE, eol=EOL)
}

gct_subset<-function(gctfile,genelist,filename){
  library(CePa)
  data<-data.frame(read.gct(gctfile))
  sub_data<-data[genelist,]
  sub_data$Gene<-rownames(sub_data)
  sub_data$Desc<-rownames(sub_data)
  n<-ncol(sub_data)
  sub_data<-sub_data[,c(n-1,n,1:(n-2))]
  write.gct(sub_data,filename)
}

phyper_sm_gene_list <- function (Allgenes,Gene1, Gene2, verySig=TRUE ,lt=TRUE){
  # This function is the same is phyper, just allows for more sensible input values
  #Allgenes is list of all unique genes,Gene1 and Gene2 are list of genes to compare
  total<-union(Gene1,Gene2)
  overlap<-intersect(Gene1,Gene2)
  q = length(overlap)
  m = length(Gene1)
  n = length(Allgenes)-length(Gene1)
  k = length(Gene2)
  prob = phyper(q, m, n, k, log.p = verySig, lower.tail=lt)
  if (verySig) return(-prob)
  return(1-prob)
}

geneOverlap_circos <- function (Gene_list){
  overlap_list_all<-list()
  gene_list<-list()
  for(i in 1:length(Gene_list)){
    gene_temp<-as.data.frame(Gene_list[[i]])
    gene_temp$x<-as.numeric(rownames(gene_temp))
    colnames(gene_temp)<-c("Gene","x")
    assign(paste0("Gene",i),gene_temp)
    gene_list[[i]]<-gene_temp
    overlap_list<-list()
    for(j in setdiff(1:length(Gene_list),i)){
      overlap_value<-intersect(Gene_list[[i]],Gene_list[[j]])
      overlap_list[[j]]<-length(overlap_value)
    }
    overlap_list_all[[i]]<-overlap_list
  }
  
  Gene_lst_data<-list()
  for (i in 1:length(gene_list))
  {
    Gene_lst_data[[i]]<-data.frame(gene_list[[i]],factor=paste("f",i,sep=""))
  }
  df <- do.call("rbind", Gene_lst_data)
  require(circlize)
  par(mar = c(1, 1, 1, 1), lwd = 0.2, cex = 0.7)
  circos.par(cell.padding = c(0, 0, 0, 0))
  circos.initialize(factors = df$factor,x=df$x)
  n_cols<-length(unique(df$factor))
  circos.trackPlotRegion(factors = df$factor, y = df$x, bg.col=rand_color(n_cols,transparency=0.1),panel.fun = function(x, y) {
    circos.axis()
  })
  for(f_num in 1:n_cols){
    for(j in setdiff(1:n_cols,i)){
      q<-overlap_list_all[[f_num]][[j]]
      circos.link(paste0("f",f_num),c(1,q), paste0("f",j), c(1,q), col = rand_color(1,transparency = 0.2))
    }
  }
  circos.clear()
}

create_geneset_from_limmares<-function(rdata,filename){
  #The results are already sorted if not sort by p value or fold change accordingly
  rdata<-rdata[order(rdata$adj.P.Val),]
  #rdata<-rdata[order(-rdata$Fold.Change..from.Tho.),]
  #rdata_up<-rdata[rdata$Fold.Change..from.Tho.>0,]
  rdata_up<-rdata[rdata$logFC>0,]
  rdata_dn<-rdata[rdata$logFC<0,]
  rdata_up_top200<-rdata_up[1:200,]
  rdata_dn_top200<-rdata_dn[1:200,]
  genes_up<-as.character(rdata_up_top200$SYMBOL)
  #genes_up<-as.character(rdata_up_top200$Ms_gene_name)
  genes_dn<-as.character(rdata_dn_top200$SYMBOL)
  l1 = c(paste0(filename,"up_genes"),"\t","Sig_genes_up","\t",paste(genes_up,collapse="\t"),"\n")
  #l1 = c(paste0(filename,"up_genes"),"\t","Sig_genes_up","\t",paste(genes_up,collapse="\t"))
  l2 = c(paste0(filename,"dn_genes"),"\t","Sig_genes_dn","\t",paste(genes_dn,collapse="\t"))
  cat(l1, file=filename, append=F, sep="")
  cat(l2, file=filename, append=T, sep="")
}
