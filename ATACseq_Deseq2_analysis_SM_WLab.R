library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("calibrate")
library("FactoMineR")
library("cowplot")
library("rGREAT")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


setwd("~/Desktop/sasi_work/ATAC_seq/External/Haining_ATAC/JG_analysis/combFiles")

#Read in counts data from ATACseq pipeline (JG)
raw.counts<-read.delim("combUnionReadsWithLabels.txt",sep="\t")

cleaned_counts<-na.omit(raw.counts)

cleaned_counts$peak_id<-paste0("Peak",seq(1:nrow(cleaned_counts)))

#Set the peak names  to be the row names for the count data
rownames(cleaned_counts) <- cleaned_counts$peak_id

peak_annot<-cleaned_counts[,c("peak_id","chr","start","end")]
cleaned_counts$peak_id<-NULL

cleaned_counts<-cleaned_counts[,-c(1:3)]

samplenames<-colnames(cleaned_counts)

samplenames_fin<-gsub("X","",samplenames)

colnames(cleaned_counts)<-samplenames_fin

##Genes with zero counts are removed

sel.rn=rowSums(cleaned_counts) != 0
cleaned_counts=cleaned_counts[sel.rn,]

cdata<-read.delim("sampleInfo_haining_nave_cl13_d8.txt",sep="\t",header=TRUE)

rownames(cdata)<-cdata[,1]
cdata[,1]<-NULL

cnts<-cleaned_counts[,as.character(rownames(cdata))]

#write.table(cnts,"ATAC_seq_counts.txt",sep="\t")
#Create a Dataset for further downstream analysis
cds <- DESeqDataSetFromMatrix(cnts,cdata, ~condition)

#Diiferential binding analysis

cds_deseq<-DESeq(cds)

deseq_result<-results(cds_deseq,contrast = c('condition','nave_cl13_d8','arm_d8'))

deseq_result<-data.frame(deseq_result)

deseq_result_annot<-merge(peak_annot,deseq_result,by.x="peak_id",by.y=0)

#write.table(peak_annot[,c(2:4,1)],"all_peaks.bed",sep="\t",row.names = F,col.names = F)

write.table(deseq_result_annot[deseq_result_annot$padj < 0.05 & deseq_result_annot$log2FoldChange > 1,c(2:4,1)],"nave_cl13_d8_vs_arm_d8_UP_peaks.bed",sep="\t",row.names = F,col.names = F)

write.table(deseq_result_annot[deseq_result_annot$padj < 0.05 & deseq_result_annot$log2FoldChange < -1,c(2:4,1)],"nave_cl13_d8_vs_arm_d8_Down_peaks.bed",sep="\t",row.names = F,col.names = F)


write.table(deseq_result_annot[deseq_result_annot$padj < 0.05 & deseq_result_annot$log2FoldChange > 2,c(2:4,1)],"nave_cl13_d8_vs_arm_d8_UP_peaks_2logs.bed",sep="\t",row.names = F,col.names = F)

write.table(deseq_result_annot[deseq_result_annot$padj < 0.05 & deseq_result_annot$log2FoldChange < -2,c(2:4,1)],"nave_cl13_d8_vs_arm_d8_Down_peaks_2logs.bed",sep="\t",row.names = F,col.names = F)


write.table(deseq_result_annot,"Deseq2_all_results_nave_cl13_d8_vs_arm_d8.txt",sep="\t",row.names = T,col.names = T)


## Examine plot of p-values
hist(deseq_result$padj, breaks=50, col="grey")


maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=peak_id, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 6000, 6000, pointsize=20,res = 300,units = "px")
maplot(deseq_result_annot, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=peak_id, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 6000, 6000, pointsize=20,res = 300,units = "px")
volcanoplot(deseq_result_annot, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-6, 6))
dev.off()


RC<-rlogTransformation(cds,blind = FALSE)

DF<-assay(RC)

res.pca = PCA(t(DF), scale.unit=TRUE, graph=F)

pca_coords_df<-data.frame(res.pca$ind$coord)

pca_coords_df$group<-cdata[,1]

pca_coords_df$sample<-rownames(cdata)

theme_custom<-theme(legend.position = "none",axis.title=element_text(face = "bold.italic", color = "black"),axis.text.x = element_text(size=12,angle = 45, hjust = 1,face="bold"),plot.title = element_text(lineheight=.8,hjust=0.5, face="bold"),axis.text.y = element_text(size=12,face="bold"))+ theme(panel.border = element_rect(linetype = "solid",colour = "black"))

p_pca12<-ggplot(pca_coords_df, aes(x=Dim.1, y=Dim.2,color=group,label=sample),height=600, width=600)+theme_linedraw()+theme_custom+geom_text(check_overlap = FALSE,size=3)

p_pca13<-ggplot(pca_coords_df, aes(x=Dim.1, y=Dim.3,color=group,label=sample),height=600, width=600)+theme_linedraw()+geom_text(check_overlap = FALSE,size=3)+theme_custom

p_pca23<-ggplot(pca_coords_df, aes(x=Dim.2, y=Dim.3,color=group,label=sample),height=600, width=600)+theme_linedraw()+geom_text(check_overlap = FALSE,size=3)+theme_custom

legend <- get_legend(ggplot(pca_coords_df, aes(x=Dim.1, y=Dim.2,color=group,label=sample),height=600, width=600)+theme_linedraw()+geom_text())

legend_plot<-ggdraw(legend)

png(filename = "PCA_3_dimensions_Haining_ATACseq.png",res = 300,units="px",width = 3000,height = 3000)
multiplot(p_pca12, p_pca13, p_pca23,legend_plot, cols=2)

dev.off()

write.table(peak_annot[,c(2:4,1)],"Tcell_low_high_all_peaks.bed",sep="\t",row.names = F,col.names = F)

write.table(deseq_result_annot[deseq_result_annot$padj < 0.05 & deseq_result_annot$log2FoldChange > 1,c(2:4,1)],"Tcell_low_UP_peaks.bed",sep="\t",row.names = F,col.names = F)

write.table(deseq_result_annot[deseq_result_annot$padj < 0.05 & deseq_result_annot$log2FoldChange < -1,c(2:4,1)],"Tcell_low_Down_peaks.bed",sep="\t",row.names = F,col.names = F)

write.table(deseq_result_annot,"Deseq2_all_results.txt",sep="\t",row.names = T,col.names = T)


all_peaks_GREAT_annot<-read.delim("20171102-public-3.0.0-spioQ8-mm10-all-region.txt",header = F,skip = 1)

deseq_result_annot2<-merge(deseq_result_annot,all_peaks_GREAT_annot,by.x="peak_id",by.y="V1")

write.table(deseq_result_annot2,"Deseq2_all_results_GREAT_associations.txt",sep="\t",row.names = T,col.names = T)

#Annotation & enrichment using GREAT

bed_up<-read.delim("UP_peaks.bed",header = T)

bed_down<-read.delim("Down_peaks.bed",header = T)


bg_bed<-read.delim("all_peaks.bed",header = F)

#Example GREAT submission
# set.seed(123)
# bed = circlize::generateRandomBed(nr = 1000, nc = 0)
# bed[1:2, ]

job = submitGreatJob(bg_bed,species="mm10",request_interval = 100)

res = plotRegionGeneAssociationGraphs(job)

df<-data.frame(res)
write.csv(df,file = "region_gene_Association_all_peaks.csv")



job_up = submitGreatJob(bed_up,species="mm10",bgChoice ="data",bg=bg_bed,request_interval = 100)

job_down = submitGreatJob(bed_down,species="mm10",bgChoice ="data",bg=bg_bed,request_interval = 100)

png(filename = "Peaks_up_GREAT_annotation.png",res = 300,units="px",width = 5000,height = 3000)
par(mfrow = c(1, 3))
res_up = plotRegionGeneAssociationGraphs(job_up)
dev.off()

df_up<-data.frame(res_up)

write.csv(df_up,file = "region_gene_Association_sig_peaks_up.csv")

png(filename = "Peaks_down_GREAT_annotation.png",res = 300,units="px",width = 5000,height = 3000)
par(mfrow = c(1, 3))
res_down = plotRegionGeneAssociationGraphs(job_down)
dev.off()

df_down<-data.frame(res_down)
write.csv(df_down,file = "region_gene_Association_sig_peaks_down.csv")


tb_up = getEnrichmentTables(job_up,ontology=c(availableOntologies(job_up)))

for ( z in 1:length(tb_up)){
  erich_df<-data.frame(tb_up[[z]])
  write.csv(erich_df,file = paste(gsub(" ","_",names(tb_up)[z]),"peaks_up.csv"))
}

tb_down = getEnrichmentTables(job_down,ontology=c(availableOntologies(job_down)))

for ( z in 1:length(tb_down)){
  erich_df<-data.frame(tb_down[[z]])
  write.csv(erich_df,file = paste(gsub(" ","_",names(tb_down)[z]),"peaks_down.csv"))
}
