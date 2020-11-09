module load r r-ape r-scales r-ggplot2 r-gplots r-deseq2
#####################################################
################# File Discription ##################
#####################################################

# 'totalReadCounts.txt' -  counts of total gene expression, genes (row) X 26 RNA-seq libraries (column)
# 'MaxxaAlleleReadCounts.txt' - for 12 F1 libraries (column), allele-specific read count assigned to Maxxa alleles were presented for these genes contain diagnostic SNPs.
# 'TX2094AlleleReadCounts.txt' - same as above but for TX2094 allele.
# 'TM1.nuclNmt.transcript.len.txt' - gene length for all genes


#####################################################
########## Finalize the sample inclusion #############
#####################################################

# load total counts
count <-read.table("totalReadCounts.txt", header=TRUE, sep="\t")
# requiring a read depth of 5 per genes, how many expressed genes in each accession?
getTotalExpressed=function(x,depth=5)
{
  length(which(rowSums(x)>=depth))
}
getTotalExpressed(count)
getTotalExpressed(count) # 71065
getTotalExpressed(count[,grep("Maxxa",names(count))]) # 66205
getTotalExpressed(count[,grep("TX2094",names(count))]) # 66745
getTotalExpressed(count[,grep("F1",names(count))]) # 69565
getTotalExpressed(count[,grep("F1.*[1-3]$",names(count))]) # 64836
getTotalExpressed(count[,grep("F1.*[4-6]$",names(count))]) # 66412
# Union of expressed fiber genes
length(which(rowSums(count[,grep("Maxxa",names(count))])>=5 | rowSums(count[,grep("TX2094",names(count))])>=5 | rowSums(count[,grep("F1.*[1-3]$",names(count))])>=5 |rowSums(count[,grep("F1.*[4-6]$",names(count))])>=5))
# 69532
length(which(rowSums(count[,grep("Maxxa",names(count))])>=5 | rowSums(count[,grep("TX2094",names(count))])>=5 | rowSums(count[,grep("F1.*[1-3]$",names(count))])>=5 |rowSums(count[,grep("F1.*[4-6]$",names(count))])>=5))/nrow(count)
# 0.9224565
length(which(rowSums(count[,grep("Maxxa",names(count))])>=5 & rowSums(count[,grep("TX2094",names(count))])>=5 & rowSums(count[,grep("F1.*[1-3]$",names(count))])>=5 &rowSums(count[,grep("F1.*[4-6]$",names(count))])>=5))
# overlap set 62379



######## This appears to be Jing modifying an info file (sample, size, accession, treatment) which I already have
#Do not ever for any reason EVER EVER use rdata objects they fucking suck and are stupid and dumb and don't even fucking work those fuckers
## load total counts from bam idxstat
#load("counts.rdata")->l
#l
rinfo<- read.table("renameInfo.txt", header=TRUE)
total<- read.table("total.txt", header=TRUE)
###### I am not getting rdata objects to work. I have save the following as their own objects. These were generated in PCAanalysis.r
## "info"       "total"      "count.rld"  "count.log2"
##This is weeding out lines from rinfo based on accession. Jing used to filter by accession, treatment, etc. I have a file called renameInfo that has only the samples I want
#select<-which( (rinfo$Accession=="Maxxa"|rinfo$Accession=="TX2094"|rinfo$Accession=="F1") )
#total<-total[,select]
#info<-info[select,]
#info$ID <- gsub("Yuc", "TX2094",paste(info$genome,info$dpa,info$rep, sep="."))
## need to match the sample orders
#names(total)
#info$ID
#total<-total[,match(gsub("x$","",names(count)),info$ID)]
#info<- info[match(gsub("x$","",names(count)),info$ID),]
######## This appears to be Jing modifying an info file (sample, size, accession, treatment) which I already have

# double check

#verify that rinfo rownames are the same as column names in count
rownames(rinfo) == gsub("x$","",names(count))
rownames(rinfo)<-names(count)
#This line excludes parts of samples name. It doesn't seem to be needed for my purposes
#rinfo$exclude<-factor(paste0(rinfo$genome,gsub(".*[0-9]","_",rownames(rinfo))))
# match gene order
total<-total[rownames(count),]
unique(rownames(total)==rownames(count))
# not exactly the same, because hylite missed the last gene, but proceed with PCA plots
#colSums(total)-colSums(count)

## deseq normalization
library(DESeq2)
dds <- DESeqDataSetFromMatrix( countData = count, colData = rinfo, design = ~Treatment)
# rlog transformation, note that with defalt blind=TRUE, design is actually ~1
rld <- rlog(dds)
count.rld <- as.data.frame(assay(rld))
names(count.rld) == names(total)
names(count.rld) == rownames(rinfo)

# also possible to perform custom transformation:
dds <- estimateSizeFactors(dds)
# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1), colData=colData(dds))
count.log2 <- as.data.frame(assay(se))

# RPKM transformation, Reads Per Kilobase of transcript per Million mapped reads.
# RPKM= 10^9 * number of mapped reads in gene
#             ----------------------------------------------------
#             total reads *  exon length
## rpkm normalization
len  <-read.table("TM1.nuclNmt.transcript.len.txt",  sep="\t")
geneLen <- len$V2
names(geneLen) <- len$V1
geneLen <- geneLen[rownames(count)]
rpkm <-sweep(sweep(count*10^9, 2 ,colSums(count), FUN="/"), 1, geneLen, FUN ="/" )
rpkm.log2 <-log2(rpkm+1)

# my custom function
library(ggplot2)
library(scales)
library(ape)
plotGrouping <- function(norm_log, color, shape, text, tip, save = "plotGrouping.pdf"){
  # norm<-sweep(total,2,info$lib_size,"/")*10^6
  # norm_log <- log2(norm+1)
  pca=prcomp(t(norm_log))
  dat = as.data.frame(pca$x)
  proportion<-summary(pca)$importance[2,1:2]*100
  proportion<-paste0(names(proportion)," (", proportion, "%)")
  p<-ggplot(aes(PC1, PC2, color=color,    shape=shape),data=dat) + geom_point() +xlab(proportion[1]) + ylab(proportion[2])
  pdf(save)
  print( p + geom_text(aes_string(x = "PC1", y = "PC2", label = text), color="grey", hjust = 0, nudge_x = 0.09) )
  
  hc<-hclust( dist(t(norm_log)) )
  tre <- as.phylo(hc)
  tre$tip.label <- as.character(tip)
  # try to match ggplot color: library(scale);
  # show_col(col4<-hue_pal()(4))
  tipCol <- color
  levels(tipCol) <- hue_pal()(nlevels(tipCol))
  plot(tre,tip.col =as.character(tipCol),  type="unrooted",cex=0.6, no.margin=TRUE)
  plot(tre,tip.col =as.character(tipCol), type="fan",cex=0.6, no.margin=TRUE)
  dev.off()   }

# plots
plotGrouping(count.rld, color=rinfo$Treatment, shape=rinfo$Accession, tip=rownames(rinfo), text=rinfo$sample, save = "plotGrouping.rld.pdf")
plotGrouping(count.log2, color=rinfo$Treatment, shape=rinfo$Accession, tip=rownames(rinfo), text=rinfo$sample, save = "plotGrouping.log2.pdf")
plotGrouping(rpkm.log2, color=rinfo$Treatment, shape=rinfo$Accession, tip=rownames(rinfo), text=rinfo$sample, save = "plotGrouping.log2rpkm.pdf")


#####################################################
################# DGE for all genes #################
#####################################################
pairwiseDE<-function(dds, contrast,savePath)
{
  # DE analysis
  print(contrast)
  ddsPW <-dds[,dds$condition %in% contrast]
  ddsPW$condition<-droplevels(ddsPW$condition)
  res <- results(DESeq(ddsPW), contrast = c("condition",contrast))  # make sure in order
  print( summary(res,alpha=.05) ) # print results
  write.table(res, file=paste(savePath,"DE/",paste(contrast, collapse="vs"),".txt", sep=""), sep="\t")
}
getSig<-function(res,fc.threshold=0,direction=NULL){
  sig<- res[res$padj<0.05 & !is.na(res$padj) & abs(res$log2FoldChange)>=fc.threshold,]
  if(is.null(direction)){
    n<-nrow(sig)
  }else if(direction=="up"){
    n<-nrow(sig[sig$log2FoldChange>0,])
  }else if(direction=="down"){
    n<-nrow(sig[sig$log2FoldChange<0,])
  }
  return(n)
}

#new dataframe from old dataframe based on rownames from rinfo
countT<-count[,!grepl("x$",rownames(rinfo))]
coldataT <- rinfo[!grepl("x$",rownames(rinfo)),]
#creates new column called condition which is a combination of accession and daylength treatment
coldataT$condition<-factor(paste(coldataT$Accession,coldataT$Treatment, sep='.'))
#double checks that names match
names(countT)==rownames(coldataT)
coldataT$lib_size - colSums(countT)  # slight difference
coldataT$lib_size <- colSums(countT) # modify

# differential expression to get expression divergence
library(DESeq2)
dds <- DESeqDataSetFromMatrix( countData = countT, colData = coldataT, design = ~ condition)
# pairwise deseq workflow
batch<- rbind(
  c("Maxxa.SD7", "TX2094.SD7" ),
  c("Maxxa.SD7", "F1.SD7"),
  c("TX2094.SD7", "F1.SD7"),
  
  c("Maxxa.SD5", "TX2094.SD5" ),
  c("Maxxa.SD5", "F1.SD5"),
  c("TX2094.SD5", "F1.SD5"),
  
  c("Maxxa.LD9", "TX2094.LD9" ),
  c("Maxxa.LD9", "F1.LD9"),
  c("TX2094.LD9", "F1.LD9"),
  
  c("Maxxa.LD7", "TX2094.LD7" ),
  c("Maxxa.LD7", "F1.LD7"),
  c("TX2094.LD7", "F1.LD7"),
  
  c("Maxxa.SD7", "Maxxa.SD5" ),
  c("TX2094.SD7", "TX2094.SD5" ),
  c("F1.SD7", "F1.SD5" ),

  c("Maxxa.LD7", "Maxxa.LD9" ),
  c("TX2094.LD7", "TX2094.LD9" ),
  c("F1.LD7", "F1.LD9" )  )

# make a "DE" folder
system("mkdir DE")
apply(batch,1,function(x) pairwiseDE(dds,x,savePath = ""))

# analyze MxT and TxM
#assign rep numbers
coldataT2<-coldataT
# a vector that indicates replicate 
coldataT2$rep <- c(1,2,3,4,5,6,1,2,3,4,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,1,2,3,3,1,2,3,1,2,3,1,2,2,3,1,2,3,1,2,3,1,2,3)
#sorts F1 samples by maternal parent
mt<-which(coldataT2$Accession=="F1" & coldataT2$rep %in% 1:3)
tm<-which(coldataT2$Accession=="F1" & coldataT2$rep %in% 4:6)
coldataT2$condition<-as.character(coldataT2$condition)
#Renames F1 as MxT or TxM based on parentage
coldataT2$condition[mt]<-gsub("F1","MxT",coldataT2$condition[mt])
coldataT2$condition[tm]<-gsub("F1","TxM",coldataT2$condition[tm])
coldataT2$condition <- as.factor(coldataT2$condition)
dds2 <- DESeqDataSetFromMatrix( countData = countT, colData = coldataT2, design = ~ condition)

#new pairwise based on separating F1s based on parentage 
batch2<- rbind(
  c("MxT.SD7", "TxM.SD7" ),
  c("Maxxa.SD7", "MxT.SD7"),
  c("Maxxa.SD7", "TxM.SD7"),
  c("TX2094.SD7", "MxT.SD7"),
  c("TX2094.SD7", "TxM.SD7"),
  
  c("MxT.SD5", "TxM.SD5" ),
  c("Maxxa.SD5", "MxT.SD5"),
  c("Maxxa.SD5", "TxM.SD5"),
  c("TX2094.SD5", "MxT.SD5"),
  c("TX2094.SD5", "TxM.SD5"),
  
  c("MxT.LD7", "TxM.LD7" ),
  c("Maxxa.LD7", "MxT.LD7"),
  c("Maxxa.LD7", "TxM.LD7"),
  c("TX2094.LD7", "MxT.LD7"),
  c("TX2094.LD7", "TxM.LD7"),
  
  c("MxT.LD9", "TxM.LD9" ),
  c("Maxxa.LD9", "MxT.LD9"),
  c("Maxxa.LD9", "TxM.LD9"),
  c("TX2094.LD9", "MxT.LD9"),
  c("TX2094.LD9", "TxM.LD9"),
  
  c("MxT.LD7", "MxT.LD9"),
  c("TxM.LD7", "TxM.LD9"),
  
  c("MxT.SD7", "MxT.SD5" ),
  c("TxM.SD7", "TxM.SD5" ))

apply(batch2,1,function(x) pairwiseDE(dds2,x,savePath = ""))


# print out results in the order of batch
fileL<- list.files("DE")
fileL<-paste0(batch[,1],"vs",batch[,2],".txt")
fileL<-c(fileL, paste0(batch2[,1],"vs",batch2[,2],".txt"))
sigT<-c("sample 1","sample 2","DE (q<0.05)","1>2","2>1")
for(file in fileL)
{
  res<-read.table(paste0("DE/",file),sep="\t",header=TRUE)
  sigRes <- c(unlist(strsplit(gsub(".txt","",file),split="vs") ), getSig(res),getSig(res,direction="up"),getSig(res,direction="down"))
  sigT<-rbind(sigT,sigRes)
}
T<-as.data.frame(sigT[-1,], row.names=FALSE)
names(T)<-sigT[1,]
print(T)

library(gplots)
pdf("DE/checkDE.pdf")
textplot(T)
mtext("DE analysis result summary")
dev.off()

# rename this "DE" folder as "DE0"
#system("mv DE/ DE-66610-mar2019/")

# save total datasets
write.table(countT, file="Genes.raw.count.txt", sep="\t")
write.table(coldataT2, file="Genes.sample.info.txt", sep="\t")

plotGrouping <- function(norm_log, color, shape, text, tip, save = "plotGrouping.pdf"){
  # norm<-sweep(total,2,info$lib_size,"/")*10^6
  # norm_log <- log2(norm+1)
  pca=prcomp(t(norm_log))
  dat = as.data.frame(pca$x)
  proportion<-summary(pca)$importance[2,1:2]*100
  proportion<-paste0(names(proportion)," (", proportion, "%)")
  p<-ggplot(aes(PC1, PC2, color=color,    shape=shape),data=dat) + geom_point() +xlab(proportion[1]) + ylab(proportion[2])
  pdf(save)
  print( p + geom_text(aes_string(x = "PC1", y = "PC2", label = text), color="grey", hjust = 0, nudge_x = 0.09) )
  
  hc<-hclust( dist(t(norm_log)) )
  tre <- as.phylo(hc)
  tre$tip.label <- as.character(tip)
  # try to match ggplot color: library(scale);
  # show_col(col4<-hue_pal()(4))
  tipCol <- color
  levels(tipCol) <- hue_pal()(nlevels(tipCol))
  plot(tre,tip.col =as.character(tipCol),  type="unrooted",cex=0.6, no.margin=TRUE)
  plot(tre,tip.col =as.character(tipCol), type="fan",cex=0.6, no.margin=TRUE)
  dev.off()   }

# normalization with rpkm, len.txt comes out of hylite
len  <-read.table("TM1.nuclNmt.transcript.len.txt",  sep="\t")
geneLen <- len$V2
names(geneLen) <- len$V1
geneLen <- geneLen[rownames(countT)]
rpkm <-sweep(sweep(countT*10^9, 2 ,coldataT$lib_size, FUN="/"), 1, geneLen, FUN ="/" )
rpkm.log2 <-log2(rpkm+1)
plotGrouping(rpkm.log2, color=coldataT$Treatment, shape=coldataT$Accession, tip=rownames(coldataT),text=coldataT$condition, save = "plotGrouping.66610.log2rpkm.pdf")

# save total datasets
write.table(rpkm, file="Genes66610.log2rpkm.txt", sep="\t")
write.table(geneLen, file="Genes66610.gene.length.txt", sep="\t")

####################waiting to hear from jing
############# read annotation_info.txt from lss/research/jfw-lab/Projects/ChenCottons/##whichever version you are using
############x<-read.table(,sep="\t",comment.char="", quote="\"",header=TRUE)
############
############annot<-x[grep("[.]1$",x$transcriptName),]
############dim(annot)
############annot<-annot[,c(3,5:13)]
############desc<-annot[,c("transcriptName","arabi.defline")]
############y<-read.table("mt.annotation.txt",header=TRUE,sep="\t",comment.char="", quote="\"")
############y<-y[,1:2]
############names(y)<-names(desc)
############desc<-rbind(desc,y[,1:2])
############z<-merge(desc, annot, by="transcriptName",all.x=TRUE)
############rownames(z)<-z$transcriptName
############annotation<-z[rownames(rpkm),]
############write.table(annotation, file="Genes66610.transcrip.description.txt", sep="\t")
############

#####################################################
############# Finalize gene inclusion ###############
#####################################################

# import allelic read count tables
maxxa  <-read.table("MaxxaAlleleReadCounts.txt", header=TRUE, sep="\t")
tx2094 <-read.table("TX2094AlleleReadCounts.txt", header=TRUE, sep="\t")
# check column order
names(maxxa)[-1] == names(countT)[1:12]
names(maxxa) == names(tx2094)
# combine
#Add .maxxa to maxxa names and .tx2094 to tx2094 names and then merge files
names(maxxa)[-1] <- paste0(names(maxxa)[-1],".maxxa")
names(tx2094)[-1] <- paste0(names(tx2094)[-1],".tx2094")
allele <- merge(maxxa, tx2094, by ="gene")

# prepare for DEseq2
countA <- allele[,-1]
rownames(countA) <- allele$gene
#taking the F1 rows (in this case 23) from coldataT and duplicating them so one set will be maxxa alleles and the other set tx2094 alleles
coldataA <- rbind(coldataT[1:23,], coldataT[1:23,])
coldataA$allele <-gsub(".*[.]", "",names(countA))
coldataA$condition <-paste0(coldataA$allele, ".", coldataA$condition)
# check percentage of allele read count in total
c(colSums(maxxa[,-1]),colSums(tx2094[,-1]))/coldataA$lib_size # ~7.6%
allele_lib_size <- c(colSums(maxxa[,-1]),colSums(tx2094[,-1]))

# add countN
names(countT) <- paste0(names(countT),".total")
rownames(coldataT2) <- names(countT)
coldataT2$allele <-"total"

# allele + total
countAll <- cbind(countA, countT[rownames(countA),])
coldataAll <- rbind(coldataA[1:5], coldataT)
coldataAll$condition2 <-coldataAll$condition
coldataAll$rep <- c(1,2,3,4,5,6,1,2,3,4,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,6,1,2,3,4,5,6,1,2,3,4,5,6,
	1,2,3,4,5,6,1,2,3,4,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,1,2,3,3,1,2,3,1,2,3,1,2,2,3,1,2,3,1,2,3,1,2,3)
coldataAll$condition2[coldataAll$Accession=="F1"&coldataAll$rep %in%1:3]<-gsub("F1","MxT",coldataAll$condition2[coldataAll$Accession=="F1"&coldataAll$rep %in%1:3])
coldataAll$condition2[coldataAll$Accession=="F1"&coldataAll$rep %in%4:6]<-gsub("F1","TxM",coldataAll$condition2[coldataAll$Accession=="F1"&coldataAll$rep %in%4:6])

#Establishing Jing's function called plotGrouping# my custom function
library(ggplot2)
library(scales)
library(ape)
plotGrouping <- function(norm_log, color, shape, text, tip, save = "plotGrouping.pdf"){
  # norm<-sweep(total,2,info$lib_size,"/")*10^6
  # norm_log <- log2(norm+1)
  pca=prcomp(t(norm_log))
  dat = as.data.frame(pca$x)
  proportion<-summary(pca)$importance[2,1:2]*100
  proportion<-paste0(names(proportion)," (", proportion, "%)")
  p<-ggplot(aes(PC1, PC2, color=color,    shape=shape),data=dat) + geom_point() +xlab(proportion[1]) + ylab(proportion[2])
  pdf(save)
  print( p + geom_text(aes_string(x = "PC1", y = "PC2", label = text), color="grey", hjust = 0, nudge_x = 0.09) )
  
  hc<-hclust( dist(t(norm_log)) )
  tre <- as.phylo(hc)
  tre$tip.label <- as.character(tip)
  # try to match ggplot color: library(scale);
  # show_col(col4<-hue_pal()(4))
  tipCol <- color
  levels(tipCol) <- hue_pal()(nlevels(tipCol))
  plot(tre,tip.col =as.character(tipCol),  type="unrooted",cex=0.6, no.margin=TRUE)
  plot(tre,tip.col =as.character(tipCol), type="fan",cex=0.6, no.margin=TRUE)
  dev.off()   }


# normalization with rpkm
len<- read.table("TM1.nuclNmt.transcript.len.txt",  sep="\t")
geneLen <- len$V2
names(geneLen) <- len$V1
geneLen <- geneLen[rownames(countAll)]
rpkm <-sweep(sweep(countAll*10^9, 2 ,coldataAll$lib_size, FUN="/"), 1, geneLen, FUN ="/" )
rpkm.log2 <-log2(rpkm+1)
plotGrouping(rpkm.log2, color=coldataAll$Treatment, shape=coldataAll$allele, tip=rownames(coldataAll),text=coldataAll$Accession, save = "plotGrouping.log2rpkm.pdf")


#####Need annotation info from Jing
# read annotation
#x<-read.table("Ghirsutum_458_v1.1.annotation_info.txt",sep="\t",comment.char="", quote="\"",header=TRUE)
#annot<-x[grep("[.]1$",x$transcriptName),]
#dim(annot)
#annot<-annot[,c(3,5:13)]
#desc<-annot[,c("transcriptName","arabi.defline")]
#y<-read.table("mt.annotation.txt",header=TRUE,sep="\t",comment.char="", quote="\"")
#y<-y[,1:2]
#names(y)<-names(desc)
#desc<-rbind(desc,y[,1:2])
#z<-merge(desc, annot, by="transcriptName",all.x=TRUE)
#rownames(z)<-z$transcriptName
#annotation<-z[rownames(rpkm),]
#
#write.table(countAll, file="Genes27816.raw.count.txt", sep="\t")
#write.table(coldataAll, file="Genes27816.sample.info.txt", sep="\t")
#write.table(rpkm, file="Genes27816.rpkm.txt", sep="\t")
#write.table(geneLen, file="Genes27816.gene.length.txt", sep="\t")
#write.table(annotation, file="Genes27816.transcrip.description.txt", sep="\t")



#####################################################
############### Differential expression #############
#####################################################
#Check that rownames in coldataAll match column names in countAll
rownames(coldataAll)==names(countAll)
names(countAll)<-rownames(coldataAll)
coldataAll$condition <- as.factor(coldataAll$condition)

# differential expression to get allelic expression divergence
library(DESeq2)
dds <- DESeqDataSetFromMatrix( countData = countAll, colData = coldataAll, design = ~ condition)


# if using pre-existing size factors; not necessary
sizeFactors(dds) <- coldataAll$lib_size/mean(coldataAll$lib_size)

# pairwise deseq workflow, c(maxxa, tx2094), ratio calculated as tx2094/maxxa
# UPDATE WITH CURRENT CONDITION NAMES
batch<- rbind(
  c("maxxa.F1.LD7", "tx2094.F1.LD7" ),
  c("Maxxa.LD9", "tx2094.F1.LD9" ),
  c("Maxxa.LD7", "TX2094.LD7" ),
  c("Maxxa.LD9", "TX2094.LD9" ),
  c("TX2094.LD7", "F1.LD7"),
  c("TX2094.LD9", "F1.LD9"),
  c("Maxxa.LD7", "F1.LD7"),
  c("Maxxa.LD9", "F1.LD9"), 
  c("maxxa.F1.SD7", "tx2094.F1.SD7" ),
  c("maxxa.F1.SD5", "tx2094.F1.SD5" ),
  c("Maxxa.SD7", "TX2094.SD7" ),
  c("Maxxa.SD5", "TX2094.SD5" ),
  c("TX2094.SD7", "F1.SD7"),
  c("TX2094.SD5", "F1.SD5"),
  c("Maxxa.SD7", "F1.SD7"),
  c("Maxxa.SD5", "F1.SD5") )
#mkdir DE if it hasn't already been made
system("mkdir DE")
#run pairwise using combos from pairwise
apply(batch,1,function(x) pairwiseDE(dds,x,savePath = ""))

### Get allelic expression divergence, test for B
BLD7<-read.table("DE/maxxa.F1.LD7vstx2094.F1.LD7.txt", header=TRUE, sep="\t")
BLD9<-read.table("DE/maxxa.F1.LD9vstx2094.F1.LD9.txt", header=TRUE, sep="\t")
BSD7<-read.table("DE/maxxa.F1.SD7vstx2094.F1.SD7.txt", header=TRUE, sep="\t")
BSD5<-read.table("DE/maxxa.F1.SD5vstx2094.F1.SD5.txt", header=TRUE, sep="\t")
### p score less than 0.05 and is not NA
dim(BLD7[BLD7$padj<0.05 & !is.na(BLD7$padj),]) # 4388
dim(BLD9[BLD9$padj<0.05 & !is.na(BLD9$padj),]) # 4780
dim(BSD7[BSD7$padj<0.05 & !is.na(BSD7$padj),]) # 1789
dim(BSD5[BSD5$padj<0.05 & !is.na(BSD5$padj),]) # 5824
### total unique p<0.05 across a day length treatment.
length(unique(c(which(BLD7$padj<0.05 & !is.na(BLD7$padj)),which(BLD9$padj<0.05 & !is.na(BLD9$padj))))) #6178
length(unique(c(which(BSD7$padj<0.05 & !is.na(BSD7$padj)),which(BSD5$padj<0.05 & !is.na(BSD5$padj))))) #6472


#### Get parental expression divergence, test for A=0
## read in each treatment table for the crosses 
ALD7<-read.table("DE/Maxxa.LD7vsTX2094.LD7.txt", header=TRUE, sep="\t")
ALD9<-read.table("DE/Maxxa.LD9vsTX2094.LD9.txt", header=TRUE, sep="\t")
ASD7<-read.table("DE/Maxxa.SD7vsTX2094.SD7.txt", header=TRUE, sep="\t")
ASD5<-read.table("DE/Maxxa.SD5vsTX2094.SD5.txt", header=TRUE, sep="\t")


dim(ALD7[ALD7$padj<0.05 & !is.na(ALD7$padj),]) # 1730
dim(ALD9[ALD9$padj<0.05 & !is.na(ALD9$padj),]) # 3119
dim(ASD7[ASD7$padj<0.05 & !is.na(ASD7$padj),]) # 4180
dim(ASD5[ASD5$padj<0.05 & !is.na(ASD5$padj),]) # 3191

length(unique(c(which(ALD7$padj<0.05 & !is.na(ALD7$padj)),which(ALD9$padj<0.05 & !is.na(ALD9$padj))))) #3847
length(intersect(which(ASD7$padj<0.05 & !is.na(ASD7$padj)),which(ASD5$padj<0.05 & !is.na(ASD5$padj)))) #1360

# analyze MxT and TxM
coldataAll2<-coldataAll
coldataAll2$condition <- coldataAll$condition2
coldataAll2$condition<-as.factor(coldataAll2$condition)
dds2 <- DESeqDataSetFromMatrix( countData = countAll, colData = coldataAll2, design = ~ condition)
rld <- rlog(dds2)
dists <- dist(t(assay(rld)))
# pairwise deseq workflow, c(maxxa, tx2094), ratio calculated as tx2094/maxxa
batchMTLD<- rbind(
  c("maxxa.MxT.LD7", "tx2094.MxT.LD7" ),
  c("maxxa.MxT.LD9", "tx2094.MxT.LD9" ),
  c("TX2094.LD7", "MxT.LD7"),
  c("TX2094.LD9", "MxT.LD9"),
  c("Maxxa.LD7", "MxT.LD7"),
  c("Maxxa.LD9", "MxT.LD9") )
batchTMLD<- rbind(
  c("maxxa.TxM.LD7", "tx2094.TxM.LD7" ),
  c("maxxa.TxM.LD9", "tx2094.TxM.LD9" ),
  c("TX2094.LD7", "TxM.LD7"),
  c("TX2094.LD9", "TxM.LD9"),
  c("Maxxa.LD7", "TxM.LD7"),
  c("Maxxa.LD9", "TxM.LD9") )
batchffLD<-rbind(
  c("MxT.LD7", "TxM.LD7"),
  c("MxT.LD9", "TxM.LD9"),
  c("maxxa.MxT.LD7", "maxxa.TxM.LD7"),
  c("maxxa.MxT.LD9", "maxxa.TxM.LD9"),
  c("tx2094.MxT.LD7", "tx2094.TxM.LD7"),
  c("tx2094.MxT.LD9", "tx2094.TxM.LD9"))
  batchMTSD<- rbind(
  c("maxxa.MxT.SD7", "tx2094.MxT.SD7" ),
  c("maxxa.MxT.SD5", "tx2094.MxT.SD5" ),
  c("TX2094.SD7", "MxT.SD7"),
  c("TX2094.SD5", "MxT.SD5"),
  c("Maxxa.SD7", "MxT.SD7"),
  c("Maxxa.SD5", "MxT.SD5") )
batchTMSD<- rbind(
  c("maxxa.TxM.SD7", "tx2094.TxM.SD7" ),
  c("maxxa.TxM.SD5", "tx2094.TxM.SD5" ),
  c("TX2094.SD7", "TxM.SD7"),
  c("TX2094.SD5", "TxM.SD5"),
  c("Maxxa.SD7", "TxM.SD7"),
  c("Maxxa.SD5", "TxM.SD5") )
batchffSD<-rbind(
  c("MxT.SD7", "TxM.SD7"),
  c("MxT.SD5", "TxM.SD5"),
  c("maxxa.MxT.SD7", "maxxa.TxM.SD7"),
  c("maxxa.MxT.SD5", "maxxa.TxM.SD5"),
  c("tx2094.MxT.SD7", "tx2094.TxM.SD7"),
  c("tx2094.MxT.SD5", "tx2094.TxM.SD5"))

### Jing's pairwise DEfunction
pairwiseDE<-function(dds, contrast,savePath)
{
  # DE analysis
  print(contrast)
  ddsPW <-dds[,dds$condition %in% contrast]
  ddsPW$condition<-droplevels(ddsPW$condition)
  res <- results(DESeq(ddsPW), contrast = c("condition",contrast))  # make sure in order
  print( summary(res,alpha=.05) ) # print results
  write.table(res, file=paste(savePath,"DE/",paste(contrast, collapse="vs"),".txt", sep=""), sep="\t")
}
getSig<-function(res,fc.threshold=0,direction=NULL){
  sig<- res[res$padj<0.05 & !is.na(res$padj) & abs(res$log2FoldChange)>=fc.threshold,]
  if(is.null(direction)){
    n<-nrow(sig)
  }else if(direction=="up"){
    n<-nrow(sig[sig$log2FoldChange>0,])
  }else if(direction=="down"){
    n<-nrow(sig[sig$log2FoldChange<0,])
  }
  return(n)
}

apply(batchMTLD,1,function(x) pairwiseDE(dds2,x,savePath = ""))
apply(batchTMLD,1,function(x) pairwiseDE(dds2,x,savePath = ""))
apply(batchffLD,1,function(x) pairwiseDE(dds2,x,savePath = ""))
apply(batchMTSD,1,function(x) pairwiseDE(dds2,x,savePath = ""))
apply(batchTMSD,1,function(x) pairwiseDE(dds2,x,savePath = ""))
apply(batchffSD,1,function(x) pairwiseDE(dds2,x,savePath = ""))


### Get allelic expression divergence, test for B=0
## mt samples are maxxa x tx 2094 and tm samples are tx2094 x maxxa parentage
B.mtLD7<-read.table("DE/maxxa.MxT.LD7vstx2094.MxT.LD7.txt", header=TRUE, sep="\t")
B.mtLD9<-read.table("DE/maxxa.MxT.LD9vstx2094.MxT.LD9.txt", header=TRUE, sep="\t")
B.mtSD7<-read.table("DE/maxxa.MxT.SD7vstx2094.MxT.SD7.txt", header=TRUE, sep="\t")
B.mtSD5<-read.table("DE/maxxa.MxT.SD5vstx2094.MxT.SD5.txt", header=TRUE, sep="\t")
dim(B.mtLD7[B.mtLD7$padj<0.05 & !is.na(B.mtLD7$padj),]) #2931 6
dim(B.mtLD9[B.mtLD9$padj<0.05 & !is.na(B.mtLD9$padj),]) #2882 6
dim(B.mtSD7[B.mtSD7$padj<0.05 & !is.na(B.mtSD7$padj),]) #398 6
dim(B.mtSD5[B.mtSD5$padj<0.05 & !is.na(B.mtSD5$padj),]) #2514 6

B.tmLD7<-read.table("DE/maxxa.TxM.LD7vstx2094.TxM.LD7.txt", header=TRUE, sep="\t")
B.tmLD9<-read.table("DE/maxxa.TxM.LD9vstx2094.TxM.LD9.txt", header=TRUE, sep="\t")
B.tmSD7<-read.table("DE/maxxa.TxM.SD7vstx2094.TxM.SD7.txt", header=TRUE, sep="\t")
B.tmSD5<-read.table("DE/maxxa.TxM.SD5vstx2094.TxM.SD5.txt", header=TRUE, sep="\t")
dim(B.tmLD7[B.tmLD7$padj<0.05 & !is.na(B.tmLD7$padj),]) #1425 6
dim(B.tmLD9[B.tmLD9$padj<0.05 & !is.na(B.tmLD9$padj),]) #1941 6
dim(B.tmSD7[B.tmSD7$padj<0.05 & !is.na(B.tmSD7$padj),]) #218 6
dim(B.tmSD5[B.tmSD5$padj<0.05 & !is.na(B.tmSD5$padj),]) #2854 6


# important: A10, B10, B.mt10, B.tm10, A20, B20, B.mt20, B.tm20,
save(ALD7, ALD9, ASD7, ASD5, BLD7, BLD9, BSD7, BSD5, B.mtLD7, B.mtLD9, B.mtSD7, B.mtSD5, B.tmLD7 , B.tmLD9, B.tmSD7,  B.tmSD5, file="cistrans.rdata")

system("mv DE/ DE-27816-mar2019/")



#####################################################
############### Visualize DE results ################
#####################################################

#fileL<- list.files("DE27816")
fileL<-paste0(batch[,1],"vs",batch[,2],".txt")
fileL<-c(fileL, paste0(batchMTLD[,1],"vs",batchMTLD[,2],".txt"))
fileL<-c(fileL, paste0(batchTMLD[,1],"vs",batchTMLD[,2],".txt"))
fileL<-c(fileL, paste0(batchffLD[,1],"vs",batchffLD[,2],".txt"))
fileL<-c(fileL, paste0(batchMTSD[,1],"vs",batchMTSD[,2],".txt"))
fileL<-c(fileL, paste0(batchTMSD[,1],"vs",batchTMSD[,2],".txt"))
fileL<-c(fileL, paste0(batchffSD[,1],"vs",batchffSD[,2],".txt"))


sigT<-c("sample 1","sample 2","DE (q<0.05)","1>2","2>1")
for(file in fileL)
{
  res<-read.table(paste0("DE-27816-mar2019/",file),sep="\t",header=TRUE)
  sigRes <- c(unlist(strsplit(gsub(".txt","",file),split="vs") ), getSig(res),getSig(res,direction="up"),getSig(res,direction="down"))
  sigT<-rbind(sigT,sigRes)
}
T<-as.data.frame(sigT[-1,], row.names=FALSE)
names(T)<-sigT[1,]
write.table(T,file="DE-27816-mar2019/Genes27816.DEsummary.txt",sep="\t")

library(gplots)
pdf("DE-27816-mar2019/checkDE.pdf")
textplot(T)
mtext("DE analysis result summary")


##### Check DE results with different cutoffs
id<-rownames(BLD7)

fc <- seq(1,3,by=0.1)
getSig<-function(res,fc.threshold=0){no<-nrow(res[res$padj<0.05 & !is.na(res$padj) & abs(res$log2FoldChange)>=fc.threshold,]);return(no)}
AB<-data.frame(
  foldChangeCutoff =fc,
  B.mtLD7=sapply(log2(fc),function(x)getSig(B.mtLD7,x)),
  B.mtLD9=sapply(log2(fc),function(x)getSig(B.mtLD9,x)),
  B.mtSD7=sapply(log2(fc),function(x)getSig(B.mtSD7,x)),
  B.mtSD5=sapply(log2(fc),function(x)getSig(B.mtSD5,x)),
  B.tmLD7=sapply(log2(fc),function(x)getSig(B.tmLD7,x)),
  B.tmLD9=sapply(log2(fc),function(x)getSig(B.tmLD9,x)),
  B.tmSD7=sapply(log2(fc),function(x)getSig(B.tmSD7,x)),
  B.tmSD5=sapply(log2(fc),function(x)getSig(B.tmSD5,x)),
  BLD7=sapply(log2(fc),function(x)getSig(BLD7,x)),
  BLD9=sapply(log2(fc),function(x)getSig(BLD9,x)),
  BSD7=sapply(log2(fc),function(x)getSig(BSD7,x)),
  BSD5=sapply(log2(fc),function(x)getSig(BSD5,x)),
  ALD7=sapply(log2(fc),function(x)getSig(ALD7,x)),
  ALD9=sapply(log2(fc),function(x)getSig(ALD9,x)),
  ASD7=sapply(log2(fc),function(x)getSig(ASD7,x)),
  ASD5=sapply(log2(fc),function(x)getSig(ASD5,x))
)
library(RColorBrewer)
n<-ncol(AB)-1
col=brewer.pal(n,"Dark2")
plot(AB[,2]~foldChangeCutoff,data=AB, type="l", col=col[1], ylim=c(1000,7000), ylab="DE genes", main="DE analysis with different cutoffs")
for(i in 1:n)
{
  lines(AB[,i+1]~foldChangeCutoff, data=AB, type="l", col=col[i])
}
#legend("topright",legend=names(AB)[2:7],col=col,pch=21)
text(1.1, AB[1,2:(n+1)]+100, names(AB)[2:(n+1)],)

par(mfrow=c(2,2))

#Long Day 7am
# compare log2FoldChange distribution
# No apparant difference between A and B
plot(density(na.omit(ALD7[id,"log2FoldChange"])), main="log2 fold change, 10 dpa")
lines(density(na.omit(BLD7[id,"log2FoldChange"])),col="purple")
lines(density(na.omit(B.tmLD7[id,"log2FoldChange"])),col="blue")
lines(density(na.omit(B.mtLD7[id,"log2FoldChange"])),col="green")
abline(v=0,col="light grey")
legend("topright",col=c("black","purple","blue","green"),legend=c("A","B","B - TxM","B - MxT"),pch="-")
# compare standard errors
# Higher standard error in B than A
plot( ALD7[id,"lfcSE"], BLD7[id,"lfcSE"],pch=".",main="Standard error, Long day 7am",xlab="A",ylab="B")
lines(c(0,1),c(0,1),col="blue")
plot( ALD7[id,"lfcSE"], B.mtLD7[id,"lfcSE"],pch=".",main="Standard error, Long day 7am",xlab="A",ylab="B - MxT")
lines(c(0,1),c(0,1),col="blue")
plot( ALD7[id,"lfcSE"], B.tmLD7[id,"lfcSE"],pch=".",main="Standard error, Long day 7am",xlab="A",ylab="B - TxM")
lines(c(0,1),c(0,1),col="blue")

#Long day 9pm 
# compare log2FoldChange distribution
# No apparant difference between A and B
plot(density(na.omit(ALD9[id,"log2FoldChange"])), main="log2 fold change, 20 dpa")
lines(density(na.omit(BLD9[id,"log2FoldChange"])),col="purple")
lines(density(na.omit(B.tmLD9[id,"log2FoldChange"])),col="blue")
lines(density(na.omit(B.mtLD9[id,"log2FoldChange"])),col="green")
abline(v=0,col="light grey")
legend("topright",col=c("black","purple","blue","green"),legend=c("A","B","B - TxM","B - MxT"),pch="-")
# compare standard errors
# Higher standard error in B than A
plot( ALD9[id,"lfcSE"], BLD9[id,"lfcSE"],pch=".",main="Standard error, Long day 9pm",xlab="A",ylab="B")
lines(c(0,1),c(0,1),col="blue")
plot( ALD9[id,"lfcSE"], B.mtLD9[id,"lfcSE"],pch=".",main="Standard error, Long day 9pm",xlab="A",ylab="B - MxT")
lines(c(0,1),c(0,1),col="blue")
plot( ALD9[id,"lfcSE"], B.tmLD9[id,"lfcSE"],pch=".",main="Standard error, Long day 9pm",xlab="A",ylab="B - TxM")
lines(c(0,1),c(0,1),col="blue")

#Short Day 5pm
# compare log2FoldChange distribution
# No apparant difference between A and B
plot(density(na.omit(ASD5[id,"log2FoldChange"])), main="log2 fold change, 10 dpa")
lines(density(na.omit(BSD5[id,"log2FoldChange"])),col="purple")
lines(density(na.omit(B.tmSD5[id,"log2FoldChange"])),col="blue")
lines(density(na.omit(B.mtSD5[id,"log2FoldChange"])),col="green")
abline(v=0,col="light grey")
legend("topright",col=c("black","purple","blue","green"),legend=c("A","B","B - TxM","B - MxT"),pch="-")
# compare standard errors
# Higher standard error in B than A
plot( ASD5[id,"lfcSE"], BSD5[id,"lfcSE"],pch=".",main="Standard error, Short day 5pm",xlab="A",ylab="B")
lines(c(0,1),c(0,1),col="blue")
plot( ASD5[id,"lfcSE"], B.mtSD5[id,"lfcSE"],pch=".",main="Standard error, Short day 5pm",xlab="A",ylab="B - MxT")
lines(c(0,1),c(0,1),col="blue")
plot( ASD5[id,"lfcSE"], B.tmSD5[id,"lfcSE"],pch=".",main="Standard error, Short day 5pm",xlab="A",ylab="B - TxM")
lines(c(0,1),c(0,1),col="blue")

#Short day 7am 
# compare log2FoldChange distribution
# No apparant difference between A and B
plot(density(na.omit(ASD7[id,"log2FoldChange"])), main="log2 fold change, 20 dpa")
lines(density(na.omit(BSD7[id,"log2FoldChange"])),col="purple")
lines(density(na.omit(B.tmSD7[id,"log2FoldChange"])),col="blue")
lines(density(na.omit(B.mtSD7[id,"log2FoldChange"])),col="green")
abline(v=0,col="light grey")
legend("topright",col=c("black","purple","blue","green"),legend=c("A","B","B - TxM","B - MxT"),pch="-")
# compare standard errors
# Higher standard error in B than A
plot( ASD7[id,"lfcSE"], BSD7[id,"lfcSE"],pch=".",main="Standard error, Short day 7am",xlab="A",ylab="B")
lines(c(0,1),c(0,1),col="blue")
plot( ASD7[id,"lfcSE"], B.mtSD7[id,"lfcSE"],pch=".",main="Standard error, Short day 7am",xlab="A",ylab="B - MxT")
lines(c(0,1),c(0,1),col="blue")
plot( ASD7[id,"lfcSE"], B.tmSD7[id,"lfcSE"],pch=".",main="Standard error, Short day 7am",xlab="A",ylab="B - TxM")
lines(c(0,1),c(0,1),col="blue")

### Check standard error of A and B
# volcano plot
plotVolvano<-function(res, title)
{
  plot(res[id,"log2FoldChange"], -log2(res[id,"padj"]), main=title, xlab="log2FoldChange", ylab="-log2padj",pch=".",ylim=c(0,200))
  abline(h=-log2(0.05))
}
plotVolvano(ALD7, "A, LD7")
plotVolvano(BLD7, "B, F1 LD7")
plotVolvano(B.mtLD7, "B, MxT LD7")
plotVolvano(B.tmLD7, "B, TxM LD7")
plotVolvano(ALD9, "A, LD9")
plotVolvano(BLD9, "B, F1 LD9")
plotVolvano(B.mtLD9, "B, MxT LD9")
plotVolvano(B.tmLD9, "B, TxM LD9")
plotVolvano(ASD7, "A, SD7")
plotVolvano(BSD7, "B, F1 SD7")
plotVolvano(B.mtSD7, "B, MxT SD7")
plotVolvano(B.tmSD7, "B, TxM SD7")
plotVolvano(ASD5, "A, SD5")
plotVolvano(BSD5, "B, F1 SD5")
plotVolvano(B.mtSD5, "B, MxT SD5")
plotVolvano(B.tmSD5, "B, TxM SD5")


# compare log2 Fold change
plot( ALD7[id,"log2FoldChange"], BLD7[id,"log2FoldChange"],pch=".",main="log2 Fold Change, LD7",xlab="A",ylab="B - F1")
lines(c(-6,6),c(-6,6),col="blue")
plot( ALD7[id,"log2FoldChange"], B.mtLD7[id,"log2FoldChange"],pch=".",main="log2 Fold Change, LD7",xlab="A",ylab="B - MxT")
lines(c(-6,6),c(-6,6),col="blue")
plot( ALD7[id,"log2FoldChange"], B.tmLD7[id,"log2FoldChange"],pch=".",main="log2 Fold Change, LD7",xlab="A",ylab="B - TxM")
lines(c(-6,6),c(-6,6),col="blue")
plot( B.mtLD7[id,"log2FoldChange"], B.tmLD7[id,"log2FoldChange"],pch=".",main="log2 Fold Change, LD7",xlab="B - MxT",ylab="B - TxM")
lines(c(-6,6),c(-6,6),col="blue")


plot( ALD9[id,"log2FoldChange"], BLD9[id,"log2FoldChange"],pch=".",main="log2 Fold Change, LD9",xlab="A",ylab="B - F1")
lines(c(-6,6),c(-6,6),col="blue")
plot( ALD9[id,"log2FoldChange"], B.mtLD9[id,"log2FoldChange"],pch=".",main="log2 Fold Change, LD9",xlab="A",ylab="B - MxT")
lines(c(-6,6),c(-6,6),col="blue")
plot( ALD9[id,"log2FoldChange"], B.tmLD9[id,"log2FoldChange"],pch=".",main="log2 Fold Change, LD9",xlab="A",ylab="B - TxM")
lines(c(-6,6),c(-6,6),col="blue")
plot( B.mtLD9[id,"log2FoldChange"], B.tmLD9[id,"log2FoldChange"],pch=".",main="log2 Fold Change, LD9",xlab="B - MxT",ylab="B - TxM")
lines(c(-6,6),c(-6,6),col="blue")

plot( ASD7[id,"log2FoldChange"], BSD7[id,"log2FoldChange"],pch=".",main="log2 Fold Change, SD7",xlab="A",ylab="B - F1")
lines(c(-6,6),c(-6,6),col="blue")
plot( ASD7[id,"log2FoldChange"], B.mtSD7[id,"log2FoldChange"],pch=".",main="log2 Fold Change, SD7",xlab="A",ylab="B - MxT")
lines(c(-6,6),c(-6,6),col="blue")
plot( ASD7[id,"log2FoldChange"], B.tmSD7[id,"log2FoldChange"],pch=".",main="log2 Fold Change, SD7",xlab="A",ylab="B - TxM")
lines(c(-6,6),c(-6,6),col="blue")
plot( B.mtSD7[id,"log2FoldChange"], B.tmSD7[id,"log2FoldChange"],pch=".",main="log2 Fold Change,SD7",xlab="B - MxT",ylab="B - TxM")
lines(c(-6,6),c(-6,6),col="blue")


plot( ASD5[id,"log2FoldChange"], BSD5[id,"log2FoldChange"],pch=".",main="log2 Fold Change, SD5",xlab="A",ylab="B - F1")
lines(c(-6,6),c(-6,6),col="blue")
plot( ASD5[id,"log2FoldChange"], B.mtSD5[id,"log2FoldChange"],pch=".",main="log2 Fold Change, SD5",xlab="A",ylab="B - MxT")
lines(c(-6,6),c(-6,6),col="blue")
plot( ASD5[id,"log2FoldChange"], B.tmSD5[id,"log2FoldChange"],pch=".",main="log2 Fold Change, SD5",xlab="A",ylab="B - TxM")
lines(c(-6,6),c(-6,6),col="blue")
plot( B.mtSD5[id,"log2FoldChange"], B.tmSD5[id,"log2FoldChange"],pch=".",main="log2 Fold Change, SD5",xlab="B - MxT",ylab="B - TxM")
lines(c(-6,6),c(-6,6),col="blue")
dev.off()


#####################################################
############### Cis/Trans Analysis ##################
#####################################################


### comparing A-B= 0 is tricky, both are log2FoldChange and its standard error lfcse
# maybe I can compare with t test
#### T test from means and standard errors ####
# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# se1, se2: the sample standard errors
# se1 <- s1/sqrt(n)
# m0: the null value for the difference in means to be tested for. Default is 0.
# equal.variance: whether or not to assume equal variance. Default is FALSE.
t.test2 <- function(m1,m2,se1,se2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE )
  {
    # se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    se <- sqrt( (se1^2) + (se2^2) )
    # welch-satterthwaite df
    # df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    df <- ( (se1^2 + se2^2)^2 )/( (se1^2)^2/(n1-1) + (se2^2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    # se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
    df <- n1+n2-2
  }
  t <- (m1-m2-m0)/se
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat)
}

x1 = rnorm(3)
x2 = rnorm(3)
# you'll find this output agrees with that of t.test when you input x1,x2
t.test(x1,x2)
t.test2( mean(x1),  mean(x2), sd(x1)/sqrt(3), sd(x2)/sqrt(3), 3,3)

### Make a categorization table for each F1 condition: TxM10, MxT10, TxM20, MxT20,
criteria<- as.data.frame(rbind(c("A!=0;B!=0;A=B", "1.Cis only"),
                               c("A!=0;B=0;A!=B", "2.Trans only"),
                               c("A!=0;B!=0;A!=B", "Cis+Trans"),
                               c("A=0;B!=0;A!=B", "5.Compensatory"),
                               c("A=0;B=0;A=B", "6.Conserved") ))
names(criteria) <- c("class","category")

# A.res <-A10
# B.res <-B10

classCisTrans<-function(A.res, B.res, A.n, B.n, log2fc.threshold=0,plotTitle=NULL)
{
  # A = log2(TX2094/Maxxa), cis + trans
  A <- A.res[,c("log2FoldChange", "lfcSE", "padj")]
  names(A) <- c("A", "A.SE", "A.padj")
  # A = log2(F1_t/F1_m), cis
  B <- B.res[,c("log2FoldChange", "lfcSE", "padj")]
  names(B) <- c("B", "B.SE", "B.padj")
  
  A<-A[rownames(B),]
  table <- cbind(A,B)
  table$AminusB <- table$A - table$B
  table$AminusB.pvalue <- apply(table,1,function(x) t.test2(m1=x[1],m2=x[4],se1=x[2], se2=x[5], n1=A.n, n2=B.n)["p-value"])
  
  table$cisNtrans <- ifelse(abs(table$A)>=log2fc.threshold & table$A.padj<0.05 & !is.na(table$A.padj), "A!=0", "A=0")
  table$cis <- ifelse(abs(table$B)>=log2fc.threshold & table$B.padj<0.05 & !is.na(table$B.padj), "B!=0", "B=0")
  table$trans <- ifelse(abs(table$AminusB)>=log2fc.threshold & table$AminusB.pvalue<0.05, "A!=B", "A=B")
  table$class <- paste(table$cisNtrans,table$cis,table$trans,sep=";")
  table$category <- as.character(criteria$category[ match(table$class,criteria$class)])
  table$category[is.na(table$category)] <- "7.Ambiguous"
  table$category[ table$category=="Cis+Trans" & table$B*table$AminusB >0 ] <- "3.Cis+Trans: enhancing"
  table$category[ table$category=="Cis+Trans" & table$B*table$AminusB <0 ] <- "4.Cis+Trans: compensating"
  
  colors <- c("red","blue","purple","brown","green","black","grey")
  if(!is.null(plotTitle)){
    p<- ggplot( table, aes(x=A, y=B, color=category)) + geom_point(alpha=0.8) +  xlab("Cis + Trans") + ylab("Cis") + ggtitle(plotTitle) + scale_color_manual(values=colors) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +  geom_abline(intercept = 0, slope = 1) + theme_bw()
    # p2<-ggplot( table, aes( factor(category), fill=category)) + geom_bar(stat="count") + scale_fill_manual(values=colors)
    print(p)
  }
  return(table)
}
resLD7 <- classCisTrans(A.res = ALD7, B.res = BLD7, A.n = 3, B.n=6, log2fc.threshold=0)
res.mtLD7 <- classCisTrans(A.res = ALD7, B.res = B.mtLD7, A.n = 3, B.n=3, log2fc.threshold=0)
res.tmLD7 <- classCisTrans(A.res = ALD7, B.res = B.tmLD7, A.n = 3, B.n=3, log2fc.threshold=0)
resLD9 <- classCisTrans(A.res = ALD9, B.res = BLD9, A.n = 3, B.n=6, log2fc.threshold=0)
res.mtLD9 <- classCisTrans(A.res = ALD9, B.res = B.mtLD9, A.n = 3, B.n=3, log2fc.threshold=0)
res.tmLD9 <- classCisTrans(A.res = ALD9, B.res = B.tmLD9, A.n = 3, B.n=3, log2fc.threshold=0)
resSD7 <- classCisTrans(A.res = ASD7, B.res = BSD7, A.n = 3, B.n=6, log2fc.threshold=0)
res.mtSD7 <- classCisTrans(A.res = ASD7, B.res = B.mtSD7, A.n = 3, B.n=3, log2fc.threshold=0)
res.tmSD7 <- classCisTrans(A.res = ASD7, B.res = B.tmSD7, A.n = 3, B.n=3, log2fc.threshold=0)
resSD5 <- classCisTrans(A.res = ASD5, B.res = BLD9, A.n = 3, B.n=6, log2fc.threshold=0)
res.mtSD5 <- classCisTrans(A.res = ASD5, B.res = B.mtSD5, A.n = 3, B.n=3, log2fc.threshold=0)
res.tmSD5 <- classCisTrans(A.res = ASD5, B.res = B.tmSD5, A.n = 3, B.n=3, log2fc.threshold=0)

plotCisTrans<-function(table, plotTitle="")
{
  colors <- c("red","blue","purple","brown","green","black","grey")
  p<- ggplot( table, aes(x=A, y=B, color=category)) + geom_point(alpha=0.8) +  xlab("Cis + Trans") + ylab("Cis") + ggtitle(plotTitle) + scale_color_manual(values=colors) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +  geom_abline(intercept = 0, slope = 1) + theme_bw()
  # p2<-ggplot( table, aes( factor(category), fill=category)) + geom_bar(stat="count") + scale_fill_manual(values=colors)
  print(p)
}


sumT<-as.data.frame(matrix(unlist(lapply(list(res.mtLD7,res.tmLD7,res.mtLD9,res.tmLD9,resLD7,resLD9), function(x)table(x$category)) ),ncol=6))
rownames(sumT)<-names(table(resLD9$category))
names(sumT)<-c("LD7 MxT","LD7 TxM","LD9 MxT","LD9 TxM","LD7 F1","LD9 F1")
print(sumT)

sumT2<-as.data.frame(matrix(unlist(lapply(list(res.mtSD7,res.tmSD7,res.mtSD5,res.tmSD5,resSD7,resSD5), function(x)table(x$category)) ),ncol=6))
rownames(sumT2)<-names(table(resLD9$category))
names(sumT2)<-c("SD7 MxT","SD7 TxM","SD5 MxT","SD5 TxM","SD7 F1","SD5 F1")
print(sumT2)

pdf("plotCistransLongDay.pdf")
textplot(sumT,cex=0.6)
plotCisTrans(res.mtLD7, "MxT LD7")
plotCisTrans(res.mtLD9, "MxT LD9")
plotCisTrans(res.tmLD7, "TxM LD7")
plotCisTrans(res.tmLD9, "TxM LD9")
plotCisTrans(resLD7, "F1 LD7")
plotCisTrans(resLD9, "F1 LD9")
dev.off()

pdf("plotCistransShortDay.pdf")
textplot(sumT2,cex=0.6)
plotCisTrans(res.mtSD7, "MxT SD7")
plotCisTrans(res.mtSD5, "MxT SD5")
plotCisTrans(res.tmSD7, "TxM SD7")
plotCisTrans(res.tmSD5, "TxM SD5")
plotCisTrans(resSD7, "F1 SD7")
plotCisTrans(resSD5, "F1 SD5")
dev.off()


#save(A10, B10, B.mt10, B.tm10, A20, B20, B.mt20, B.tm20, res10, res.mt10, res.tm10, res20, res.mt20, res.tm20, sumT, file="cistrans.rdata")

# plot percentage 
library(reshape2);
library(plyr)
library(ggplot2)
colors <- c("red","blue","purple","brown","green","black","grey")
# plot 7 categories
T=sumT[,1:4]; 
T$category=rownames(T)
dat=melt(T,id.var="category")
p1<-ggplot(data=dat, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=colors)+
  theme_minimal()
df <- ddply(dat, "variable",  transform, label_ypos=rev(cumsum(rev(value))))
df <- ddply(df, "variable",  transform, t=round(value/sum(value)*100,1))
df$value2=paste0(df$value," (",df$t, "%)")
df$value2[grep("^[1-4]",df$category)]=""
p2<-ggplot(data=df, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=value2), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_manual(values=colors)+
  theme_minimal()
# plot 4 categories
T=sumT[1:4,1:4]; 
T$category=rownames(T)
dat=melt(T,id.var="category")
df <- ddply(dat, "variable",  transform, label_ypos=rev(cumsum(rev(value))))
df <- ddply(df, "variable",  transform, t=round(value/sum(value)*100,1))
df$value2=paste0(df$value," (",df$t, "%)")
p3<-ggplot(data=df, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=colors)+
  theme_minimal()
p4<-ggplot(data=df, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=value2), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_manual(values=colors)+
  theme_minimal()
pdf("plotCistransLongDay.barsummary.pdf")
p1;p2;p3;p4
dev.off()

T=sumT2[1:4,1:4]; 
T$category=rownames(T)
dat=melt(T,id.var="category")
df <- ddply(dat, "variable",  transform, label_ypos=rev(cumsum(rev(value))))
df <- ddply(df, "variable",  transform, t=round(value/sum(value)*100,1))
df$value2=paste0(df$value," (",df$t, "%)")
p3<-ggplot(data=df, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=colors)+
  theme_minimal()
p4<-ggplot(data=df, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=value2), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_manual(values=colors)+
  theme_minimal()
pdf("plotCistransShortDay.barsummary.pdf")
p1;p2;p3;p4
dev.off()

write.table(sumT, file ="cistransLongDay.summary.txt",sep="\t")
write.table(sumT2, file ="cistransShortDay.summary.txt",sep="\t")
write.table(resLD7, file ="cistrans.F1_LD7.txt",sep="\t")
write.table(res.mtLD7, file ="cistrans.MxT_LD7.txt",sep="\t")
write.table(res.tmLD7, file ="cistrans.TxM_LD7.txt",sep="\t")
write.table(resLD9, file ="cistrans.F1_LD9.txt",sep="\t")
write.table(res.mtLD9, file ="cistrans.MxT_LD9.txt",sep="\t")
write.table(res.tmLD9, file ="cistrans.TxM_LD9.txt",sep="\t")
write.table(resSD5, file ="cistrans.F1_SD5.txt",sep="\t")
write.table(res.mtSD5, file ="cistrans.MxT_SD5.txt",sep="\t")
write.table(res.tmSD5, file ="cistrans.TxM_SD5.txt",sep="\t")
write.table(resSD7, file ="cistrans.F1_SD7.txt",sep="\t")
write.table(res.mtSD7, file ="cistrans.MxT_SD7.txt",sep="\t")
write.table(res.tmSD7, file ="cistrans.TxM_SD7.txt",sep="\t")

###################################################################
## Change fold change cutoff and inspect cis trans summary table ##
###################################################################

library(RColorBrewer)
col=brewer.pal(7,"Dark2")
n<-length(id)
tests<-rbind(c("ALD7","B.mtLD7","MxTLD7"),c("ALD7","B.tmLD7","TxMLD7"),c("ALD9","B.mtLD9","MxTLD9"),c("ALD9","B.tmLD9","TxMLD9"), c("ALD7","BLD7","F1_LD7"),c("ALD9","BLD9","F1_LD9"))
resultT<-list()
for(i in 1:4){
  res<-sapply(log2(fc),function(x){tt<-classCisTrans(A.res = get(tests[i,1]), B.res = get(tests[i,2]), 3,3,log2fc.threshold=x ); return(table(tt$category))} )
  colnames(res)<-fc
  resultT[[tests[i,3]]]<-res
}
resultT
for(i in 5:6){
  res<-sapply(log2(fc),function(x){tt<-classCisTrans(A.res = get(tests[i,1]), B.res = get(tests[i,2]), 3,6,log2fc.threshold=x ); return(table(tt$category))} )
  colnames(res)<-fc
  resultT[[tests[i,3]]]<-res
}

library(RColorBrewer)
col=brewer.pal(7,"Dark2")
n<-length(id)
tests2<-rbind(c("ASD7","B.mtSD7","MxTSD7"),c("ASD7","B.tmSD7","TxMSD7"),c("ASD5","B.mtSD5","MxTSD5"),c("ASD5","B.tmSD5","TxMSD5"), c("ASD7","BSD7","F1_SD7"),c("ASD5","BSD5","F1_SD5"))
resultT2<-list()
for(i in 1:4){
  res<-sapply(log2(fc),function(x){tt<-classCisTrans(A.res = get(tests2[i,1]), B.res = get(tests2[i,2]), 3,3,log2fc.threshold=x ); return(table(tt$category))} )
  colnames(res)<-fc
  resultT2[[tests2[i,3]]]<-res
}
resultT2
for(i in 5:6){
  res<-sapply(log2(fc),function(x){tt<-classCisTrans(A.res = get(tests2[i,1]), B.res = get(tests2[i,2]), 3,6,log2fc.threshold=x ); return(table(tt$category))} )
  colnames(res)<-fc
  resultT2[[tests2[i,3]]]<-res
}

colors <- c("red","blue","purple","brown","green","black","grey")
pdf("checkCistranLongDay.by.cutoff.pdf")
noCut<-cbind(resultT[[1]][,1],resultT[[2]][,1],resultT[[3]][,1],resultT[[4]][,1],resultT[[5]][,1],resultT[[6]][,1])
colnames(noCut)<-names(resultT)
barplot(noCut,col=col, xlab=names(resultT), xlim=c(0,10),legend=rownames(noCut),las=2, main="fold change cutoff: 1")
for(i in 1:6){
  res<-resultT[[i]]
  barplot(res,xlab=fc,col=col, main=tests[i,3],las=2)
}
for(i in 1:6){
  res<-resultT[[i]][1:5,]
  barplot(res,xlab=fc,col=col[1:5], main=tests[i,3],las=2)
}
dev.off()

colors <- c("red","blue","purple","brown","green","black","grey")
pdf("checkCistranShortDay.by.cutoff.pdf")
noCut<-cbind(resultT2[[1]][,1],resultT2[[2]][,1],resultT2[[3]][,1],resultT2[[4]][,1],resultT2[[5]][,1],resultT2[[6]][,1])
colnames(noCut)<-names(resultT2)
barplot(noCut,col=col, xlab=names(resultT2), xlim=c(0,10),legend=rownames(noCut),las=2, main="fold change cutoff: 1")
for(i in 1:6){
  res<-resultT2[[i]]
  barplot(res,xlab=fc,col=col, main=tests[i,3],las=2)
}
for(i in 1:6){
  res<-resultT2[[i]][1:5,]
  barplot(res,xlab=fc,col=col[1:5], main=tests[i,3],las=2)
}
dev.off()

