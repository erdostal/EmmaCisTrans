#####################################################
############### SNP categorization ##################
#####################################################

x<-read.table(grep("hyliteResults.snp.txt",list.files(),value=TRUE),header=TRUE,sep="\t")
dim(x) # total snps 1062283
#Creates a table of occurrences between F1 x Maxxa, Maxxa x TX2094, and F1 x TX2094
type <- as.data.frame(ftable(xtabs(~F1+Maxxa+TX2094, data=x)))
#Remove instances where nothing matched
type<-type[type$Freq>0,]
dim(type)  #45 categories
type$note <- ''
##if value is TRUE in the following lines, assign the corresponding label
#if all columns have -1 mark it as poor coverage
type$note[type$F1==-1 | type$Maxxa ==-1 | type$TX2094 == -1] <- "poor coverage"
type$note[ type$Maxxa == type$TX2094 & type$Maxxa==1 & type$F1==1 ] <- "common"
type$note[ type$Maxxa == type$TX2094 & type$Maxxa==1 & type$F1=="1,0" ] <- "common"
type$note[ type$Maxxa == type$TX2094 & type$Maxxa==1 & (type$F1==0 |type$F1=="0,0" ) ] <- "parents specific"
type$note[ type$Maxxa==1 & type$TX2094 ==0 & (type$F1==0 |type$F1=="0,0") ] <- "Maxxa unique1"
type$note[ type$Maxxa==0 & type$TX2094 ==1 & (type$F1==0 |type$F1=="0,0") ] <- "TX2094 unique1"
type$note[ type$TX2094==0 & type$Maxxa==0 ] <- "F1 unique1"
type$note[ type$Maxxa ==1 & type$TX2094==0 & (type$F1=="1,0" )] <- "allele diagnostic"
type$note[ type$Maxxa ==0 & type$TX2094==1 & (type$F1=="1,0" )] <- "allele diagnostic"
type$note[ type$Maxxa ==1 & type$TX2094==0 & (type$F1==1)] <- "allele diagnostic: maxxa allele in F1"
type$note[ type$Maxxa ==0 & type$TX2094==1 & (type$F1==1)] <- "allele diagnostic: tx2094 allele in F1"

bg<-aggregate(type$Freq,list(type$note),sum)
total<-sum(type$Freq)#1062283
bg$percentage <- bg$x/total
bg
Group.1      x   percentage
#1                      allele diagnostic 166364 0.1566098676
#2  allele diagnostic: maxxa allele in F1    370 0.0003483064
#3 allele diagnostic: tx2094 allele in F1   1374 0.0012934406
#4                                 common  96845 0.0911668548
#5                             F1 unique1 636643 0.5993158132
#6                          Maxxa unique1  11227 0.0105687467
#7                       parents specific    396 0.0003727820
#8                          poor coverage 142838 0.1344632268
#9                         TX2094 unique1   6226 0.0058609617
type<-type[order(type$note),]
write.table(type, "SNPtypes.txt",row.names=FALSE,sep="\t")

# what genes contain disgnostic SNPs
x$note <- ''
x$note[x$F1==-1 | x$Maxxa ==-1 | x$TX2094 == -1] <- "poor coverage"
x$note[ x$Maxxa == x$TX2094 & x$Maxxa==1 & x$F1==1 ] <- "common"
x$note[ x$Maxxa == x$TX2094 & x$Maxxa==1 & x$F1=="1,0" ] <- "common"
x$note[ x$Maxxa == x$TX2094 & x$Maxxa==1 & (x$F1==0 |x$F1=="0,0" ) ] <- "parents specific"
x$note[ x$Maxxa==1 & x$TX2094 ==0 & (x$F1==0 |x$F1=="0,0") ] <- "Maxxa unique1"
x$note[ x$Maxxa==0 & x$TX2094 ==1 & (x$F1==0 |x$F1=="0,0") ] <- "TX2094 unique1"
x$note[ x$Maxxa == x$TX2094 & x$Maxxa==0 ] <- "F1 unique1"
x$note[ x$Maxxa ==1 & x$TX2094==0 & (x$F1=="1,0" | x$F1==1)] <- "allele diagnostic: Maxxa"
x$note[ x$Maxxa ==0 & x$TX2094==1 & (x$F1=="1,0" | x$F1==1)] <- "allele diagnostic: tx2094"

length(unique(x$GENE))
# No. of genes contain SNPs 61696
select <- grep("diagnostic",x$note)
length(unique(x$GENE[select]))
# No. of genes contain diagnostic SNPs 39872

#####################################################
################# SNPs on genes #####################
#####################################################

y<-read.table(grep("snp.summary.txt",list.files(),value=TRUE),header=TRUE,sep="\t")
nrow(y) # total No. of genes 75377
y$snps <- rowSums(y[,-1])
table(y$snp>0) #61696 contain SNPs
table(y$F1.TX2094>0 | y$F1.Maxxa>0)  # 39872 contains SNPs diagonostic to alleles in F1
genes <- as.character( y$GENE[y$F1.TX2094>0 | y$F1.Maxxa>0] )
length(genes) #39872

#####################################################
######### Allelic specific expression ###############
#####################################################

fileL<-grep(".read.summary.txt", list.files(), value=TRUE)
# keep allele-specific counts to 4 data frames
alleleM  <- data.frame(gene=y$GENE)    # Maxxa
alleleMn <- data.frame(gene=y$GENE)    # Maxxa + Maxxa.N; N refers to F1 unique SNPs
alleleT  <- data.frame(gene=y$GENE)    # TX2094
alleleTn <- data.frame(gene=y$GENE)    # TX2094
readT <- 0
for(file in fileL){
  ac <- read.table(file,header=TRUE,sep="\t")
  tag<-gsub("-",".",gsub(".*F1[.]|.read.summary.txt","",file))
  alleleM[,tag]   <- ac$Maxxa
  alleleMn[,tag]  <- ac$Maxxa + ac$Maxxa.N
  alleleT[,tag]   <- ac$TX2094
  alleleTn[,tag]  <- ac$TX2094 + ac$TX2094.N
  readT=readT+sum(ac[,-1])
}

readT #432427971
sum(alleleM[,-1]) # No. of reads containing only Maxxa diagnostic SNPs 25813229
sum(alleleMn[,-1]) # No. of reads containing Maxxa and F1 diagnostic SNPs ** 28456815
sum(alleleTn[,-1]) # No. of reads containing TX2094 and F1 diagnostic SNPs ** 30834357
sum(alleleT[,-1]) # No. of reads containing TX2094 diagnostic SNPs 27829776
# aggregate per F1 sample
for(file in fileL){
  ac <- read.table(file,header=TRUE,sep="\t")
  tag<-gsub("-",".",gsub(".*F1[.]|.read.summary.txt","",file))
  if(!exists("F1readSummary")){F1readSummary = c(tag,colSums(ac[,-1]))}else{F1readSummary = rbind(F1readSummary,c(tag,colSums(ac[,-1])))}
}
rownames(F1readSummary) = F1readSummary[,1]
F1readSummary<-as.data.frame(F1readSummary[,-1])
write.table(F1readSummary,file="F1readSummary.txt",sep="\t")

# check genes without allelic counts
check<-alleleMn$gene[rowSums(cbind(alleleMn[,-1],alleleTn[,-1]))>0]
length(check) # No. of genes containing allelic counts 41705
length(check)>length(genes) # bigger than final reported number
length(intersect(check,genes)) # 39872
length(setdiff(check,genes)) # 1833
length(setdiff(genes,check)) # 0
# check results mostly agree with genes containning disgnostic snps, the discrepancy probably came from MASKED snps, which means low coverage SNP site in some accession but probably allowing read assignment in some other accessions
# SO just focus on the common set
genes_diagnostic_with_alleles <- intersect(check,genes)
# save(alleleM, alleleMn, alleleT, alleleTn, genes, check, genes_diagnostic_with_alleles, file="Allelic_read_count.Rdata" )


#####################################################
################ Total Expression ###################
#####################################################

# import read count table:  .expression.txt file
x<-read.table(grep("expression.txt",list.files(),value=TRUE), header=TRUE, sep="\t")
dd<-which(duplicated(x$GENE))
count <- x[-1]
rownames(count) <- x$GENE
names(count) <- gsub("^F1[.]|^TX2094[.]|^Maxxa[.]","",names(count))

#####################################################
#################### Make files  ####################
#####################################################
write.table(count, file = "totalReadCounts.txt", sep="\t")
write.table(alleleMn[alleleMn$gene %in% genes_diagnostic_with_alleles,], row.names=FALSE, file = "MaxxaAlleleReadCounts.txt", sep="\t")
write.table(alleleTn[alleleTn$gene %in% genes_diagnostic_with_alleles,], row.names=FALSE, file = "TX2094AlleleReadCounts.txt", sep="\t")

q()
n

#####################################################
################# File Discription ##################
#####################################################

# 'totalReadCounts.txt' -  counts of total gene expression, genes (row) X 26 RNA-seq libraries (column)
# 'MaxxaAlleleReadCounts.txt' - for 12 F1 libraries (column), allele-specific read count assigned to Maxxa alleles were presented for these genes contain diagnostic SNPs.
# 'TX2094AlleleReadCounts.txt' - same as above but for TX2094 allele.
# 'TM1.nuclNmt.transcript.len.txt' - gene length for all genes
# 'TM1.nuclNmt.transcript.len.diagnostic.txt' - gene length for genes with diagnostic allele-specifc SNPs