# EmmaCisTrans
Working with Cis/Trans data per Jing's instruction


# Set up Reference
## locate reference transcriptome on server
cd ~/jfw-lab/Projects/ChenCottons/Genomes/Gossypium_hirsutum/V2.0

## how many primary transcripts? 
75,376
zcat Ghirsutumv2.1.primaryTrs.fa.gz | grep -c '>'

## rename mt genes
cp mt.transcript.fasta mt.transcript.rename.fasta 
sed -i 's/lcl.*gene=/mt_/g' AD1mitochondria.fasta 
sed -i 's/\].*//g' AD1mitochondria.fasta 

## cat references into one file
cat Ghirsutumv2.1.primaryTrs.fa.gz AD1mitochondria.fasta > TM1.nuclNmt.transcript.fasta

## load bowtie2 Version: 2.2.6 
module load bowtie2
## build bowtie reference
bowtie2-build TM1.nuclNmt.transcript.fasta TM1new

