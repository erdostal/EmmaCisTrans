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

## load bwa
module load bwa
## build reference index
bwa index TM1.nuclNmt.transcript.fasta

## run slurm array 
sbatch /work/LAS/jfw-lab/erdostal/CisTrans/CisTrans.CEGedit.slurm

## Generate count files from bams
for j in $(ls sortbam/|grep 'bam$'); do samtools idxstat sortbam/$j > counts/${j%%_.sort.bam}.counts.txt;done

## Check quality of count files
for m in $(ls *counts.txt); do echo $m; cut $m -f2 |awk '{total = total + $1}END{print total}'; done

## Count raw reads from fastqs
for j in $(ls RNAseq/*.fq.gz); do paste <(echo $j) <(zcat $j |  awk 'END{ print NR/4 }') >> your.file.txt; done

## Hylite run  
### slurm file = hylite_v05182020
### hylite protocol = cistrans_protocol

## Processing output from hylite
### cistrans.pre.r
### hylitePCA.r

## Based on PCAs Jing has suggested the following
LD7-Maxxa-L2 is very likely LD9-Maxxa, correct it and give it a new name. SD7-F1-5, LD7-F1-4NF and LD7-F1-6 may be removed. But before that, probably first check alleic read counts. So far, PCAs are done at whole gene level, I  think the script will also check alleic level grouping

