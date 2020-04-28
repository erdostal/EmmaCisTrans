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
htseq-count -f bam --idattr Parent --stranded=no sortbam/SDPF-Maxxa-S1-202_FCH35H3BBXX-HKRDCOTwkbEACFRAAPEI-202_L3_.sort.bam Tx-JGI_G.hirsutum_v1.1.gene_exons.gff3.gz > counts/SDPFtest

## Check quality of count files
for m in $(ls *counts.txt); do echo $m; cut $m -f2 |awk '{total = total + $1}END{print total}'; done




