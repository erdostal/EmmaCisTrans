#!/bin/bash
#SBATCH --job-name=array
#SBATCH --constraint=AVX2
#SBATCH --array=1-169%15
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=70
#SBATCH --mem=200G 
#SBATCH --time=4-02:30:00
#SBATCH --output=job.A2.%A_%a.out 
#SBATCH --error=job.A2.%A_%a.err
#SBATCH -J "CisTrans.trimArray"  
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=erdostal@iastate.edu 
module load trimmomatic
module load bwa
module load samtools
DIR=/work/LAS/jfw-lab/erdostal/CisTrans/RNAseq
tDIR=/work/LAS/jfw-lab/erdostal/CisTrans/trimmed
ref=/work/LAS/jfw-lab/erdostal/CisTrans/TM1.nuclNmt.transcript.fasta 
thr=70
file1=$(ls -1 $DIR/*1.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
file2=$(ls -1 $DIR/*2.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
name=$(basename $file1 1.fq.gz) 

trimmomatic PE -threads $thr $file1 $file2 $name.trimmed.1.fq.gz $tDIR/$name.trimmed.2.fq.gz ILLUMINACLIP:/work/LAS/jfw-lab/erdostal/CisTrans/Adapters.fa:2:30:15 LEADING:28 TRAILING:28 SLIDINGWINDOW:8:28 SLIDINGWINDOW:1:10 MINLEN:75 TOPHRED33
bwa mem -M -K 10000000 -R "@RG\tID:$name \tSM:$name \tPL:ILLUMINA" -t $thr $ref $DIR/$file1 $DIR/$file2 > $name.sam &> $name.bwa.timelog

samtools view -@ 2 -u -F 4 $name.sam | samtools sort -m 3G -o $name.sort.bam &> $name.bam.err
samtools index $name.sort.bam
ml 
