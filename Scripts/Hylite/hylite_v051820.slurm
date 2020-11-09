#!/bin/bash
#SBATCH --constraint=AVX2
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 36   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH --mem=500G 
#SBATCH --time=30-02:30:00
#SBATCH --output=job.%j.out.txt
#SBATCH --error=job.%j.err.txt
#SBATCH -J "CisTransHylite"  
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=hugj2006@iastate.edu 
module load python/3.6.3-u4oaxsb samtools/1.9-k6deoga py-scipy/1.1.0-py3-pzig4lr
#python3 /work/LAS/jfw-lab/HyLiTE-2.0.1-py3-none-any.whl -v -S -f hylite/cistrans_protocol.txt -r TM1.nuclNmt.transcript.fasta -n hyliteResults 2>hyliteResults.log
samtools mpileup -B -A -Q 0 -d 1000000 -f TM1.nuclNmt.transcript.fasta -b /work/LAS/jfw-lab/erdostal/CisTrans/hyliteResults/hyliteResults.txt > /work/LAS/jfw-lab/erdostal/CisTrans/hyliteResults/hyliteResults.pileup
python3 /work/LAS/jfw-lab/HyLiTE-2.0.1-py3-none-any.whl -v -f hylite/cistrans_protocol.txt -p hyliteResults/hyliteResults.pileup -r TM1.nuclNmt.transcript.fasta -n hyliteResults 2>hyliteResults.2.log

