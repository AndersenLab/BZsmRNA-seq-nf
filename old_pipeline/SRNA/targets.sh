#!/bin/bash

#SBATCH --nodes=1   
#SBATCH --cpus-per-task=8
#SBATCH --error=/exports/people/andersenlab/mzy668/logs/job.%J.err
#SBATCH --output=/exports/people/andersenlab/mzy668/logs/job.%J.out
#SBATCH -J mz-N2unique
#SBATCH --workdir=/lscr2/andersenlab/mzy668/celegans/SRNA/11-Targets
#SBATCH --mem=8096

bwa aln -o 0 -n 3 -t 8 ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa upstream_QTL_f.fa > US_t3.sai
bwa samse ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa US_t3.sai upstream_QTL_f.fa > US_t3.sam
samtools view -bS US_t3.sam > US_t3.unsorted.bam
samtools flagstatUS_t3.unsorted.bam
samtools sort -@ 8 -o US_t3.bam US_t3.unsorted.bam
samtools index -b US_t3.bam
samtools view US_t3.bam > US_t3.txt 

grep -c "^Chr" US_t3.txt 
grep -c "XA:Z" US_t3.txt 
grep "XA:Z" US_t3.txt  > US_t3f.txt

rm *.sai
rm *.sam
rm *.unsorted.bam

bwa aln -o 0 -n 3 -t 8 ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa Seq_N2-unique.fa > N2-unique_t3.sai
bwa samse ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa N2-unique_t3.sai Seq_N2-unique.fa > N2-unique_t3.sam
samtools view -bS N2-unique_t3.sam > N2-unique_t3.unsorted.bam
samtools flagstat N2-unique_t3.unsorted.bam
samtools sort -@ 8 -o N2-unique_t3.bam N2-unique_t3.unsorted.bam
samtools index -b N2-unique_t3.bam
samtools view N2-unique_t3.bam > N2-unique_t3.txt 

grep -c "^Chr" US_t3.txt 
grep -c "XA:Z" US_t3.txt 
grep "XA:Z" N2-unique_t3.txt  > N2-unique_t3f.txt

rm *.sai
rm *.sam
rm *.unsorted.bam

bwa aln -o 0 -n 3 -t 8 ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa Seq_CB-unique.fa > CB-unique_t3.sai
bwa samse ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa CB-unique_t3.sai Seq_CB-unique.fa > CB-unique_t3.sam
samtools view -bS CB-unique_t3.sam > CB-unique_t3.unsorted.bam
samtools flagstat CB-unique_t3.unsorted.bam
samtools sort -@ 8 -o CB-unique_t3.bam CB-unique_t3.unsorted.bam
samtools index -b CB-unique_t3.bam
samtools view CB-unique_t3.bam > CB-unique_t3.txt 

grep -c "^Chr" US_t3.txt 
grep -c "XA:Z" US_t3.txt 
grep "XA:Z" CB-unique_t3.txt  > CB-unique_t3f.txt

rm *.sai
rm *.sam
rm *.unsorted.bam

rm *.bam
rm *.bam.bai