#!/bin/bash

#SBATCH --nodes=1   
#SBATCH --cpus-per-task=6
#SBATCH --error=/exports/people/andersenlab/mzy668/logs/job.%J.err
#SBATCH --output=/exports/people/andersenlab/mzy668/logs/job.%J.out
#SBATCH -J mz-bwa-align
#SBATCH --workdir=/lscr2/andersenlab/mzy668/celegans/SRNA/4-BWA
#SBATCH --mem=8096

bwa aln -o 0 -n 0 -t 6 ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa ../2-FastQ_q1/N2_q1.fq.gz > N2.sai
bwa samse ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa N2.sai ../2-FastQ_q1/N2_q1.fq.gz > N2.sam
samtools view -bS N2.sam > N2.unsorted.bam
samtools flagstat N2.unsorted.bam
samtools sort -@ 6 -o N2.bam N2.unsorted.bam
samtools index -b N2.bam

bwa aln -o 0 -n 0 -t 6 ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa ../2-FastQ_q1/N2-P_q1.fq.gz > N2P.sai
bwa samse ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa N2P.sai ../2-FastQ_q1/N2-P_q1.fq.gz > N2P.sam
samtools view -bS N2P.sam > N2P.unsorted.bam
samtools flagstat N2P.unsorted.bam
samtools sort -@ 6 -o N2P.bam N2P.unsorted.bam
samtools index -b N2P.bam

bwa aln -o 0 -n 0 -t 6 ../../CB4856_genome/Caenorhabditis_elegans_CB4856_v1.0.fa ../2-FastQ_q1/CB_q1.fq.gz > CB.sai
bwa samse ../../CB4856_genome/Caenorhabditis_elegans_CB4856_v1.0.fa CB.sai ../2-FastQ_q1/CB_q1.fq.gz > CB.sam
samtools view -bS CB.sam > CB.unsorted.bam
samtools flagstat CB.unsorted.bam
samtools sort -@ 6 -o CB.bam CB.unsorted.bam
samtools index -b CB.bam

bwa aln -o 0 -n 0 -t 6 ../../CB4856_genome/Caenorhabditis_elegans_CB4856_v1.0.fa ../2-FastQ_q1/CB-P_q1.fq.gz > CBP.sai
bwa samse ../../CB4856_genome/Caenorhabditis_elegans_CB4856_v1.0.fa CBP.sai ../2-FastQ_q1/CB-P_q1.fq.gz > CBP.sam
samtools view -bS CBP.sam > CBP.unsorted.bam
samtools flagstat CBP.unsorted.bam
samtools sort -@ 6 -o CBP.bam CBP.unsorted.bam
samtools index -b CBP.bam

rm *.sai
rm *.sam
rm *.unsorted.bam