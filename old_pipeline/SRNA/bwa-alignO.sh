#!/bin/bash

#SBATCH --nodes=1   
#SBATCH --cpus-per-task=10
#SBATCH --error=/exports/people/andersenlab/mzy668/logs/job.%J.err
#SBATCH --output=/exports/people/andersenlab/mzy668/logs/job.%J.out
#SBATCH -J mz-bwa-align
#SBATCH --workdir=/lscr2/andersenlab/mzy668/celegans/SRNA/6-BWAO
#SBATCH --mem=20096

bwa aln -o 0 -n 0 -t 10 ../../CB4856_genome/Caenorhabditis_elegans_CB4856_v1.0.fa ../2-FastQ_q1/N2_q1.fq.gz > N2-CB.sai
bwa samse ../../CB4856_genome/Caenorhabditis_elegans_CB4856_v1.0.fa N2-CB.sai ../2-FastQ_q1/N2_q1.fq.gz > N2-CB.sam
samtools view -bS N2-CB.sam > N2-CB.unsorted.bam
samtools flagstat N2-CB.unsorted.bam
samtools sort -@ 10 -o N2-CB.bam N2-CB.unsorted.bam
samtools index -b N2-CB.bam

bwa aln -o 0 -n 0 -t 10 ../../CB4856_genome/Caenorhabditis_elegans_CB4856_v1.0.fa ../2-FastQ_q1/N2-P_q1.fq.gz > N2P-CB.sai
bwa samse ../../CB4856_genome/Caenorhabditis_elegans_CB4856_v1.0.fa N2P-CB.sai ../2-FastQ_q1/N2-P_q1.fq.gz > N2P-CB.sam
samtools view -bS N2P-CB.sam > N2P-CB.unsorted.bam
samtools flagstat N2P-CB.unsorted.bam
samtools sort -@ 10 -o N2P-CB.bam N2P-CB.unsorted.bam 
samtools index -b N2P-CB.bam

bwa aln -o 0 -n 0 -t 10 ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa ../2-FastQ_q1/CB_q1.fq.gz > CB-N2.sai
bwa samse ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa CB-N2.sai ../2-FastQ_q1/CB_q1.fq.gz > CB-N2.sam
samtools view -bS CB-N2.sam > CB-N2.unsorted.bam
samtools flagstat CB-N2.unsorted.bam
samtools sort -@ 10 -o CB-N2.bam CB-N2.unsorted.bam
samtools index -b CB-N2.bam

bwa aln -o 0 -n 0 -t 10 ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa ../2-FastQ_q1/CB-P_q1.fq.gz > CBP-N2.sai
bwa samse ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa CBP-N2.sai ../2-FastQ_q1/CB-P_q1.fq.gz > CBP-N2.sam
samtools view -bS CBP-N2.sam > CBP-N2.unsorted.bam
samtools flagstat CBP-N2.unsorted.bam
samtools sort -@ 10 -o CBP-N2.bam CBP-N2.unsorted.bam
samtools index -b CBP-N2.bam

rm *.sai
rm *.sam
rm *.unsorted.bam