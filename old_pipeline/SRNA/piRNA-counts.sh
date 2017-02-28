#!/bin/bash

#SBATCH --nodes=1   
#SBATCH --cpus-per-task=8
#SBATCH --error=/exports/people/andersenlab/mzy668/logs/job.%J.err
#SBATCH --output=/exports/people/andersenlab/mzy668/logs/job.%J.out
#SBATCH -J mz-piRNAcounts
#SBATCH --workdir=/lscr2/andersenlab/mzy668/celegans/SRNA/10-piRNAcounts
#SBATCH --mem=8096

samtools view -h ../4-BWA/N2.bam | awk '$1 ~ /^@/ || length($10) == 21' | samtools view -bS -> N2-21.unsorted.bam
samtools sort -@ 8 -o N2-21.bam N2-21.unsorted.bam
samtools index -b N2-21.bam
samtools view -b N2-21.bam IV:13500000-17200000 > N2-21pi.unsorted.bam 
samtools sort -@ 8 -o N2-21pi.bam N2-21pi.unsorted.bam 
samtools index -b N2-21pi.bam
samtools flagstat N2-21pi.bam 

rm *.unsorted.bam
rm N2-21.bam
rm N2-21.bam.bai

samtools view -h ../6-BWAO/CB-N2.bam | awk '$1 ~ /^@/ || length($10) == 21' | samtools view -bS -> CB-N2-21.unsorted.bam
samtools sort -@ 8 -o CB-N2-21.bam CB-N2-21.unsorted.bam
samtools index -b CB-N2-21.bam
samtools view -b CB-N2-21.bam IV:13500000-17200000 > CB-N2-21pi.unsorted.bam 
samtools sort -@ 8 -o CB-N2-21pi.bam CB-N2-21pi.unsorted.bam 
samtools index -b CB-N2-21pi.bam
samtools flagstat CB-N2-21pi.bam 

rm *.unsorted.bam
rm CB-N2-21.bam
rm CB-N2-21.bam.bai

bedtools multicov -s -bams  N2-21pi.bam CB-N2-21pi.bam -bed ../../N2_genome/WS255.piRNA_exons.sorted.bed > BWA-counts-piRNAs.bed