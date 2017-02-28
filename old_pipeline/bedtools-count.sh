#!/bin/bash

#SBATCH --nodes=1   
#SBATCH --cpus-per-task=6
#SBATCH --error=/exports/people/andersenlab/mzy668/logs/job.%J.err
#SBATCH --output=/exports/people/andersenlab/mzy668/logs/job.%J.out
#SBATCH -J bedtools-count
#SBATCH --workdir=/lscr2/andersenlab/mzy668/celegans/Counts
#SBATCH --mem=8096

bedtools multicov -s -bams ../RNA/BWA/N2.bam ../RNA/BWA/CB-N2.bam -bed ../N2_genome/WS255.transcripts_exons.sorted.bed > BWA-counts-RNA.bed

bedtools multicov -S -bams ../SRNA/4-BWA/N2P-22G.bam ../SRNA/6-BWAO/CBP-N2-22G.bam -bed ../N2_genome/WS255.transcripts_exons.sorted.bed > BWA-counts-SRNA.bed

	