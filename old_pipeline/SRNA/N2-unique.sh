#!/bin/bash

#SBATCH --nodes=1   
#SBATCH --cpus-per-task=8
#SBATCH --error=/exports/people/andersenlab/mzy668/logs/job.%J.err
#SBATCH --output=/exports/people/andersenlab/mzy668/logs/job.%J.out
#SBATCH -J mz-N2unique
#SBATCH --workdir=/lscr2/andersenlab/mzy668/celegans/SRNA/7-N2unique
#SBATCH --mem=8096

samtools view -h ../4-BWA/N2.bam | awk '$1 ~ /^@/ || length($10) == 21' | samtools view -bS -> N2-21.unsorted.bam
samtools sort -@ 8 -o N2-21.bam N2-21.unsorted.bam
samtools index -b N2-21.bam
samtools view -b N2-21.bam IV:13500000-17200000 > N2-21pi.unsorted.bam 
samtools sort -@ 8 -o N2-21pi.bam N2-21pi.unsorted.bam 
samtools index -b N2-21pi.bam
samtools flagstat N2-21pi.bam 
bedtools bamtofastq -i N2-21pi.bam -fq N2-21pi.fq

rm *.unsorted.bam

bwa aln -o 0 -n 0 -t 8 ../../CB4856_genome/Caenorhabditis_elegans_CB4856_v1.0.fa N2-21pi.fq > N2pi-CB.sai
bwa samse ../../CB4856_genome/Caenorhabditis_elegans_CB4856_v1.0.fa N2pi-CB.sai N2-21pi.fq > N2pi-CB.sam
samtools view -bS N2pi-CB.sam > N2pi-CB.unsorted.bam
samtools flagstat N2pi-CB.unsorted.bam
samtools sort -@ 8 -o N2pi-CB.bam N2pi-CB.unsorted.bam
samtools index -b N2pi-CB.bam

rm *.sai
rm *.sam
rm *.unsorted.bam

samtools view -b -f 4 N2pi-CB.bam > N2pi-CB.unmapped.bam
bedtools bamtofastq -i N2pi-CB.unmapped.bam -fq N2pi-CB.unmapped.fq

bwa aln -o 0 -n 0 -t 8 ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa N2pi-CB.unmapped.fq > N2-unique.sai
bwa samse ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa N2-unique.sai N2pi-CB.unmapped.fq > N2-unique.sam
samtools view -bS N2-unique.sam > N2-unique.unsorted.bam
samtools flagstat N2-unique.unsorted.bam
samtools sort -@ 8 -o N2-unique.bam N2-unique.unsorted.bam
samtools index -b N2-unique.bam

rm *.sai
rm *.sam
rm *.unsorted.bam

bedtools bamtobed -i N2-unique.bam >N2-unique_temp1.bed 
cut -f1,2,3,6 N2-unique_temp1.bed  > N2-unique_temp2.bed
cat N2-unique_temp2.bed | uniq -c | awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $1}' > N2-unique_temp3.bed
cat N2-unique_temp3.bed | awk '{ if ($5 >= 2) print}' > N2-unique.bed
python N2_seqextract.py N2-unique.bed ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa > Seq_N2-unique.bed

rm *temp*bed