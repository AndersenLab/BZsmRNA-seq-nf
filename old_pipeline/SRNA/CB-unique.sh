#!/bin/bash

#SBATCH --nodes=1   
#SBATCH --cpus-per-task=8
#SBATCH --error=/exports/people/andersenlab/mzy668/logs/job.%J.err
#SBATCH --output=/exports/people/andersenlab/mzy668/logs/job.%J.out
#SBATCH -J mz-CBunique
#SBATCH --workdir=/lscr2/andersenlab/mzy668/celegans/SRNA/8-CBunique
#SBATCH --mem=8096

samtools view -h ../4-BWA/CB.bam | awk '$1 ~ /^@/ || length($10) == 21' | samtools view -bS -> CB-21.unsorted.bam
samtools sort -@ 8 -o CB-21.bam CB-21.unsorted.bam
samtools index -b CB-21.bam
samtools view -b CB-21.bam IV:13500000-17200000 > CB-21pi.unsorted.bam 
samtools sort -@ 8 -o CB-21pi.bam CB-21pi.unsorted.bam 
samtools index -b CB-21pi.bam
samtools flagstat CB-21pi.bam 
bedtools bamtofastq -i CB-21pi.bam -fq CB-21pi.fq

rm *.unsorted.bam

bwa aln -o 0 -n 0 -t 8 ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa CB-21pi.fq > CBpi-N2.sai
bwa samse ../../N2_genome/c_elegans.PRJNA13758.WS255.genomic.fa CBpi-N2.sai CB-21pi.fq > CBpi-N2.sam
samtools view -bS CBpi-N2.sam > CBpi-N2.unsorted.bam
samtools flagstat CBpi-N2.unsorted.bam
samtools sort -@ 8 -o CBpi-N2.bam CBpi-N2.unsorted.bam
samtools index -b CBpi-N2.bam

rm *.sai
rm *.sam
rm *.unsorted.bam

samtools view -b -f 4 CBpi-N2.bam > CBpi-N2.unmapped.bam
bedtools bamtofastq -i CBpi-N2.unmapped.bam -fq CBpi-N2.unmapped.fq

bwa aln -o 0 -n 0 -t 8 ../../CB4856_genome/Caenorhabditis_elegans_CB4856_v1.0.fa CBpi-N2.unmapped.fq > CB-unique.sai
bwa samse ../../CB4856_genome/Caenorhabditis_elegans_CB4856_v1.0.fa CB-unique.sai CBpi-N2.unmapped.fq > CB-unique.sam
samtools view -bS CB-unique.sam > CB-unique.unsorted.bam
samtools flagstat CB-unique.unsorted.bam
samtools sort -@ 8 -o CB-unique.bam CB-unique.unsorted.bam
samtools index -b CB-unique.bam

rm *.sai
rm *.sam
rm *.unsorted.bam

bedtools bamtobed -i CB-unique.bam >CB-unique_temp1.bed 
cut -f1,2,3,6 CB-unique_temp1.bed  > CB-unique_temp2.bed
cat CB-unique_temp2.bed | uniq -c | awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $1}' > CB-unique_temp3.bed
cat CB-unique_temp3.bed | awk '{ if ($5 >= 2) print}' > CB-unique.bed
python CB_seqextract.py CB-unique.bed ../../CB4856_genome/Caenorhabditis_elegans_CB4856_v1.0.fa > Seq_CB-unique.bed

rm *temp*bed