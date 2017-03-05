#!/usr/bin/env nextflow

// Edit nextflow.configuration!

data_location=config.data_location
large_core=config.large_core
small_core=config.small_core

// ** - Recurse through subdirectories
fq_set = Channel.fromPath(data_location + "**/*.fastq.gz")
                .map { n -> [ n.getName(), n ] }


reference="WS255"
N2_project="PRJNA13758"
N2_prefix="ftp://ftp.wormbase.org/pub/wormbase/releases/${reference}/species/c_elegans/${N2_project}/"
CB_project="PRJNA275000"
CB_prefix="ftp://ftp.wormbase.org/pub/wormbase/releases/${reference}/species/c_elegans/${CB_project}/"

process fetch_reference {

    publishDir "output/", mode: 'copy', pattern: 'meta.pipeline.txt'
    
    output:
        file("N2_geneset.gtf.gz") into N2_geneset_gtf
        file("N2_reference.fa.gz") into N2_reference
        file("N2_annotations.gff3.gz") into N2_annotations
        file("CB_reference.fa.gz") into CB_reference
        file("meta.pipeline.txt")

    """
        echo '${N2_project}' > meta.pipeline.txt
        echo '${CB_project}' >> meta.pipeline.txt
        echo '${reference}' >> meta.pipeline.txt
        # Download N2 gene set
        curl ${N2_prefix}/c_elegans.${N2_project}.${reference}.canonical_geneset.gtf.gz > N2_geneset.gtf.gz

        #Download N2 reference genome
        curl ${N2_prefix}/c_elegans.${N2_project}.${reference}.genomic.fa.gz > N2_reference.fa.gz

        #Download N2 gff file
        curl ${N2_prefix}/c_elegans.${N2_project}.${reference}.annotations.gff3.gz > N2_annotations.gff3.gz

        #Download CB4856 reference genome
        curl ${CB_prefix}/c_elegans.${CB_project}.${reference}.genomic.fa.gz > CB_reference.fa.gz
        
    """
}
N2_geneset_gtf.into { N2_geneset_gtf_stringtie }
N2_reference.into { N2_reference_BWA; N2_reference_seqextract; N2_reference_blast }
CB_reference.into { CB_reference_BWA; CB_reference_seqextract; CB_reference_blast }

process bwa_index_N2 {

    cpus small_core

    tag { name }

    input:
        file("N2_reference.fa.gz") from N2_reference_BWA

    output:
       file "N2_reference.*" into N2_bwaindex

    """
        zcat N2_reference.fa.gz > N2_reference.fa
        bwa index N2_reference.fa   
    """
}
N2_bwaindex.into { N2_bwaindex_1; N2_bwaindex_21 }

//CB genome headers require clipping of extraneous chromosome length information
process bwa_index_CB {

    cpus small_core

    tag { name }

    input:
        file("CB_reference.fa.gz") from CB_reference_BWA

    output:
       file "CB_reference.*" into CB_bwaindex

    """ 
        zcat CB_reference.fa.gz | awk '{print \$1; }' > CB_reference.fa
        bwa index CB_reference.fa
    """
}
CB_bwaindex.into { CB_bwaindex_1; CB_bwaindex_21 }


adapters = file("auxillary/TruSeq3-SE.fa")

process trimmomatic {

    cpus small_core

    tag { name }

    input:
        set val(name), file(reads) from fq_set

    output:
        file(name_out) into trimmed_reads

    script:
    name_out = name.replace('.fastq.gz', '_trim.fq.gz')

    """
        trimmomatic SE -phred33 -threads ${small_core} ${reads} ${name_out} ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 
    """
}

process align {

    cpus large_core

    tag { reads }

    input:
        file reads from trimmed_reads
        file CB_bwaindex from CB_bwaindex_1.first()
        file N2_bwaindex from N2_bwaindex_1.first()

    output:
        set val(strain_id), val(phosphate_id), file("${fa_prefix}.bam"), file("${fa_prefix}.bam.bai") into bwa_bams

    script:
        m = reads =~ /\w+-([^_P]{2}).(P{0,1})_.*/
        strain_id = m[0][1]
        phosphate_id = m[0][2]
        fa_prefix = reads[0].toString() - ~/(_trim)(\.fq\.gz)$/

        println strain_id
        println phosphate_id

        if (phosphate_id == "P")
            """
            bwa aln -o 0 -n 0 -t ${large_core} N2_reference.fa ${reads} > ${fa_prefix}.sai
            bwa samse N2_reference.fa ${fa_prefix}.sai ${reads} > ${fa_prefix}.sam
            samtools view -bS ${fa_prefix}.sam > ${fa_prefix}.unsorted.bam
            samtools flagstat ${fa_prefix}.unsorted.bam
            samtools sort -@ ${large_core} -o ${fa_prefix}.bam ${fa_prefix}.unsorted.bam
            samtools index -b ${fa_prefix}.bam
            """
        else if (strain_id == "N2" && phosphate_id == "")
            """
            bwa aln -o 0 -n 0 -t ${large_core} N2_reference.fa ${reads} > ${fa_prefix}.sai
            bwa samse N2_reference.fa ${fa_prefix}.sai ${reads} > ${fa_prefix}.sam
            samtools view -bS ${fa_prefix}.sam > ${fa_prefix}.unsorted.bam
            samtools flagstat ${fa_prefix}.unsorted.bam
            samtools sort -@ ${large_core} -o ${fa_prefix}.bam ${fa_prefix}.unsorted.bam
            samtools index -b ${fa_prefix}.bam
            """
        else if (strain_id == "CB" && phosphate_id == "")
            """
            bwa aln -o 0 -n 0 -t ${large_core} CB_reference.fa ${reads} > ${fa_prefix}.sai
            bwa samse CB_reference.fa ${fa_prefix}.sai ${reads} > ${fa_prefix}.sam
            samtools view -bS ${fa_prefix}.sam > ${fa_prefix}.unsorted.bam
            samtools flagstat ${fa_prefix}.unsorted.bam
            samtools sort -@ ${large_core} -o ${fa_prefix}.bam ${fa_prefix}.unsorted.bam
            samtools index -b ${fa_prefix}.bam
            """
        else
            """
            """
}

bwa_bams.into { P_unfiltered; N2_unfiltered; CB_unfiltered ; N2CB_unfiltered }

//it[0] = sample_id, it[1] = phosphate_id
P_unfiltered.filter{ it[1] == "P"}.into {P_filtered} 

//N2_unfiltered.filter{ it[1] == "" && it[0] == "N2"}.into {N2_filtered}
//CB_unfiltered.filter{ it[1] == "" && it[0] == "CB"}.into {CB_filtered}
//Could combine these and split later, they are being fed to common processes
N2CB_unfiltered.filter{ it[1] == ""}.into {N2CB_filtered}

////////////////////////
// P-treated smRNA -> 22G (filter) -> quantify abundance (antisense to features) 
////////////////////////

process _22G_filter_bams {

    tag { sample_id }

    input:
        set val(strain_id), val(phosphate_id), file(sample_bam), file(sample_bai) from P_filtered

    output:
        set val(sample_id), file("${fa_prefix}-22G.bam"), file("${fa_prefix}-22G.bam.bai") into _22g_bams

    script:
        fa_prefix = sample_bam[0].toString() - ~/(\.bam)$/
        m = fa_prefix =~ /\w+-([^_]+)_.*/
        sample_id = m[0][1]

    """
        samtools view -h ${fa_prefix}.bam | awk '\$1 ~ /^@/ || length(\$10) == 22 && \$10 ~ /^G/' | samtools view -bS -> ${fa_prefix}-22G.bam
        samtools index -b ${fa_prefix}-22G.bam
    """

}

_22g_joint_bams = _22g_bams.groupTuple()

process _22G_join_bams {

    publishDir "output/22Gbam", mode: 'copy'

    cpus small_core

    tag { sample_id }

    input:
        set val(sample_id), file(bam), file(bai) from _22g_joint_bams

    output:
        set val(sample_id), file("${sample_id}.bam"), file("${sample_id}.bam.bai") into _22g_merged_bams

    """
        samtools merge -f ${sample_id}.unsorted.bam ${bam}
        samtools sort -@ ${small_core} -o ${sample_id}.bam ${sample_id}.unsorted.bam
        samtools flagstat ${sample_id}.bam
        samtools index -b ${sample_id}.bam

    """
}

//Fix to check for anti-sense
process stringtie_22G_counts {

    publishDir "output/22Gexpression", mode: 'copy'

    cpus small_core

    tag { sample_id }

    input:
        set val(sample_id), file(bam), file(bai) from _22g_merged_bams
        file("geneset.gtf.gz") from N2_geneset_gtf_stringtie.first()

    output:
        file("${sample_id}/*") into stringtie_exp

    """ 
        zcat geneset.gtf.gz > geneset.gtf
        stringtie -p ${small_core} -G geneset.gtf -A ${sample_id}/${sample_id}_abund.tab -e -B -o ${sample_id}/${sample_id}_expressed.gtf ${bam}
    """
}

prepDE = file("scripts/prepDE.py")

process stringtie_22G_table_counts {

    echo true

    publishDir "output/22Gdiffexp", mode: 'copy'

    cpus small_core

    tag { sample_id }

    input:
        val(sample_file) from stringtie_exp.toSortedList()

    output:
        file ("gene_count_matrix.csv") into gene_count_matrix
        file ("transcript_count_matrix.csv") into transcript_count_matrix

    """
        for i in ${sample_file.flatten().join(" ")}; do
            bn=`basename \${i}`
            full_path=`dirname \${i}`
            sample_name=\${full_path##*/}
            echo "\${sample_name} \${i}"
            mkdir -p expression/\${sample_name}
            ln -s \${i} expression/\${sample_name}/\${bn}
        done;
        python ${prepDE} -i expression -l 22 -g gene_count_matrix.csv -t transcript_count_matrix.csv

    """
}




////////////////////////
// Find N2 and CB unique 21-mers (non-phosphatase-treated samples)
////////////////////////

// Filter aligned bams (N2->N2, CB->CB) to aligned 21-mers 
process N2CB_filter_bams {

    tag { sample_bam }

    input:
        set val(strain_id), val(phosphate_id), file(sample_bam), file(sample_bai) from N2CB_filtered


    output:
        set val(strain_id), val(sample_id), val(fa_prefix), file("${fa_prefix}-21.bam"), file("${fa_prefix}-21.bam.bai") into _21_bams

    script:
        fa_prefix = sample_bam[0].toString() - ~/(\.bam)$/
        m = fa_prefix =~ /\w+-([^_]+)_.*/
        sample_id = m[0][1]

    """
        samtools view -h ${fa_prefix}.bam | awk '\$1 ~ /^@/ || length(\$10) == 21' | samtools view -bS -> ${fa_prefix}-21.bam
        samtools index -b ${fa_prefix}-21.bam
    """

}

// Extract aligned 21-mers as fastq
process N2CB_extract_21mers {

    tag { bam }

    input:
        set val(strain_id), val(sample_id), val(fa_prefix), file(bam), file(bai) from _21_bams

    output:
        set val(strain_id), val(sample_id), val(fa_prefix), file("${fa_prefix}-21.fq") into _21_fqs

    """
        bedtools bamtofastq -i ${bam} -fq ${fa_prefix}-21.fq
    """
}

// Align 21-mers to opposite strain (N2->CB, CB->N2), retain those that do not align -> bams containing strain-unique 21mers
process N2CB_unique_21mers {

    cpus large_core

    tag { reads }

    input:
        set val(strain_id), val(sample_id), val(fa_prefix), file(reads) from _21_fqs
        file CB_bwaindex from CB_bwaindex_21.first()
        file N2_bwaindex from N2_bwaindex_21.first()
        
    output:
        set val(sample_id), file("${fa_prefix}-unique.bam"), file("${fa_prefix}-unique.bam.bai") into _21unique_bams

    script:

        if (strain_id == "N2")
            """
            bwa aln -o 0 -n 0 -t ${large_core} CB_reference.fa ${reads} > ${fa_prefix}-opp.sai
            bwa samse CB_reference.fa ${fa_prefix}-opp.sai ${reads} > ${fa_prefix}-opp.sam
            samtools view -bS ${fa_prefix}-opp.sam > ${fa_prefix}-opp.unsorted.bam
            samtools flagstat ${fa_prefix}-opp.unsorted.bam
            samtools sort -@ ${large_core} -o ${fa_prefix}-opp.bam ${fa_prefix}-opp.unsorted.bam
            samtools index -b ${fa_prefix}-opp.bam

            samtools view -b -f 4 ${fa_prefix}-opp.bam > ${fa_prefix}-opp.unmapped.bam
            bedtools bamtofastq -i ${fa_prefix}-opp.unmapped.bam -fq ${fa_prefix}-opp.unmapped.fq

            bwa aln -o 0 -n 0 -t ${large_core} N2_reference.fa ${fa_prefix}-opp.unmapped.fq > ${fa_prefix}-unique.sai
            bwa samse N2_reference.fa ${fa_prefix}-unique.sai ${fa_prefix}-opp.unmapped.fq > ${fa_prefix}-unique.sam
            samtools view -bS ${fa_prefix}-unique.sam > ${fa_prefix}-unique.unsorted.bam
            samtools flagstat ${fa_prefix}-unique.unsorted.bam
            samtools sort -@ ${large_core} -o ${fa_prefix}-unique.bam ${fa_prefix}-unique.unsorted.bam
            samtools index -b ${fa_prefix}-unique.bam

            """
        else if (strain_id == "CB")
            """
            bwa aln -o 0 -n 0 -t ${large_core} N2_reference.fa ${reads} > ${fa_prefix}-opp.sai
            bwa samse N2_reference.fa ${fa_prefix}-opp.sai ${reads} > ${fa_prefix}-opp.sam
            samtools view -bS ${fa_prefix}-opp.sam > ${fa_prefix}-opp.unsorted.bam
            samtools flagstat ${fa_prefix}-opp.unsorted.bam
            samtools sort -@ ${large_core} -o ${fa_prefix}-opp.bam ${fa_prefix}-opp.unsorted.bam
            samtools index -b ${fa_prefix}-opp.bam

            samtools view -b -f 4 ${fa_prefix}-opp.bam > ${fa_prefix}-opp.unmapped.bam
            bedtools bamtofastq -i ${fa_prefix}-opp.unmapped.bam -fq ${fa_prefix}-opp.unmapped.fq

            bwa aln -o 0 -n 0 -t ${large_core} CB_reference.fa ${fa_prefix}-opp.unmapped.fq > ${fa_prefix}-unique.sai
            bwa samse CB_reference.fa ${fa_prefix}-unique.sai ${fa_prefix}-opp.unmapped.fq > ${fa_prefix}-unique.sam
            samtools view -bS ${fa_prefix}-unique.sam > ${fa_prefix}-unique.unsorted.bam
            samtools flagstat ${fa_prefix}-unique.unsorted.bam
            samtools sort -@ ${large_core} -o ${fa_prefix}-unique.bam ${fa_prefix}-unique.unsorted.bam
            samtools index -b ${fa_prefix}-unique.bam
            
            """
        else
            """
            """
}



//Merge bams containing species-unique 21mers by sample_id
//Todo: change order later: merge first, then check for unique
_21unique_joint_bams = _21unique_bams.groupTuple()
process N2CB_21unique_join_bams {

    publishDir "output/bam_21unique", mode: 'copy'

    cpus small_core

    tag { sample_id }

    input:
        set val(sample_id), file(bam), file(bai) from _21unique_joint_bams

    output:
        set val(sample_id), file("${sample_id}-unique.bam"), file("${sample_id}-unique.bam.bai") into _21unique_merged_bams

    """
        samtools merge -f ${sample_id}-unique.unsorted.bam ${bam}
        samtools sort -@ ${small_core} -o ${sample_id}-unique.bam ${sample_id}-unique.unsorted.bam
        samtools flagstat ${sample_id}-unique.bam
        samtools index -b ${sample_id}-unique.bam

    """
}

//For each sample, create bed file with locations of unique 21mers and their counts
process N2CB_21unique_bed {

    //publishDir "output/bed_21unique", mode: 'copy'

    cpus small_core

    tag { sample_id }

    input:
        set val(sample_id), file(bam), file(bai) from _21unique_merged_bams

    output:
        set val(strain_id), val(sample_id), file("${sample_id}-unique.bed") into _21unique_beds

    script:
        m = sample_id =~ /(.{2}).*/
        strain_id = m[0][1]

    //minium count per sample = 2 (change if necessary)
    """
        bedtools bamtobed -i ${bam} > temp.bed
        cut -f1,2,3,6 temp.bed  > temp2.bed
        cat temp2.bed | uniq -c | awk '{print \$2 "\t" \$3 "\t" \$4 "\t" \$5 "\t" \$1 "\t" "${sample_id}"}' > temp3.bed
        cat temp3.bed | awk '{ if (\$5 >= 2) print}' > ${sample_id}-unique.bed

    """
}


//For each strain, retain only strain-specific 21mers that are present in all replicate samples
_21unique_joint_beds = _21unique_beds.groupTuple() //group replicate samples by strain

process N2CB_21unique_bed_strainmerge {

    publishDir "output/bed_21unique", mode: 'copy'

    cpus small_core

    tag { sample_id }

    input:
        set val(strain_id), val(sample_id), file(bed) from _21unique_joint_beds

    output:
        file("${strain_id}-unique.bed") into dead_end
        set val(strain_id), file ("${strain_id}-shared_unique.bed") into _21unique_shared

    script:

    """
        cat ${bed} > temp.bed 
        sortBed -i temp.bed > ${strain_id}-unique.bed
        cat ${strain_id}-unique.bed | awk '{print \$1 "\t" \$2 "\t" \$3 "\t" \$4}' | uniq -c | awk '{ if (\$1 == 4) print \$2 "\t" \$3 "\t" \$4 "\t" \$5 }' > ${strain_id}-shared_unique.bed

    """
}
_21unique_shared.into { N2_21unique_shared_unfiltered ; CB_21unique_shared_unfiltered }

//it[0] = sample_id, it[1] = phosphate_id
N2_21unique_shared_unfiltered.filter{ it[0] == "N2"}.into {N2_21unique_shared} 
CB_21unique_shared_unfiltered.filter{ it[0] == "CB"}.into {CB_21unique_shared} 


//For each strain, extract sequences of replicate-shared 21mers 
N2seqextract = file("scripts/N2_seqextract.py")
process N2_21unique_seqextract {

    publishDir "output/seq_21unique", mode: 'copy'

    cpus small_core

    tag { strain_id }

    input:
        set val(strain_id), file (bed) from N2_21unique_shared
        file("N2_reference.fa.gz") from N2_reference_seqextract

    output:
        file ("${strain_id}_unique.fa") into N2_21unique_fasta


    """
        zcat N2_reference.fa.gz > N2_reference.fa
        python ${N2seqextract} ${bed} N2_reference.fa > ${strain_id}_unique.fa

    """
}

CBseqextract = file("scripts/CB_seqextract.py")
process CB_21unique_seqextract {

    publishDir "output/seq_21unique", mode: 'copy'

    cpus small_core

    tag { strain_id }

    input:
        set val(strain_id), file (bed) from CB_21unique_shared
        file("CB_reference.fa.gz") from CB_reference_seqextract

    output:
        file ("${strain_id}_unique.fa") into CB_21unique_fasta


    """
        zcat CB_reference.fa.gz | awk '{print \$1; }' > CB_reference.fa
        python ${CBseqextract} ${bed} CB_reference.fa > ${strain_id}_unique.fa

    """
}



//Index N2 BLAST database 
process blast_formatdb {

    cpus small_core

    input:
        file("N2_reference.fa.gz") from N2_reference_blast

    output:
       file "N2*" into N2_blastindex

    """
        zcat N2_reference.fa.gz > N2_reference.fa
        makeblastdb -in N2_reference.fa -dbtype nucl -out N2ref
    """
}
N2_blastindex.into { N2_blastindex_1 ; N2_blastindex_2 }


//Blast N2 unique 21mers
process blast_N2_21 {

    publishDir "output/blast", mode: 'copy'

    cpus small_core

    input:
        file N2_blastindices from N2_blastindex_1.first()
        file("N2_unique.fa") from N2_21unique_fasta

    output:
        file "N2_blast.txt" into N2_blastout

    """
        blastn -db N2ref -query N2_unique.fa -outfmt 6 -task blastn-short -qcov_hsp_perc 80 -out N2_blast.txt  
    """
}

//Blast CB unique 21mers
process blast_CB_21 {

    publishDir "output/blast", mode: 'copy'

    cpus small_core

    input:
        file N2_blastindices from N2_blastindex_2.first()
        file("CB_unique.fa") from CB_21unique_fasta

    output:
        file "CB_blast.txt" into CB_blastout

    """
        blastn -db N2ref -query CB_unique.fa -outfmt 6 -task blastn-short -qcov_hsp_perc 80 -out CB_blast.txt   
    """
}





