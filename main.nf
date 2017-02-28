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
//geneset_gtf.into { geneset_hisat; geneset_stringtie }


process bwa_index_N2 {

    cpus small_core

    tag { name }

    input:
        file("N2_reference.fa.gz") from N2_reference

    output:
       file "N2_reference.*" into N2_bwaindex

    """
        gzcat N2_reference.fa.gz > N2_reference.fa
        bwa index N2_reference.fa   
    """
}

process bwa_index_CB {

    cpus small_core

    tag { name }

    input:
        file("CB_reference.fa.gz") from CB_reference

    output:
       file "CB_reference.*" into CB_bwaindex

    """ 
        gzcat CB_reference.fa.gz > CB_reference.fa
        bwa index CB_reference.fa
    """
}

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
        trimmomatic SE -phred33 -threads ${small_core} ${reads} ${name_out} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 
    """
}

process align {

    cpus large_core

    tag { reads }

    input:
        file reads from trimmed_reads
        file CB_bwaindex from CB_bwaindex.first()
        file N2_bwaindex from N2_bwaindex.first()

    output:
        set val(sample_id), val(phosphate_id), file("${fa_prefix}.bam"), file("${fa_prefix}.bam.bai") into bwa_bams

    script:
        m = reads =~ /\w+-([^_P]{2}).(P{0,1})_.*/
        sample_id = m[0][1]
        phosphate_id = m[0][2]
        fa_prefix = reads[0].toString() - ~/(_trim)(\.fq\.gz)$/

        println sample_id
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
        else if (sample_id == "N2" && phosphate_id == "")
            """
            bwa aln -o 0 -n 0 -t ${large_core} N2_reference.fa ${reads} > ${fa_prefix}.sai
            bwa samse N2_reference.fa ${fa_prefix}.sai ${reads} > ${fa_prefix}.sam
            samtools view -bS ${fa_prefix}.sam > ${fa_prefix}.unsorted.bam
            samtools flagstat ${fa_prefix}.unsorted.bam
            samtools sort -@ ${large_core} -o ${fa_prefix}.bam ${fa_prefix}.unsorted.bam
            samtools index -b ${fa_prefix}.bam
            """
        else if (sample_id == "CB" && phosphate_id == "")
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

bwa_bams.into { P_unfiltered; N2_unfiltered; CB_unfiltered }

P_unfiltered.filter{ it[1] == "P"}.into {P_filtered} //it[0] = sample_id
N2_unfiltered.filter{ it[1] == "" && it[0] == "N2"}.into {N2_filtered}
CB_unfiltered.filter{ it[1] == "" && it[0] == "CB"}.into {CB_filtered}


process _22Gexpression {

    input:
        set val(sample_id), val(phosphate_id), file(sample_bam), file(sample_bai) from P_filtered


    """
        samtools view -h ${sample_bam} | awk '\$1 ~ /^@/ || length(\$10) == 22 && \$10 ~ /^G/' | samtools view -bS -> ${sample_bam}-22G.bam
        samtools index -b ${sample_bam}-22G.bam
    """


}


