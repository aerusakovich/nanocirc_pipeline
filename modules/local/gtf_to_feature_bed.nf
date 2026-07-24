process GTF_TO_FEATURE_BED {
    tag "gtf"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pybedtools:0.12.0--py39h475c85d_0' :
        'quay.io/biocontainers/pybedtools:0.12.0--py39h475c85d_0' }"

    input:
    path gtf

    output:
    path "gene.bed",     emit: gene_bed
    path "exon.bed",     emit: exon_bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk 'BEGIN{OFS="\\t"} !/^#/ && \$3=="gene" {
        n="."
        for(i=9;i<=NF;i++){if(\$i=="gene_id"){n=\$(i+1);gsub(/"/, "", n);break}}
        print \$1,\$4-1,\$5,n,0,\$7
    }' ${gtf} | sort -k1,1 -k2,2n > gene.bed

    awk 'BEGIN{OFS="\\t"} !/^#/ && \$3=="exon" {
        n="."
        for(i=9;i<=NF;i++){if(\$i=="gene_id"){n=\$(i+1);gsub(/"/, "", n);break}}
        print \$1,\$4-1,\$5,n,0,\$7
    }' ${gtf} | sort -k1,1 -k2,2n > exon.raw.bed

    # Different transcripts of the same gene can have overlapping exons.
    # classify_types() adds up the overlap bp for every matching exon.bed
    # line. Without merging, a locus with many exon-annotated isoforms
    # gets its bases counted more than once. That pushes exon_frac above
    # 1.0 and can wrongly turn an eicirna (part exon, part intron) into an
    # ecirna. Merge same-strand overlapping intervals first, so each base
    # is counted only once.
    bedtools merge -s -i exon.raw.bed -c 6 -o distinct \\
        | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3,".",0,\$4}' > exon.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>/dev/null | head -1 || awk -W version 2>/dev/null | head -1 || echo "unknown")
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}
