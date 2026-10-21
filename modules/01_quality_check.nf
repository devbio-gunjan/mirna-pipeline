process QUALITY_CHECK{

    errorStrategy 'terminate'

    publishDir "${params.outdir}/quality_check", mode: 'copy'

    container "${projectDir}/container/fastqc.sif"

    tag "${sample_id}"

    cpus 2
    memory '4 GB'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*_fastqc.zip"), emit: zip
    path("*_fastqc.html"), emit: html
    
    script:
    """
    fastqc --thread ${task.cpus} ${reads}
    """

}
