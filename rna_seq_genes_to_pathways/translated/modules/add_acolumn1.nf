nextflow.enable.dsl=2

process ADD_ACOLUMN1 {
    
    container "python:3.7.16"
    publishDir "${params.outdir}/add_acolumn1"

    input:
    path input_file
    path script

    output:
    path "${input_file.simpleName}_extracol.txt", emit: outFile12

    script:
    """
    pip install numpy
    python \
    ${script} \
    ${input_file} \
    ${input_file.simpleName}_extracol.txt \
    "bool(c8<0.01) and bool(abs(c4)>0.58)" \
    9 \
    "int,str,str,float,float,float,float,float,float" \
    """

}
