process setupNCBIcogs {
    container "${params.container__anvio}"
    label "cpu_high"
    publishDir "${params.output_folder}"
    
    output:
    file "COGS_DIR.tar"

    """#!/bin/bash
set -e

anvi-setup-ncbi-cogs --num-threads ${task.cpus} --cog-data-dir COGS_DIR --just-do-it --reset
tar cvf COGS_DIR.tar COGS_DIR
    """
}
