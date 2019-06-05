#!/usr/bin/env nextflow

// Read in the file listing the genomes to analyze
genome_ch = Channel
    .from(file("${params.sample_sheet}").readLines())
    .map { line -> file(line) }

process makeGenomeDB {
    container "meren/anvio:5.5"
    cpus 4
    memory "8 GB"
    
    input:
    file fasta from genome_ch
    
    output:
    file "*db" into genomeDB_ch

    afterScript "rm -rf *"

    """
#!/bin/bash

fasta=${fasta}

# Decompress the FASTA if it is compressed
gzip -t \$fasta && gunzip \$fasta && fasta=\$(echo \$fasta | sed 's/.gz//')

# The file ending must be "fa" or "fasta"
if [[ \$fasta =~ ".fa" ]] || [[ \$fasta =~ ".fasta" ]]; then
    pass
else
    mv \$fasta \$fasta.fasta
    fasta=\$fasta.fasta
fi

anvi-script-FASTA-to-contigs-db \$fasta
    """
}

process setupNCBIcogs {
    container "meren/anvio:5.5"
    cpus 4
    memory "8 GB"
    
    output:
    file "COGS_DIR.tar" into anvio_cogs_tar

    afterScript "rm -rf *"

    """
#!/bin/bash

anvi-setup-ncbi-cogs --num-threads 4 --cog-data-dir COGS_DIR --just-do-it

tar cvf COGS_DIR.tar COGS_DIR
    """
}

process annotateGenes {
    container "meren/anvio:5.5"
    cpus 4
    memory "8 GB"
    
    input:
    file db from genomeDB_ch
    file anvio_cogs_tar
    
    output:
    file "${db}" into annotatedDB_ch

    afterScript "rm -rf *"

    """
#!/bin/bash

tar xvf ${anvio_cogs_tar}

anvi-run-ncbi-cogs -c "${db}" --num-threads 4 --cog-data-dir COGS_DIR
    """
}

process combineGenomes {
    container "meren/anvio:5.5"
    cpus 4
    memory "8 GB"
    
    input:
    file db_list from annotatedDB_ch.collect()
    
    output:
    file "COMBINED-GENOMES.db" into combinedDB_ch

    afterScript "rm -rf *"

    """
#!/bin/bash

echo name\\tcontigs_db_path > external-genomes.txt

for db in ${db_list}; do
    echo \$(echo \$db | sed 's/.db//' | sed 's/.fasta//' | sed 's/.fna//' | sed 's/.gz//')\\t\$db >> external-genomes.txt
done

anvi-gen-genomes-storage -e external-genomes.txt \
                         -o COMBINED-GENOMES.db

    """
}

