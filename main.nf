#!/usr/bin/env nextflow

// Read in the file listing the genomes to analyze
Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header:true)
    .map { row -> tuple(row.name, file(row.genome)) }
    .set { genome_ch }

params.output_name = "COMBINED_GENOMES"
params.output_folder = "./"

params.min_occurrence = "1"
params.minbit = "0.5"
params.distance = "euclidean"
params.linkage = "ward"
params.mcl_inflation = "2"

process makeGenomeDB {
    container "meren/anvio:5.5"
    cpus 4
    memory "8 GB"
    
    input:
    set name, file(fasta) from genome_ch
    
    output:
    set name, file("*db") into genomeDB_ch

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
    set name, file(db) from genomeDB_ch
    file anvio_cogs_tar
    
    output:
    file "${db}" into annotatedDB
    file "${db}.txt" into layer_txt

    afterScript "rm -rf *"

    """
#!/bin/bash

tar xvf ${anvio_cogs_tar}

anvi-run-ncbi-cogs -c "${db}" --num-threads 4 --cog-data-dir COGS_DIR

# Link the name to the database
echo \$name\\t${db} > ${db}.txt
    """
}

process combineGenomes {
    container "meren/anvio:5.5"
    cpus 4
    memory "8 GB"
    publishDir "${params.output_folder}"
    
    input:
    file db_list from annotatedDB.collect()
    file txt_list from layer_txt.collect()
    
    output:
    file "${params.output_name}.db" into combinedDB

    afterScript "rm -rf *"

    """
#!/bin/bash

echo -e "name\\tcontigs_db_path" > external-genomes.txt

for fp in ${txt_list}; do cat \$fp; done >> external-genomes.txt

anvi-gen-genomes-storage -e external-genomes.txt \
                         -o ${params.output_name}.db

    """
}

process panGenomeAnalysis {
    container "meren/anvio:5.5"
    cpus 4
    memory "8 GB"
    
    input:
    file combinedDB
    val output_name from params.output_name
    val min_occurrence from params.min_occurrence
    val minbit from params.minbit
    val distance from params.distance
    val linkage from params.linkage
    val mcl_inflation from params.mcl_inflation
    
    output:
    file "${params.output_name}-PAN.db" into panGenome

    afterScript "rm -rf *"

    """
#!/bin/bash

anvi-pan-genome -g ${combinedDB} \
                --project-name ${output_name} \
                --output-dir ./ \
                --num-threads 4 \
                --use-ncbi-blast \
                --min-occurrence ${min_occurrence} \
                --minbit ${minbit} \
                --distance ${distance} \
                --linkage ${linkage} \
                --mcl-inflation ${mcl_inflation}

    """
}

process addMetadata {
    container "meren/anvio:5.5"
    cpus 1
    memory "2 GB"
    publishDir "${params.output_folder}"
    
    input:
    file panGenome
    file sample_sheet from file("${params.sample_sheet}")
    
    output:
    file "${panGenome}"

    afterScript "rm -rf *"

    """
#!/bin/bash

# Strip out the genome file path
cat ${sample_sheet} | cut -f 2- > TEMP && mv TEMP ${sample_sheet}

anvi-import-misc-data ${sample_sheet} \
                      -p ${pan_genome} \
                      --target-data-table layers

    """
}

