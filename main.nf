#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process parseSampleSheet {
    container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
    label "io_limited"
    
    input:
    path sample_sheet_csv
    
    output:
    path "${sample_sheet_csv}"

    """#!/usr/bin/env python3
import pandas as pd

print("Reading in ${sample_sheet_csv}")

df = pd.read_csv("${sample_sheet_csv}", sep=",")

print("Read in %d rows and %d columns" % (df.shape[0], df.shape[1]))

for k in ["name", "genome"]:
    assert k in df.columns.values, "Must provide a column '%s' in the sample sheet" % k

# Strip away all whitespace and carriage returns
df = df.applymap(str).applymap(lambda s: s.strip())

print("Writing out sanitized sample sheet")
df.to_csv("${sample_sheet_csv}", index=None, sep="\\t")
print("Done")
    """
}

process makeGenomeDB {
    container "${params.container__anvio}"
    label "mem_medium"
    
    input:
    tuple val(name), path(fasta)
    
    output:
    tuple val(name), path("${name}.db")

    """#!/bin/bash
set -e

fasta=${fasta}
# Decompress the FASTA if it is compressed
gzip -t \$fasta && gunzip \$fasta && fasta=\$(echo \$fasta | sed 's/.gz//')
# The file ending must be "fa", "fsa", or "fasta"
if [[ \$fasta =~ ".fa" ]] || [[ \$fasta =~ ".fsa" ]] || [[ \$fasta =~ ".fasta" ]]; then
    echo "No need to rename \$fasta"
else
    mv \$fasta \$fasta.fasta
    fasta=\$fasta.fasta
fi

# Reformat the FASTA to sanitize deflines
anvi-script-reformat-fasta --seq-type NT --simplify-names -l 0 -o \$fasta.clean.fasta \$fasta

# Make the genome database
anvi-gen-contigs-database -f \$fasta.clean.fasta -n ${name} -o ${name}.db
    """
}

process setupNCBIcogs {
    container "${params.container__anvio}"
    label "cpu_high"
    
    output:
    path "COGS_DIR.tar"

    """#!/bin/bash
set -e

anvi-setup-ncbi-cogs --num-threads ${task.cpus} --cog-data-dir COGS_DIR --just-do-it --reset
tar cvf COGS_DIR.tar COGS_DIR
    """
}


process annotateGenes {
    container "${params.container__anvio}"
    label "cpu_high"
    
    input:
    tuple val(name), path(db)
    each path(anvio_cogs_tar)
    
    output:
    path "${db}"

    """#!/bin/bash
set -e

tar xvf ${anvio_cogs_tar}
anvi-run-ncbi-cogs -c "${db}" --num-threads ${task.cpus} --cog-data-dir COGS_DIR
    """
}

process linkGeneName {
    container "${params.container__anvio}"
    label "mem_medium"
    
    input:
    tuple val(name), path(db)
    
    output:
    path "${db}.txt"

    """#!/bin/bash
set -e

# Link the name to the database
echo -e ${name},${db} | tr ',' '\\t' > ${db}.txt
    """
}

process combineGenomes {
    container "${params.container__anvio}"
    label "mem_medium"
    publishDir "${params.output_folder}", mode: "copy", overwrite: true
    
    input:
    path db_list
    path txt_list
    
    output:
    path "${params.output_name}-GENOMES.db"
    path "external-genomes.txt" 

    """#!/bin/bash
set -e

echo -e "name\\tcontigs_db_path" > external-genomes.txt
for fp in ${txt_list}; do cat \$fp; done >> external-genomes.txt
cat external-genomes.txt
anvi-gen-genomes-storage -e external-genomes.txt \
                         -o ${params.output_name}-GENOMES.db
    """
}

process panGenomeAnalysis {
    container "${params.container__anvio}"
    label "cpu_high"
    
    input:
    path combinedDB
    
    output:
    path "${params.output_name}-PAN.db"

    """#!/bin/bash
set -e

anvi-pan-genome -g ${combinedDB} \
                --project-name ${params.output_name} \
                --output-dir ./ \
                --num-threads ${task.cpus} \
                --use-ncbi-blast \
                --min-occurrence ${params.min_occurrence} \
                --minbit ${params.minbit} \
                --distance ${params.distance} \
                --linkage ${params.linkage} \
                --mcl-inflation ${params.mcl_inflation}
    """
}

process getSequencesForGCs {
    container "${params.container__anvio}"
    label "mem_medium"
    publishDir "${params.output_folder}", mode: "copy", overwrite: true
    
    input:
    path panGenome
    path combinedDB
    
    output:
    path "${params.output_name}.gene_clusters.fastp"
    path "${params.output_name}.gene_clusters.fasta"

    """#!/bin/bash
set -e
    
anvi-get-sequences-for-gene-clusters \
    -p ${panGenome} \
    -g ${combinedDB} \
    -o ${params.output_name}.gene_clusters.fastp \
    --just-do-it

anvi-get-sequences-for-gene-clusters \
    -p ${panGenome} \
    -g ${combinedDB} \
    -o ${params.output_name}.gene_clusters.fasta \
    --report-DNA-sequences \
    --just-do-it
    """
}
    
process addMetadata {
    container "${params.container__anvio}"
    label "io_limited"
    publishDir "${params.output_folder}", mode: "copy", overwrite: true
    
    input:
    path panGenome
    path combinedDB
    path sample_sheet
    
    output:
    path "${panGenome}"

    """#!/bin/bash
set -e


# Strip out the genome file path
cat ${sample_sheet} | cut -f 2- > TEMP && mv TEMP ${sample_sheet}
echo "Printing the reformatted sample sheet:"
cat ${sample_sheet}

if (( \$(head -1 ${sample_sheet} | tr "\\t" "\\n" | wc -l ) > 1 )); then
    echo ""
    echo "Adding metadata"
    anvi-import-misc-data ${sample_sheet} \
                          -p ${panGenome} \
                          --target-data-table layers
else
    echo ""
    echo "Not adding any metadata, none provided"

fi
    """
}

process enrichFunctions{
    container "${params.container__anvio}"
    label "mem_medium"
    publishDir "${params.output_folder}", mode: "copy", overwrite: true
    
    input:
    path panGenome
    path combinedDB
    val category_name
    
    output:
    path "${params.output_name}-enriched-functions-${category_name}.txt"
    path "${params.output_name}-functions-occurrence.txt"

    script:

    if (params.gene_enrichment == false)
        """#!/bin/bash
        anvi-compute-functional-enrichment -p ${panGenome} \
                                            -g ${combinedDB} \
                                            --category-variable ${category_name} \
                                            --annotation-source COG20_FUNCTION \
                                            -o "${params.output_name}-enriched-functions-${category_name}.txt" \
                                            --functional-occurrence-table-output "${output_name}-functions-occurrence.txt"
    """

    else
        """#!/bin/bash
        anvi-compute-functional-enrichment -p ${panGenome} \
                                            -g ${combinedDB} \
                                            --category-variable ${category_name} \
                                            --annotation-source IDENTITY \
                                            --include-gc-identity-as-function \
                                            -o "${params.output_name}-enriched-functions-${category_name}.txt" \
                                            --functional-occurrence-table-output "${output_name}-functions-occurrence.txt"
        """
}

process computeANI {
    container "${params.container__anvio}"
    label "cpu_high"
    publishDir "${params.output_folder}", mode: "copy", overwrite: true
    
    input:
    path panGenome
    path combinedDB
    path genome_db_list
    path externalGenomes
    
    output:
    path "${panGenome}"
    path "ANI/*"

    """
#!/bin/bash

set -e

mkdir tmp
    
TMP=\$PWD/tmp \
TMPDIR=\$PWD/tmp \
anvi-compute-genome-similarity \
    --external-genomes ${externalGenomes} \
    --min-alignment-fraction ${params.min_alignment_fraction} \
    --output-dir ANI \
    --num-threads ${task.cpus} \
    --pan-db ${panGenome} \
    --program ${params.ani_program} \
    -T ${task.cpus}
    """
}

workflow {

    parseSampleSheet(
        file(params.sample_sheet, checkIfExists: true)
    )

    makeGenomeDB(
         parseSampleSheet
            .out
            .splitCsv(
                header:true,
                sep:"\t"
            ).map {
                row -> [
                    row.name,
                    file(row.genome, checkIfExists: true)
                ]
            }
    )

    if ( "${params.cogs_tar}" != "false" ){

        anvio_cogs_tar = file("${params.cogs_tar}", checkIfExists: true)

    } else {
        
        setupNCBIcogs()
        anvio_cogs_tar = setupNCBIcogs.out

    }

    annotateGenes(
        makeGenomeDB.out,
        anvio_cogs_tar
    )

    linkGeneName(
        makeGenomeDB.out
    )

    combineGenomes(
        annotateGenes.out.toSortedList(),
        linkGeneName.out.toSortedList()
    )

    panGenomeAnalysis(
        combineGenomes.out[0]
    )

    getSequencesForGCs(
        panGenomeAnalysis.out,
        combineGenomes.out[0]
    )

    addMetadata(
        panGenomeAnalysis.out,
        combineGenomes.out[0],
        parseSampleSheet.out
    )

    if ( params.category_name ){
        enrichFunctions(
            addMetadata.out,
            combineGenomes.out[0],
            Channel.of(params.category_name.split(","))
        )
    }

    computeANI(
        addMetadata.out,
        combineGenomes.out[0],
        makeGenomeDB
            .out
            .map { name, genome -> genome }
            .toSortedList(),
        combineGenomes.out[1]
    )
}