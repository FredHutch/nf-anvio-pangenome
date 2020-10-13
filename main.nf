#!/usr/bin/env nextflow

params.output_name = "COMBINED_GENOMES"
params.output_folder = "./"

params.min_occurrence = "1"
params.minbit = "0.5"
params.distance = "euclidean"
params.linkage = "ward"
params.mcl_inflation = "2"
params.category_name = false
params.min_alignment_fraction = 0

process parseSampleSheet {
    container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
    label "io_limited"
    
    input:
    file sample_sheet_csv from file(params.sample_sheet)
    
    output:
    file "${sample_sheet_csv}" into sample_sheet_ch, sample_sheet_to_parse

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
    container "quay.io/fhcrc-microbiome/anvio:6.2.1"
    label "mem_medium"
    
    input:
    set name, file(fasta) from sample_sheet_ch.splitCsv(header:true, sep:"\t").map { row -> tuple(row.name, file(row.genome)) }
    
    output:
    set name, file("*db") into genomeDB_ch, nameDB_ch, aniDB_ch

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
    container "quay.io/fhcrc-microbiome/anvio:6.2.1"
    label "mem_medium"
    
    output:
    file "COGS_DIR.tar" into anvio_cogs_tar

    """
#!/bin/bash
anvi-setup-ncbi-cogs --num-threads 4 --cog-data-dir COGS_DIR --just-do-it --reset
tar cvf COGS_DIR.tar COGS_DIR
    """
}

process annotateGenes {
    container "quay.io/fhcrc-microbiome/anvio:6.2.1"
    label "mem_medium"
    
    input:
    set name, file(db) from genomeDB_ch
    file anvio_cogs_tar
    
    output:
    file "${db}" into annotatedDB

    """
#!/bin/bash
tar xvf ${anvio_cogs_tar}
anvi-run-ncbi-cogs -c "${db}" --num-threads 4 --cog-data-dir COGS_DIR
    """
}

process linkGeneName {
    container "quay.io/fhcrc-microbiome/anvio:6.2.1"
    label "mem_medium"
    
    input:
    set name, file(db) from nameDB_ch
    
    output:
    file "${db}.txt" into layer_txt_for_combineGenomes

    """
#!/bin/bash
# Link the name to the database
echo -e ${name},${db} | tr ',' '\\t' > ${db}.txt
    """
}

process combineGenomes {
    container "quay.io/fhcrc-microbiome/anvio:6.2.1"
    label "mem_medium"
    publishDir "${params.output_folder}"
    
    input:
    file db_list from annotatedDB.collect()
    file txt_list from layer_txt_for_combineGenomes.collect()
    
    output:
    file "${params.output_name}-GENOMES.db" into combinedDB
    file "external-genomes.txt" into external_genomes_for_ani

    """
#!/bin/bash
echo -e "name\\tcontigs_db_path" > external-genomes.txt
for fp in ${txt_list}; do cat \$fp; done >> external-genomes.txt
cat external-genomes.txt
anvi-gen-genomes-storage -e external-genomes.txt \
                         -o ${params.output_name}-GENOMES.db
    """
}

process panGenomeAnalysis {
    container "quay.io/fhcrc-microbiome/anvio:6.2.1"
    label "mem_medium"
    
    input:
    file combinedDB
    val output_name from params.output_name
    val min_occurrence from params.min_occurrence
    val minbit from params.minbit
    val distance from params.distance
    val linkage from params.linkage
    val mcl_inflation from params.mcl_inflation
    
    output:
    file "${params.output_name}-PAN.db" into panGenome_for_addMetadata

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
    container "quay.io/fhcrc-microbiome/anvio:6.2.1"
    label "io_limited"
    publishDir "${params.output_folder}"
    
    input:
    file panGenome from panGenome_for_addMetadata
    file combinedDB
    file sample_sheet from sample_sheet_to_parse
    
    output:
    file "${panGenome}" into panGenome_for_enrichFunctions, panGenome_for_ani

    """
#!/bin/bash

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

if ( params.category_name ){
    process enrichFunctions{
        container "quay.io/fhcrc-microbiome/anvio:6.2.1"
        label "io_limited"
        publishDir "${params.output_folder}"
        
        input:
        file panGenome from panGenome_for_enrichFunctions
        file combinedDB
        val output_name from params.output_name
        val category_name from params.category_name
        
        output:
        file "${output_name}-enriched-functions-${category_name}.txt"
        file "${output_name}-functions-occurrence.txt"

        script:

        if (params.gene_enrichment == false)
            """#!/bin/bash
            anvi-get-enriched-functions-per-pan-group -p ${panGenome} \
                                                    -g ${combinedDB} \
                                                    --category ${category_name} \
                                                    --annotation-source COG_FUNCTION \
                                                    -o "${output_name}-enriched-functions-${category_name}.txt" \
                                                    --functional-occurrence-table-output "${output_name}-functions-occurrence.txt"
            """

        else
            """#!/bin/bash
            anvi-get-enriched-functions-per-pan-group -p ${panGenome} \
                                                    -g ${combinedDB} \
                                                    --category ${category_name} \
                                                    --annotation-source IDENTITY \
                                                    --include-gc-identity-as-function \
                                                    -o "${output_name}-enriched-functions-${category_name}.txt" \
                                                    --functional-occurrence-table-output "${output_name}-functions-occurrence.txt"
            """
    }
}
    process computeANI {
        container "quay.io/fhcrc-microbiome/anvio:6.2.1"
        label "mem_veryhigh"
        publishDir "${params.output_folder}"
        
        input:
        file panGenome from panGenome_for_ani
        val output_name from params.output_name
        val min_alignment_fraction from params.min_alignment_fraction
        file combinedDB
        file genome_db_list from aniDB_ch.collect()
        file externalGenomes from external_genomes_for_ani
        
        output:
        file "${panGenome}"

        """
    #!/bin/bash

    set -e
        
    anvi-compute-genome-similarity \
        --external-genomes ${externalGenomes} \
        --min-alignment-fraction ${min_alignment_fraction} \
        --output-dir ANI \
        --num-threads 4 \
        --pan-db ${panGenome}
        """
    }