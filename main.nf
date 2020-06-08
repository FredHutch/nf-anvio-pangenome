#!/usr/bin/env nextflow

params.output_name = "COMBINED_GENOMES"
params.output_folder = "./"

params.min_occurrence = "1"
params.minbit = "0.5"
params.distance = "euclidean"
params.linkage = "ward"
params.mcl_inflation = "2"

process parseSampleSheet {
    container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
    cpus 1
    memory "2 GB"
    
    input:
    file sample_sheet_csv from file(params.sample_sheet)
    
    output:
    file "${sample_sheet_csv}" into sample_sheet_ch

    afterScript "rm -rf *"

    """
#!/usr/bin/env python3
import pandas as pd

df = pd.read_csv("${sample_sheet_csv}", sep=",")

for k in ["name", "genome"]:
    assert k in df.columns.values, "Must provide a column '%s' in the sample sheet" % k

# Strip away all whitespace and carriage returns
df = df.apply(str).applymap(lambda s: s.strip())

df.to_csv("${sample_sheet_csv}", index=None, sep=",")
    """
}

process makeGenomeDB {
    container "meren/anvio:6.2"
    cpus 4
    memory "8 GB"
    
    input:
    set name, file(fasta) from sample_sheet_ch.splitCsv(header:true).map { row -> tuple(row.name, file(row.genome)) }
    
    output:
    set name, file("*db") into genomeDB_ch, nameDB_ch, aniDB_ch

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
    container "meren/anvio:6.2"
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
    container "meren/anvio:6.2"
    cpus 4
    memory "8 GB"
    
    input:
    set name, file(db) from genomeDB_ch
    file anvio_cogs_tar
    
    output:
    file "${db}" into annotatedDB

    afterScript "rm -rf *"

    """
#!/bin/bash
tar xvf ${anvio_cogs_tar}
anvi-run-ncbi-cogs -c "${db}" --num-threads 4 --cog-data-dir COGS_DIR
    """
}

process linkGeneName {
    container "meren/anvio:6.2"
    cpus 4
    memory "8 GB"
    
    input:
    set name, file(db) from nameDB_ch
    
    output:
    file "${db}.txt" into layer_txt_for_combineGenomes

    afterScript "rm -rf *"

    """
#!/bin/bash
# Link the name to the database
echo -e ${name},${db} | tr ',' '\\t' > ${db}.txt
    """
}

process combineGenomes {
    container "meren/anvio:6.2"
    cpus 4
    memory "8 GB"
    publishDir "${params.output_folder}"
    
    input:
    file db_list from annotatedDB.collect()
    file txt_list from layer_txt_for_combineGenomes.collect()
    
    output:
    file "${params.output_name}-GENOMES.db" into combinedDB
    file "external-genomes.txt" into external_genomes_for_ani

    afterScript "rm -rf *"

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
    container "meren/anvio:6.2"
    cpus 5
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
    file "${params.output_name}-PAN.db" into panGenome_for_addMetadata

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
    container "meren/anvio:6.2"
    cpus 4
    memory "2 GB"
    publishDir "${params.output_folder}"
    
    input:
    file panGenome from panGenome_for_addMetadata
    file combinedDB
    file sample_sheet from file("${params.sample_sheet}")
    
    output:
    file "${panGenome}" into panGenome_for_enrichFunctions, panGenome_for_ani


    afterScript "rm -rf *"

    """
#!/bin/bash
# Strip out the genome file path
cat ${sample_sheet} | tr '\\r' '\\n' | tr ',' '\\t' | cut -f 2- > TEMP && mv TEMP ${sample_sheet}
echo "Printing the reformatted sample sheet:"
cat ${sample_sheet}
echo ""
anvi-import-misc-data ${sample_sheet} \
                      -p ${panGenome} \
                      --target-data-table layers
    """
}

if ( params.category_name ){
    process enrichFunctions{
        container "meren/anvio:6.2"
        cpus 4
        memory "2 GB"
        publishDir "${params.output_folder}"
        
        input:
        file panGenome from panGenome_for_enrichFunctions
        file combinedDB
        val output_name from params.output_name
        val category_name from params.category_name
        
        output:
        file "${output_name}-enriched-functions-${category_name}.txt"
        file "${output_name}-functions-occurrence.txt"


        afterScript "rm -rf *"

        """
    #!/bin/bash
    anvi-get-enriched-functions-per-pan-group -p ${panGenome} \
                                              -g ${combinedDB} \
                                              --category ${category_name} \
                                              --annotation-source COG_FUNCTION \
                                              -o "${output_name}-enriched-functions-${category_name}.txt" \
                                              --functional-occurrence-table-output "${output_name}-functions-occurrence.txt"
        """
    }
}
    process computeANI {
        container "meren/anvio:6.2"
        cpus 4
        memory "16 GB"
        publishDir "${params.output_folder}"
        
        input:
        file panGenome from panGenome_for_ani
        val output_name from params.output_name
        val min_alignment_fraction from params.min_alignment_fraction
        file combinedDB
        file genome_db_list from aniDB_ch.collect()
        file externalGenomes from external_genomes_for_ani
        
        afterScript "rm -rf *"

        output:
        file "${panGenome}"

        """
    #!/bin/bash
        
        anvi-compute-ani --external-genomes ${externalGenomes} \
                     --min-alignment-fraction ${min_alignment_fraction} \
                     --output-dir ANI \
                     --num-threads 4 \
                     --pan-db ${panGenome}
        """
    }