#!/bin/bash

set -e
# Change "EXAMPLE_OUTPUT" to something that describes this group of genomes
OUTPUT_NAME=EXAMPLE_OUTPUT

NXF_VER=20.04.1 nextflow \
    -c ~/nextflow.config \
    run \
    main.nf \
    --sample_sheet tests/species_test.csv \
    --output_folder tests/output \
    --mcl_inflation 10 \
    --output_name $OUTPUT_NAME \
    --min_alignment_fraction 0 \
    --category_name species \
    -work-dir <INSERT WORK DIRECTORY> \
    --list-collections \
    -process.queue mixed \
    -resume

echo "Your genomes have been processed, and now the interactive anvi'o viewer is starting up. If this is the first time you have run this command, it may take ~5 minutes to download the anvi'o Docker image. Please wait..."

docker \
    run \
    -p 127.0.0.1:80:8080/tcp \
    -v $PWD:/share \
        meren/anvio:5.5 \
            anvi-display-pan \
            -g /share/tests/output/$OUTPUT_NAME-GENOMES.db \
            -p /share/tests/output/$OUTPUT_NAME-PAN.db
