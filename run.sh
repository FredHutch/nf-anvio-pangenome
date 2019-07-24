#!/bin/bash

set -e
# Change "EXAMPLE_OUTPUT" to something that describes this group of genomes
OUTPUT_NAME=EXAMPLE_OUTPUT


NXF_VER=19.05.0-SNAPSHOT nextflow \
    -C ~/nextflow.config \
    run \
    main.nf \
    --sample_sheet tests/species_test.csv \
    --output_folder tests/output \
    --mcl_inflation 10 \
    --output_name $OUTPUT_NAME \
    --category_name species \
    -work-dir s3://fh-pi-fredricks-d/lab/Sam_Minot/data/nextflow/work/ \
    --list-collections \
    -process.queue mixed \
    -resume

docker \
    run \
    -p 127.0.0.1:80:8080/tcp \
    -v $PWD:/share \
        meren/anvio:5.5 \
            anvi-display-pan \
            -g /share/tests/output/$OUTPUT_NAME-GENOMES.db \
            -p /share/tests/output/$OUTPUT_NAME-PAN.db