#!/bin/bash

set -e

OUTPUT_NAME=EXAMPLE_OUTPUT

nextflow \
    -C ~/nextflow.config \
    run \
    main.nf \
    --sample_sheet tests/sample_sheet.txt \
    --output_directory tests/output \
    --output_name $OUTPUT_NAME \
    -process.queue mixed \
    -resume

docker \
    run \
    -p 127.0.0.1:80:8080/tcp \
    -v $PWD:/share \
        meren/anvio:5.5 \
            anvi-display-pan \
            -g /share/$OUTPUT_NAME-GENOMES.db \
            -p /share/$OUTPUT_NAME-PAN.db
