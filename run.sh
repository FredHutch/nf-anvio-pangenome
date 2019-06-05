#!/bin/bash

set -e

nextflow \
    -C ~/nextflow.config \
    run \
    main.nf \
    --sample_sheet tests/sample_sheet.txt \
    --output_directory tests/output \
    -work-dir s3://fh-pi-fredricks-d/lab/Sam_Minot/data/nextflow/work/ \
    -process.queue mixed \
    -resume
