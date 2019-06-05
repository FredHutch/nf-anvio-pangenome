#!/bin/bash

set -e

nextflow \
    -C ~/nextflow.config \
    run \
    main.nf \
    --sample_sheet tests/sample_sheet.txt \
    --output_directory tests/output \
    -process.queue mixed \
    -resume
