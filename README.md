# nf-anvio-pangenome
Analyze a set of genomes with the [anvi'o](http://merenlab.org/software/anvio/) pangenome pipeline

The motivation here is to wrap up the set of anvi'o functions which would require
some familiarity of the command line, and isolate those within a reproducible
workflow which can be executed in a single command. The workflow can also take 
advanage of cloud computing in case the number of genomes is larger than your
computer can easily handle.

### Getting Started

The example below assumes the following:

  * [Docker Desktop](https://www.docker.com/products/docker-desktop) installed on your local computer
  * [Nextflow installed](https://nextflow.io) on your computer
  * [Nextflow configuration file](https://sciwiki.fredhutch.org/compdemos/nextflow/) at `~/nextflow.config`


### Input Data

In order to run this workflow, you need:

1. A set of genomes in FASTA format
2. A text file describing each of those genomes including a name, as well as additional optional metadata

The text file describing the genomes must be comma-delimited, at a minimum including 
`genome` and `name` (`genome` must be the first column).

NOTE: `name` must only contain alphanumeric characters (and `_`), and cannot start with a number.

```
genome,name,club,hobby
tests/genomes/GCF_000005845.2_ASM584v2_genomic.fna.gz,ASM584v2,shark,knitting
tests/genomes/GCF_000008865.2_ASM886v2_genomic.fna.gz,ASM886v2,shark,crochet
tests/genomes/GCF_000026325.1_ASM2632v1_genomic.fna.gz,ASM2632v1,shark,crochet
tests/genomes/GCF_000026345.1_ASM2634v1_genomic.fna.gz,ASM2634v1,jet,crochet
tests/genomes/GCF_000183345.1_ASM18334v1_genomic.fna.gz,ASM18334v1,jet,knitting
tests/genomes/GCF_000299455.1_ASM29945v1_genomic.fna.gz,ASM29945v1,jet,knitting
```

The additional columns to the right will be imported as metadata into the anvi'o 
pangenome viewer.

### Running the Workflow

To run the workflow, make a BASH script containing the following commands:

```
#!/bin/bash

set -e

# Get the most recent version of the workflow
nextflow pull fredhutch/nf-anvio-pangenome

# Run the workflow
nextflow \
    -C ~/nextflow.config \
    run \
    fredhutch/nf-anvio-pangenome \
    --sample_sheet sample_sheet.txt \
    --output_directory ./ \
    --output_name EXAMPLE_OUTPUT \
    -resume

```

### Visulizing the Pan-Genome

To launch the visual browser for the pan-genome, run the following command 
(replacing the database files and folder as appropriate):

```
docker \
    run \
    -p 127.0.0.1:80:8080/tcp \
    -v $PWD:/share \
        meren/anvio:5.5 \
            anvi-display-pan \
            -g /share/$OUTPUT_NAME-GENOMES.db \
            -p /share/$OUTPUT_NAME-PAN.db
```

For more details on how to navigate this visual browser, check out the amazing Anvi'o
[documentation](http://merenlab.org/2016/11/08/pangenomics-v2/.)

![Example Data](https://github.com/FredHutch/nf-anvio-pangenome/raw/master/assets/screenshot.png)
