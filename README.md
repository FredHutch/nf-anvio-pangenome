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
            `genome` must be the complete path to the file not only `tests/genomes`. For example my
            file path is  `/Users/willfrohlich/Desktop/Hutch/G.Vaginalis/tests/genomes/`.

```
genome,name,species
tests/genomes/GCA_004336685.1_ASM433668v1_genomic.fna,ATCC_14018,Gardnerella vaginalis
tests/genomes/GCA_003397665.1_ASM339766v1_genomic.fna,Ugent_09_07,Gardnerella vaginalis
tests/genomes/GCA_003397605.1_ASM339760v1_genomic.fna,Ugent_25_49,Gardnerella vaginalis
tests/genomes/GCA_003397755.1_ASM339775v1_genomic.fna,Ugent_09_01,Gardnerella vaginalis
tests/genomes/GCA_003293675.1_ASM329367v1_genomic.fna,UGent_06_41,Gardnerella leopoldii
tests/genomes/GCA_003397635.1_ASM339763v1_genomic.fna,Ugent_09_48,Gardnerella leopoldii
tests/genomes/GCA_003397585.1_ASM339758v1_genomic.fna,Ugent_18_01,Gardnerella piotii
tests/genomes/GCA_003397615.1_ASM339761v1_genomic.fna,Ugent_21_28,Gardnerella piotii
tests/genomes/GCA_003397705.1_ASM339770v1_genomic.fna,GS_9838_1,Gardnerella swidsinskii
tests/genomes/GCA_003397745.1_ASM339774v1_genomic.fna,GS_10234,Gardnerella swidsinskii
```

The additional columns to the right will be imported as metadata into the anvi'o 
pangenome viewer.

### Running the Workflow

To run the workflow, make a BASH script containing the following commands:

```
#!/bin/bash

set -e

# Change "EXAMPLE_OUTPUT" to something that describes this group of genomes
OUTPUT_NAME=EXAMPLE_OUTPUT

# Get the most recent version of the workflow
nextflow pull fredhutch/nf-anvio-pangenome

# Run the workflow
nextflow \
    -C ~/nextflow.config \
    run \
    fredhutch/nf-anvio-pangenome \
    --sample_sheet species_test.csv \
    --output_folder ./ \
    --output_name $OUTPUT_NAME \
    --min_alignment_fraction 0 \
    --category_name species \
    --mcl_inflation 10 \
    -work-dir #insert your work dir
    -process.queue mixed \
    -resume

```

### Picking Parameters

You can use the `--mcl_inflation` parameter to control how similar two genes have to be
for grouping together. The Meren lab recommends using a value of 10 for comparing closely
related genomes (for example, strains from the same species), and a value of 2 for comparing
distantly related genomes (for example, across species).

You can use the `--category_name` parameter to select which column of metadata you want 
to use to identify groups within your genomes and find functions that are enriched in those 
groups: functions that are characteristic of these genomes, and predominantly absent from 
genomes from outside this group. This data will be available in your output folder and be 
titled YOUR_PANGENOME-enriched-functions-category.txt

You can use the `--min_allignment_fraction` parameter to eliminate ANI scores between two 
genomes if the alignment fraction is less than you deem trustable. When you set a value, anvi'o 
will go through the ANI results, and set percent identity scores between two genomes to 0 if
alignment fraction *between either of them* is less than the parameter described here. The default
is 0.0, so every hit is reported, but you can choose any value between 0.0 and 1.0. 


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
            -g /share/tests/output/$OUTPUT_NAME-GENOMES.db \
            -p /share/tests/output/$OUTPUT_NAME-PAN.db
```

For more details on how to setup and navigate this visual browser, check out this [presentation](https://figshare.com/articles/Fred_Hutch_Nextflow_Anvi_o_Pangenome_Pipeline_Bacterial_Vaginosis/9273533)

Or, check out the amazing Anvi'o
[documentation](http://merenlab.org/2016/11/08/pangenomics-v2/.)

![Example Data](https://github.com/FredHutch/nf-anvio-pangenome/raw/master/assets/screenshot.png)
