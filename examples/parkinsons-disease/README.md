#  Parkinson's disease case study

An example dataset of 16S rRNA reads from 507 individuals was used to develop this workflow, including Parkinson's Disease (PD) subjects and control subjects. The data source for this project is [the Parkinson's Disease Bioproject](https://www.ncbi.nlm.nih.gov/bioproject/?term=601994) available on NCBI BioProject (ID: 601994).

## Preparing workflow inputs

Tool versions used:

- [`fasterq-dump`](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) v2.10.7
- [`R`](https://cran.r-project.org/) (tested on version 4.2.2)
- [`dplyr`](https://dplyr.tidyverse.org/) R library

### `fastq_tar`

The `fastq_tar` contains unzipped paired-end reads from all runs associated with the BioProject. To generate this file, the following steps were followed:

1. Obtain sample metadata from the NCBI

To obtain sample run metadata, [these instructions](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/#download-metadata-associated-wit) were followed using [this BioProject](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=601994). Biosample metadata was retrieved by selecting Full(text) from the dropbown menu, from [this link](https://www.ncbi.nlm.nih.gov/biosample?LinkName=bioproject_biosample&from_uid=601994).

The [biosample_result.txt](biosample_result.txt) found in this directory file was retrieved using this method.

2. Convert the project information to a TSV using the `connecting_metadata.sh` helper script.

```bash
./connecting_metadata.sh
```

This uses the `biosample_result.txt` file as input and outputs a `connecting_metadata_output.tsv` file.

3. Clean metadata and generate an accession list using the `PD_metadata_join.R` helper script.

```bash
Rscript PD_metadata_join.R
```

This uses the `connecting_metadata_output.tsv` file to generate the following files:

- `accession_list.txt`: a newline-delimited file with all SRA run IDs associated with the study.
- [`PD_metadata.csv`](PD_metadata.csv): a comma-separated file containing cleaned metadata from the BioProject. This file was used for downstream statistical analyses and data visualization in combination with the output from the Mothur workflow.

4. Download fastqs and tar FASTQ files

```bash
cat accessions_list.txt | xargs -I {} fasterq-dump {}

tar -cf PD_fastqs.tar --no-xattrs *.fastq
```

The `PD_fastqs.tar` file was used as the `fastq_tar` input to the workflow.

### `silva_ref_fasta`, `silva_ref_tax`

The SILVA ref fasta and tax files were retrieved as described in the [required inputs](../../README.md#required-inputs) section

### `oligos`

The V4 region of the 16S rRNA gene was targetted by this study. Since the reads stored in NCBI were already demultiplexed, barcodes were not provided. The contents of the oligos file was set to the following:

```txt
primer GTGCCAGCMGCCGCGGTAA GGACTACHVGGGTWTCTAAT V4
```

For more information on generating the oligos file, see [Mothur's documentation](https://mothur.org/wiki/oligos_file/).

### `ref_fasta_start`, `ref_fasta_end`

The parameters `ref_fasta_start` and `ref_fasta_end` were set to `11894` and `25319` respectively as per the [Mothur SOP](https://mothur.org/wiki/miseq_sop/) as the data used in this case study was sequenced at the V4 region of the 16S rRNA gene.

Instructions for customizing the start and end position of the reference alignment if a region other than the V4 region was sequenced are located [here](https://mothur.org/blog/2016/Customization-for-your-region/).

## Running the workflow

All default values from the workflow were accepted as they currently exist in the [Mothur workflow](../../mothur.wdl).
