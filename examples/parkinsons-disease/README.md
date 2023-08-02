#  Parkinson's disease case study

An example dataset of 16S rRNA reads from 507 individuals was used to develop this workflow, including Parkinson's Disease (PD) subjects and control subjects. The data source for this project is [the Parkinson's Disease Bioproject](https://www.ncbi.nlm.nih.gov/bioproject/?term=601994) available on NCBI BioProject (ID: 601994).

## Example metadata

While sample metadata is not explicitly required for this workflow, metadata and the associated raw data and scripts used to generate PD-related metadata from the example dataset can be found in this directory in .csv format. Two raw metadata files were obtained directly from NCBI (`SraRunInfo.csv` and `biosample_raw.txt`) and were manipulated using two scripts (`connecting_metadata.sh` and `PD_metadata_join.R`) to generate `PD_metadata.csv`. This .csv file was used for subsequent statistical analyses and data visualization.

The `connecting_metadata.sh` script can be called from the current directory using `./connecting_metadata.sh`. This will use the `biosample_result.txt` file as input and output a `connecting_metadata_output.txt` file that can then be imported into R and run through the `PD_metadata_join.R` file, which will produce a `PD_metadata.csv`. The `PD_metadata_join.R` file will also produce a file with a list of accession numbers (`accession_list.txt`) that is required to download fastq files using fasterq-dump. The full metadata file used in this case study is also [provided](PD_metadata.csv). The metadata .csv file was used solely for statistical analyses and data visualization that was run on the data output by the Mothur workflow written in Workflow Description Language.

## Preparing inputs for workflow
To download raw fastq files, fasterq-dump v2.10.7 was used along with the `accessions_list.txt` file output from `PD_metadata_join.R` using the command below. <br>
```cat accessions_list.txt | xargs -I {} fasterq-dump {}```

To tar the fastq files, the command below was used. <br>
```tar -cf example.tar --no-xattrs *.fastq```

The parameters `ref_fasta_start` and `ref_fasta_end` were set to 11894 and 25319 respectively as per the [Mothur SOP](https://mothur.org/wiki/miseq_sop/) as the data used in this case study was sequenced at the V4 region of the 16S rRNA gene. Instructions for customizing the start and end position of the reference alignment if a region other than the V4 region was sequenced are located [here](https://mothur.org/blog/2016/Customization-for-your-region/).

## Running the workflow
All default values from the workflow were accepted as they currently exist in the [Mothur workflow](../../mothur.wdl).