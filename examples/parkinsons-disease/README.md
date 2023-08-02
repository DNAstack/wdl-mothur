#  Parkinson's disease case study

An example dataset of 16S rRNA reads from 507 individuals was used to develop this workflow, including Parkinson's Disease (PD) subjects and control subjects. The data source for this project is [the Parkinson's Disease Bioproject](https://www.ncbi.nlm.nih.gov/bioproject/?term=601994) available on NCBI BioProject (ID: 601994).

## Example metadata

While sample metadata is not explicitly required for this workflow, metadata and the associated raw data and scripts used to generate PD-related metadata from the example dataset can be found in this directory in .csv format. Two raw metadata example files were obtained directly from NCBI (`SraRunInfoExample.csv` and `biosampleDataExample.txt`) and were manipulated using two scripts (`connecting_metadata.sh` and `PD_metadata_join.R`) to generate `PD_metadata.csv`. This .csv file was used for subsequent statistical analyses and data visualization.
