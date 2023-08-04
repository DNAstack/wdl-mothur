# Load in dplyr and read in output file from connecting_metadata.sh
library(dplyr)
bio <- read.table('connecting_metadata_output.txt', sep = "\t",  header=TRUE)

# Use gsub to clean up output file from connecting_metadata.sh
bio$host_age <- gsub("\\D+", "", bio$host_age)
bio$host_sex<- gsub("    /host sex=", "", bio$host_sex)
bio$health_status<- gsub ("    /health state=", "", bio$health_status)
bio$dataset<- gsub ("    /dataset=", "", bio$dataset)
bio$sequencing_pool <- gsub ("    /sequencing_pool=", "", bio$sequencing_pool)
bio$SRA<- gsub (".*SRA: ", "", bio$SRA)

# Read in csv containing SRA run info
sra<- read.csv('SraRunInfo.csv')

# Rename sample to SRA to facilitate left_join
sra <- sra %>%
  rename (SRA = Sample)

# Join sra and bio data frames and filter out dataset_1 samples that were not used in this case study 
bio_sra<- left_join(bio, sra, by = 'SRA')%>%
  select (SRA, host_age, host_sex, health_status, dataset, sequencing_pool, Run, spots, bases, spots_with_mates, avgLength, size_MB, Experiment, LibraryName, BioSample, SampleName)%>%
  dplyr::filter (dataset == 'dataset_2')

# Write joined data frame to csv
write.csv(bio_sra, 'PD_metadata.csv', row.names = FALSE)

# Create and and write single column data frame of SRR identifiers to a text file for usage with fasterq-dump
accession_list <- as.data.frame(bio_sra$Run)
write.table(accession_list, 'accession_list.txt', col.names = FALSE, row.names=FALSE)
