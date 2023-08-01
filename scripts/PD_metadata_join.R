library(dplyr)
bio <- read.table('output_cols_header.txt', sep = "\t",  header=TRUE)
bio$host_age <- gsub("\\D+", "", bio$host_age)
bio$host_sex<- gsub("    /host sex=", "", bio$host_sex)
bio$health_status<- gsub ("    /health state=", "", bio$health_status)
bio$dataset<- gsub ("    /dataset=", "", bio$dataset)
bio$sequencing_pool <- gsub ("    /sequencing_pool=", "", bio$sequencing_pool)
bio$SRA<- gsub ("SRA: ", "", bio$SRA)

sra<- read.csv('SraRunInfo.csv')
sra <- sra %>%
  rename (SRA = Sample)

bio_sra<- left_join(bio, sra, by = 'SRA')%>%
  select (SRA, host_age, host_sex, health_status, dataset, sequencing_pool, Run, spots, bases, spots_with_mates, avgLength, size_MB, Experiment, LibraryName, BioSample, SampleName)

write.csv(bio_sra, 'PD_reads_metadata.csv', row.names = FALSE)
