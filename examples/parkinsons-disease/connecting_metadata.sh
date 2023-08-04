#!/bin/bash

# grep related rows and create a new text file
grep -E 'host sex|host age|dataset|health state|SRA:' biosample_result.txt > output_grep.txt

# Manipulate file so every 6 rows go into columns - paste - - - - - -: The paste command combines lines from input files horizontally, separating them by a delimiter (default is a tab). In this case, we use six hyphens ("-") to represent six input fields.
paste - - - - - - < output_grep.txt > output_cols.tsv

# Manually create a file with column names
echo -e "SRA\t host_age\t host_sex\t health_status\t dataset\t sequencing_pool" > header.tsv

# Combine our header and output file to create a text file that can be imported into PD_metadata_join.R for further manipulation
cat header.tsv output_cols.tsv > connecting_metadata_output.tsv

# Create a line_count variable to check the size of our output metadata file (should be 1 line per individual in addition to the header)
line_count=$(wc -l < connecting_metadata_output.tsv)

# Remove intermediate files
rm header.tsv output_cols.tsv output_grep.txt

# Check the exit status and number of lines in output file and print the appropriate message
if [ $? -eq 0 ]; then
  line_count=$(wc -l < connecting_metadata_output.tsv)
  
  if [ $line_count -gt 1 ]; then
    printf "connecting_metada.sh complete!\nOutput file has $line_count lines and can be imported into PD_metadata_join.R.\n"
  else
    printf "Error: connecting_metada.sh completed successfully, but the output file has only $line_count line. Please check the script's output.\n"
  fi
else
  printf "Error: connecting_metada.sh encountered an issue during execution.\n"
fi