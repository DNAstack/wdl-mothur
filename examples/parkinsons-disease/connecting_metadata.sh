#!/bin/bash

# grep related rows and create a new text file
grep -E 'host sex|host age|dataset|health state|SRA:' biosampleDataExample.txt > output_grep.txt

# Manipulate file so every 6 rows go into columns - paste - - - - - -: The paste command combines lines from input files horizontally, separating them by a delimiter (default is a tab). In this case, we use six hyphens ("-") to represent six input fields.
paste - - - - - - < output_grep.txt > output_cols.txt

# Manually create a file with column names
echo -e "SRA\t host_age\t host_sex\t health_status\t dataset\t sequencing_pool" > header.txt

# Combine our header and output file to create a text file that can be imported into PD_metadata_join.R for further manipulation
cat header.txt output_cols.txt > output_example.txt

# Remove intermediate files
rm header.txt output_cols.txt output_grep.txt