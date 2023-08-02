# Fix line break characters to facilitate grep etc
#sed 's/\^M/\n/g' biosample_result\ \(2\).txt > biosample_result_linebreak.txt
mv biosample_result \ \(2\).txt > biosample_raw.txt
# grep related rows and create a new text file
grep -E 'host sex|host age|dataset|health state|SRA:' biosample_raw.txt > output_raw.txt
#grep -E 'host sex|host age|dataset|health state|SRA:' biosample_result_linebreak.txt > output.txt

# Remove sample identifiers except for SRA: SRS[0-9]+
sed 's/.*\(SRA: SRS[0-9]\+\)/\1/' output_raw.txt > output_cleaner.txt

# Manipulate file so every 6 rows go into columns - paste - - - - - -: The paste command combines lines from input files horizontally, separating them by a delimiter (default is a tab). In this case, we use six hyphens ("-") to represent six input fields.
paste - - - - - - < output_cleaner.txt > output_cols.txt

echo -e "SRA\t host_age\t host_sex\t health_status\t dataset\t sequencing_pool" > header.txt
cat header.txt output_cols.txt > output_cols_header.txt

grep "dataset_1" PD_reads_metadata.csv | cut -f1 -d ',' > dataset_1_files.txt
sed -e 's/^"//' -e 's/"$//' < dataset_2_files.txt > dataset_2_SraAccList.txt
