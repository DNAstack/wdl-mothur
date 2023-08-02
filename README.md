# Workflow: Process 16S rRNA reads with Mothur in Workflow Description Language

## Purpose
This repository provides instructions and a workflow in Workflow Description Language for processing 16S rRNA sequences using Mothur[1]. An example dataset of 16S rRNA reads from 507 individuals was used to develop this workflow, including Parkinson's Disease (PD) subjects and control subjects. The data source for this project is the Parkinson's Disease Bioproject available on NCBI BioProject (ID: 601994). The study associated with the data is titled "Characterizing dysbiosis of gut microbiome in PD: evidence for overabundance of opportunistic pathogens" [2].

## Example metadata
While sample metadata is not explicitly required for this workflow, metadata and the associated raw data and scripts used to generate PD-related metadata from the example dataset can be found in the [examples](examples/) subdirectory in .csv format. Two raw metadata example files were obtained directly from NCBI (`SraRunInfoExample.csv` and `biosampleDataExample.txt`) and were manipulated using two scripts (`connecting_metadata.sh` and `PD_metadata_join.R`) to generate `PD_metadata.csv`. This .csv file was used for subsequent statistical analyses and data visualization. 

## Workflow inputs
| Type | Input | Description |
|:-|:-|:-|
| File | `fastq_tar` | Compressed tarball of multiple paired-end fastq files
| String | `container_registry` | Registry that hosts workflow containers
| File | `silva_ref_fasta` | Full length sequence  database from SILVA (https://www.arb-silva.de/)
| File | `silva_ref_tax` | Full length reference taxonomy database from SILVA (https://www.arb-silva.de/)
| File | `oligos` | The oligos file provides barcodes and primers to Mothur
| String | `file_type` | Input file type (either fastq or gz)
| String | `prefix` | String of prefix for output files
| Int | `processors` | Number of processors to use in various tasks
| Int | `max_ambig` | The number of ambiguous bases permitted in the make.contigs command
| Int | `max_length` | The maximum read length permitted in the make.contigs command
| Int | `max_homop` | The maximum number of homopolymers permitted in the make.contigs command
| Int | `primer_differences` | The maximum number of differences to the primer sequence allowed in the make.contigs command
| Int | `ref_fasta_start` | Start position to trim reference fasta from in pcr.seqs command
| Int | `ref_fasta_end` | End position to trim reference fasta from in pcr.seqs command
| Int | `sequence_differences` | Maximum number of sequence differences that can be clustered together in the pre.cluster command
| String | `classify_method` | Sequence classification method in classify.seqs command. The options are: 'wang' or 'knn'.
| String | `taxon_remove` | Taxa to remove from sequence data in remove.lineage command
| String | `align_method` | Alignment method to use in make.contigs command. The options are: 'needleman' or 'gotoh'.
| String | `search_method` | Search method to use in align.seqs command. The options are: 'kmer' or 'suffix'.
| Int | `ksize` | Size of kmers to use in align.seqs command
| Boolean | `keepdots` | Parameter to keep leading and trailing .'s in reference sequence database in pcr.seqs command. Default=true.
| Int | `match` | Reward for a match in the pairwise alignment portion of the aligning procedure of align.seqs command
| Int | `mismatch` | Penalty for a mismatch in the pairwise alignment portion of the aligning procedure of align.seqs command
| Int | `gap_open` | Penalty for a gap open in the pairwise alignment portion of the aligning procedure of align.seqs command
| Int | `gap_extend` | Penalty for a gap extend in the pairwise alignment portion of the aligning procedure of align.seqs command
| Boolean | `flip` | Parameter used to specify whether or not Mothur should try the reverse complement of a sequence if the sequence falls below the threshold in align.seqs command. Default=true.
| Float | `threshold` | Parameter used to specify a cutoff below which an alignment is deemed poor and the reverse complement may be tried in align.seqs command
| Boolean | `vertical` | Parameter used to ignore/include any column that only contains gap characters in filter.seqs command. Default=true.
| String | `trump` | Parameter to remove a column if the trump character is found at that position in any sequence of the alignment in filter.seqs command
| String | `precluster_method` | Algorithm used to perform preclustering in pre.cluster command. The options are: 'simple', 'tree', 'unoise', and 'deblur'.
| String | `classify_output_format` | Format of *.tax.summary file in classify.seqs command. The options are: 'simple' and 'detail'.
| Int | `classify_seqs_cutoff` | Bootstrap value below which to remove sequences in classify.seqs command
| Boolean | `relative_abundance` | Value to indicate if summary file values should be displayed as relative abundances rather than raw abundances in classify.seqs command. Default=false.
| Boolean | `bootstrap_probabilities` | Parameter to include/exclude bootstrap probabilities from output in classify.seqs command. Default=true.
| Int | `classify_seqs_iters` | Number of iterations to perform when calculating the bootstrap confidence score for the taxonomy in classify.seqs command
| Boolean | `dereplicate_chimeras` | The dereplicate parameter can be used when checking for chimeras by group in remove.chimeras command. If the dereplicate parameter is false, then if one group finds the sequence to be chimeric, then all groups find it to be chimeric, default=f. If dereplicate=t, and then when a sequence is found to be chimeric it is removed from it’s group, not the entire dataset.
| Boolean | `remove_chimeras` | Parameter to remove the chimeras instead of just flagging them in remove.chimeras command. Default=true.
| Float | `minh_chimeras` | Mininum score to report chimera in remove.chimeras command
| Float | `mindiv_chimeras` | Minimum divergence ratio in remove.chimeras command
| Int | `mindiffs_chimeras` | Minimum number of differences in segment in remove.chimeras command
| Float | `xn_chimeras` | Weight of a no vote in remove.chimeras command
| Float | `dn_chimeras` | pseudo-count prior on number of no votes in remove.chimeras command
| Float | `cluster_cutoff` | Similarity threshold used for clustering sequences in cluster command
| String | `cluster_method` | Method for clustering sequences based on similarity in cluster command. The options are: 'opti', 'average', 'furthest', 'nearest', 'agc' 'dgc', and 'unique'.
| String | `cluster_metric` | Metric in the opticluster method in cluster command. The options are: 'mcc', 'sens', 'spec', 'tptn', 'fpfn', 'tp', 'tn', 'fp', 'fn', 'f1score', 'accuracy', 'ppv', 'npv', and 'fdr'.
| String | `cluster_initialize` | Initial randomization for the opticluster method in cluster command. The options are: 'singleton' and 'oneotu'.
| Float | `cluster_delta` | Stable value for the metric in the opticluster method in cluster command
| Int | `cluster_iters` | Maxiters for the opticluster method in cluster command
| Int | `classify_otus_cutoff` | Consensus confidence threshold for taxonomy in classify.otus command
| String | `classify_basis` | The basis parameter indicates what the summary file should represent in classify.otus command. The options are: 'sequence' and 'otu'.
| Int | `classify_printlevel` | Taxlevel of the *tax.summary file to print to in classify.otus command. The options are 1 to the max level in the final .taxonomy file; a value of -1 will print the max level within the .taxonomy file.
| Boolean | `classify_probabilities` | Parameter that shuts off the outputting of the consensus confidence results in classify.otus command. Default=true, meaning the confidence will be shown.
| Boolean | `classify_persample` | Parameter that indicates the assignment of a consensus taxonomy for each group in classify.otus command. Default=false.
| Int | `classify_threshold` | Cutoff for the taxonomy input file in classify.otus command. It’s a way to after the fact “adjust” the cutoff used in classify.seqs command without having to reclassify
| Boolean | `check_orient` | If you are running the make.contigs command with paired barcodes or primers, the checkorient parameter can be used. When checkorient=t and Mothur can’t find the barcodes and primers, it will search the reverse compliment. Default=true.
| String | `sequencing_format` | Sequencing format in make.contigs command. The options are: 'sanger', 'solexa', 'illumina1.8+' and 'illumina'.
| Int | `quality_score_threshold` | Quality score threshold in make.contigs command
| Float | `dist_cutoff` | Cutoff in dist.seqs command; if you know that you are not going to form OTUs with distances larger than 0.10, Mothur can be instructed to not save any distances larger than 0.10. This will significantly cut down on the amount of hard drive space required to store the matrix.
| Boolean | `count_ends` | Parameter to penalize gaps that occur at end of sequences in dist.seqs command. Default=true.
| String | `dist_calc` | Distance calculation method in dist.seqs command. The options are: 'onegap', 'nogaps', and 'eachgap'.
| Int | `tax_level` | Integer parameter representing the level of taxonomic classification to classify up to; a value of -1 will print the max level within the .taxonomy file.

## Workflow outputs
| Type | Output | Description |
|:-|:-|:-|
| File | `make_contigs_fasta` | .fasta file with sequences from make_contigs_scrap_fasta removed during make_contigs task
| File | `make_contigs_scrap_fasta` | Sequences that were removed from raw reads in make_contigs task       
| File | `align_seqs_summary` | Summary of sequence alignment output by align_seqs_to_silva
| File | `classify_seqs_taxonomy` | Sequence and corresponding taxonomic classification prior to the removal of unwanted lineages
| File | `remove_lineage_fasta` | .fasta file to be have singletons and doubletons removed and then to be used in subsequent ASV/OTU generation
| File | `remove_lineage_taxonomy` | Sequence and corresponding taxonomic classification with lineages removed as per the taxon_remove input parameter
| File | `fasta_singletons_and_doubletons_removed` | Final .fasta file to be used in subsequent ASV/OTU generation
| File | `count_table_singletons_and_doubletons_removed` | Final .count_table file to be used in subsequent ASV/OTU generation
| Array[File] | `OTUs` | Array of output files that correspond to OTU data to be used for downstream processing/analysis/visualization
| Array[File] | `ASVs` | Array of output files that correspond to ASV data to be used for downstream processing/analysis/visualization


## References
[1] Schloss, P. D., Westcott, S. L., Ryabin, T., Hall, J. R., Hartmann, M., Hollister, E. B., Lesniewski, R. A., Oakley, B. B., Parks, D. H., Robinson, C. J., Sahl, J. W., Stres, B., Thallinger, G. G., Van Horn, D. J., & Weber, C. F. (2009). Introducing mothur: Open-Source, Platform-Independent, Community-Supported Software for Describing and Comparing Microbial Communities. Applied and Environmental Microbiology, 75(23), 7537–7541. https://doi.org/10.1128/AEM.01541-09

[2] Wallen, Z.D., Appah, M., Dean, M.N. et al. Characterizing dysbiosis of gut microbiome in PD: evidence for overabundance of opportunistic pathogens. npj 
Parkinsons Dis. 6, 11 (2020). https://doi.org/10.1038/s41531-020-0112-6
