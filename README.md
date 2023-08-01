# Workflow: Process 16S rRNA reads with mothur

## Purpose
This repository provides instructions and a workflow for processing 16S rRNA sequences from 507 individuals, including Parkinson's Disease (PD) patients and control samples using mothur[1]. The data source for this project is the Parkinson's Disease Bioproject available on NCBI BioProject (ID: 601994). The study associated with the data is titled "Characterizing dysbiosis of gut microbiome in PD: evidence for overabundance of opportunistic pathogens" [2].

## Metadata
While sample data metadata is not explicitly required for this workflow, metadata can be found here (gs://jesse_dev/PD_metadata) across 3 files. Two files are obtained directly from NCBI and were manipulated using scripts found in the scripts subdirectory to generate the PD_reads_metadata.csv file. This .csv file was used for subsequent statistical analyses and data visualization. Examples of the two metadata files from NCBI are also present in the scripts subdirectory within this repository.

## Workflow inputs
| Type | Input | Description |
|:-|:-|:-|
| File | `fastq_tar` | Compressed tarball of multiple paired-end fastq files
| String | `container_registry` | Registry that hosts workflow containers
| File | `silva_ref_fasta` | Full length sequence  database from SILVA (https://www.arb-silva.de/)
| File | `silva_ref_tax` | Full length reference taxonomy database from SILVA (https://www.arb-silva.de/)
| File | `oligos` | The oligos file is used to provide barcodes and primers to mothur
| String | `file_type` | Input file type (either fastq or gz)
| String | `prefix` | String of prefix for output files
| Int | `processors` | Number of processors to use in various tasks
| Int | `max_ambig` | The number of ambiguous bases allowed in the make.contigs command
| Int | `max_length` | The maximum read length allowed in the make.contigs command
| Int | `max_homop` | The maximum number of homopolymers allowed in the make.contigs command
| Int | `primer_differences` | The maximum number of differences to the primer sequence allowed in the make.contigs command
| Int | `ref_fasta_start` | Start position to trim reference fasta from in pcr.seqs command
| Int | `ref_fasta_end` | End position to trim reference fasta from in pcr.seqs command
| Int | `sequence_differences` | Maximum number of sequence differences that can be clustered together in the pre.cluster command
| String | `classify_method` | Sequence classification method in classify.seqs command
| String | `taxon_remove` | Taxa to remove from sequence data in remove.lineage command
| String | `align_method` | Alignment method to use in make.contigs command
| String | `search_method` | Search method to use in align.seqs command
| Int | `ksize` | Size of kmers to use in align.seqs command
| Boolean | `keepdots` | Parameter to keep leading and trailing .'s in reference sequence database in pcr.seqs command
| Int | `match` | Reward for a match in pairwise alignment portion of the aligning procedure of align.seqs command
| Int | `mismatch` | Penalty for a mismatch in pairwise alignment portion of the aligning procedure of align.seqs command
| Int | `gap_open` | Penalty for a gap open in pairwise alignment portion of the aligning procedure of align.seqs command
| Int | `gap_extend` | Penalty for a gap extend in pairwise alignment portion of the aligning procedure of align.seqs command
| Boolean | `flip` | Parameter used to specify whether or not you want mothur to try the reverse complement of a sequence if the sequence falls below the threshold in align.seqs command
| Float | `threshold` | Parameter used to specify a cutoff at which an alignment is deemed ‘bad’ and the reverse complement may be tried in align.seqs command
| Boolean | `vertical` | Boolean parameter to ignore/include any column that only contains gap characters in filter.seqs command
| String | `trump` | Parameter to remove a column if the trump character is found at that position in any sequence of the alignment in filter.seqs command
| String | `precluster_method` | Algorithm used to perform preclustering in pre.cluster command
| String | `classify_output_format` | Format of *.tax.summary file in classify.seqs command
| Int | `classify_seqs_cutoff` | Bootstrap value below which to remove sequences in classify.seqs command
| Boolean | `relative_abundance` | Allows you to indicate if you want the summary file values to be relative abundances rather than raw abundances in classify.seqs command
| Boolean | `bootstrap_probabilities` | Parameter to include/exclude bootstrap probabilities from output in classify.seqs command
| Int | `classify_seqs_iters` | Number of iterations to do when calculating the bootstrap confidence score for your taxonomy in classify.seqs command
| Boolean | `dereplicate_chimeras` | The dereplicate parameter can be used when checking for chimeras by group in remove.chimeras command. If the dereplicate parameter is false, then if one group finds the sequence to be chimeric, then all groups find it to be chimeric, default=f. If you set dereplicate=t, and then when a sequence is found to be chimeric it is removed from it’s group, not the entire dataset.
| Boolean | `remove_chimeras` | Parameter to remove the chimeras from your files instead of just flagging them in remove.chimeras command
| Float | `minh_chimeras` | Mininum score to report chimera in remove.chimeras command
| Float | `mindiv_chimeras` | Minimum divergence ratio in remove.chimeras command
| Int | `mindiffs_chimeras` | Minimum number of differences in segment in remove.chimeras command
| Float | `xn_chimeras` | Weight of a no vote in remove.chimeras command
| Float | `dn_chimeras` | pseudo-count prior on number of no votes in remove.chimeras command
| Float | `cluster_cutoff` | Similarity threshold used for clustering sequences in cluster command
| String | `cluster_method` | Method for clustering sequences based on similarity in cluster command
| String | `cluster_metric` | Metric in the opticluster method in cluster command
| String | `cluster_initialize` | Initial randomization for the opticluster method in cluster command
| Float | `cluster_delta` | Stable value for the metric in the opticluster method in cluster command
| Int | `cluster_iters` | Maxiters for the opticluster method in cluster command
| Int | `classify_otus_cutoff` | Consensus confidence threshold for your taxonomy in classify.otus command
| String | `classify_basis` | The basis parameter allows you indicate what you want the summary file to represent in classify.otus command, options are otu and sequence in classify.otus command
| Int | `classify_printlevel` | Taxlevel of your *tax.summary file to print to in classify.otus command
| Boolean | `classify_probabilities` | Parameter that shuts off the outputting of the consensus confidence results in classify.otus command. The default is true, meaning you want the confidence to be shown.
| Boolean | `classify_persample` | Parameter that allows you to find a consensus taxonomy for each group in classify.otus command. Default=false.
| Int | `classify_threshold` | Cutoff for your taxonomy input file in classify.otus command. It’s a way to after the fact “adjust” the cutoff used in classify.seqs command without having to reclassify
| Boolean | `check_orient` | If you are running the make.contigs command with paired barcodes or primers, you can use the checkorient parameter. When checkorient=t and mothur can’t find the barcodes and primers, it will search the reverse compliment. Default=true.
| String | `sequencing_format` | Sequencing format in make.contigs command
| Int | `quality_score_threshold` | Quality score threshold in make.contigs command
| Float | `dist_cutoff` | Cutoff in dist.seqs command; if you know that you are not going to form OTUs with distances larger than 0.10, you can tell mothur to not save any distances larger than 0.10. This will significantly cut down on the amount of hard drive space required to store the matrix.
| Boolean | `count_ends` | Parameter to penalize gaps that occur at end of sequences in dist.seqs command
| String | `dist_calc` | Distance calculation method in dist.seqs command
| Int | `tax_level` | Integer parameter representing the level of taxonomic classification to classify up to (-1 prints the maxlevel)

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
