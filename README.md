# Workflow: Process 16S rRNA reads with Mothur in Workflow Description Language

## Purpose

This repository provides instructions and a workflow in Workflow Description Language for processing 16S rRNA sequences using Mothur [[1](#references)].

The Docker images used in this workflow are defined in DNAstack's [bioinformatics image repository](https://github.com/DNAstack/bioinformatics-public-docker-images/tree/main), and may be found and pulled from [DNAstack's DockerHub container registry](https://hub.docker.com/u/dnastack).

## Parkinson's Disease case study

An example dataset of 16S rRNA reads from 507 individuals was used to develop this workflow, including Parkinson's Disease (PD) subjects and control subjects. The data source for this project is [the Parkinson's Disease Bioproject](https://www.ncbi.nlm.nih.gov/bioproject/?term=601994) available on NCBI BioProject (ID: 601994). The study associated with the data is titled "Characterizing dysbiosis of gut microbiome in PD: evidence for overabundance of opportunistic pathogens" [[2](#references)].

Example metadata and scripts for processing the output of the PD case study can be found in the [examples/parkinsons-disease](examples/parkinsons-disease) directory.

## Workflow inputs

A template input file containing only required inputs can be found [here](input_template.json).

### Required inputs

| Type | Input | Description |
| :- | :- | :- |
| File | `fastq_tar` | Compressed tarball of multiple paired-end fastq files. |
| File | `silva_ref_fasta` | Full length sequence database from [SILVA](https://www.arb-silva.de/). |
| File | `silva_ref_tax` | Full length reference taxonomy database from [SILVA](https://www.arb-silva.de/). <br> The silva_ref_fasta and silva_ref_tax files were downloaded using the following command:<br> `wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138_1.tgz && tar -xzf silva.nr_v138_1.tgz && rm silva.nr_v138_1.tgz`. <br>The resulting .align and .tax files were not manipulated prior to being used in this workflow. |
| File | `oligos` | The oligos file provides barcodes and primers to Mothur. See Mothur's [documentation](https://mothur.org/wiki/oligos_file/) for information on formatting this file. 
| String | `prefix` | Prefix used for naming output files; this can be any string that the users wishes to appear prior to the rest of the file name(s). |
| Int | `ref_fasta_start` | Start position to trim reference fasta from in pcr.seqs command. This value depends on the region of the 16S rRNA gene that was sequenced. Instructions on how to customize the reference alignment based on the region sequenced can be found [here](https://mothur.org/blog/2016/Customization-for-your-region/) |
| Int | `ref_fasta_end` | End position to trim reference fasta from in pcr.seqs command. This value depends on the region of the 16S rRNA gene that was sequenced. Instructions on how to customize the reference alignment based on the region sequenced can be found [here](https://mothur.org/blog/2016/Customization-for-your-region/) |

### Optional inputs

The values defined here have defaults set that may be overridden by the user.

| Type | Input | Description |
| :- | :- | :- |
| Int | `max_ambig` | The number of ambiguous bases permitted in the make.contigs command. [0] |
| Int | `max_homop` | The maximum number of homopolymers permitted in the make.contigs command. [8] |
| Int | `max_length` | The maximum read length permitted in the make.contigs command. [300] |
| Int | `primer_differences` | The maximum number of differences to the primer sequence allowed in the make.contigs command. [4] |
| Boolean | `check_orient` | If you are running the make.contigs command with paired barcodes or primers, the `checkorient` parameter can be used. When `checkorient=t` and Mothur can’t find the barcodes and primers, it will search the reverse compliment. [true] |
| String | `sequencing_format` | Sequencing format in make.contigs command. ('sanger', 'solexa', 'illumina1.8+', 'illumina') ['illumina1.8+'] |
| Int | `quality_score_threshold` | Quality score threshold in make.contigs command. [20] |
| Boolean | `keepdots` | Parameter to keep leading and trailing `.`s in reference sequence database in pcr.seqs command. [false] |
| String | `align_method` | Alignment method to use in make.contigs command. ('needleman', 'gotoh') ['needleman'] |
| Int | `match` | Reward for a match in the pairwise alignment portion of the aligning procedure of align.seqs command. [1] |
| Int | `mismatch` | Penalty for a mismatch in the pairwise alignment portion of the aligning procedure of align.seqs command. [-1] |
| Int | `gap_open` | Penalty for a gap open in the pairwise alignment portion of the aligning procedure of align.seqs command. [-2] |
| Int | `gap_extend` | Penalty for a gap extend in the pairwise alignment portion of the aligning procedure of align.seqs command. [-1] |
| Int | `ksize` | Size of kmers to use in align.seqs command. [8] |
| String | `search_method` | Search method to use in align.seqs command. ('kmer', 'suffix') ['kmer'] |
| Boolean | `flip` | Parameter used to specify whether or not Mothur should try the reverse complement of a sequence if the sequence falls below the threshold in align.seqs command. [true] |
| Float | `threshold` | Parameter used to specify a cutoff below which an alignment is deemed poor and the reverse complement may be tried in align.seqs command. [0.5] |
| String | `trump_character` | If the trump character is found in any position in any sequence of the alignment, remove that column in the filter.seqs command. ['.'] |
| Boolean | `vertical` | Parameter used to ignore/include any column that only contains gap characters in filter.seqs command. [true] |
| Int | `sequence_differences` | Maximum number of sequence differences that can be clustered together in the pre.cluster command. [2] |
| String | `precluster_method` | Algorithm used to perform preclustering in pre.cluster command. ('simple', 'tree', 'unoise', 'deblur') ['simple'] |
| Boolean | `remove_chimeras` | Parameter to remove the chimeras instead of just flagging them in remove.chimeras command. [true] |
| Boolean | `dereplicate_chimeras` | The dereplicate parameter can be used when checking for chimeras by group in remove.chimeras command. If the dereplicate parameter is false, then if one group finds the sequence to be chimeric, then all groups find it to be chimeric, default=f. If dereplicate=t, and then when a sequence is found to be chimeric it is removed from it’s group, not the entire dataset. [true] |
| Float | `minh_chimeras` | Mininum score to report chimera in remove.chimeras command. [0.3] |
| Float | `mindiv_chimeras` | Minimum divergence ratio in remove.chimeras command. [0.5] |
| Int | `mindiffs_chimeras` | Minimum number of differences in segment in remove.chimeras command. [3] |
| Float | `xn_chimeras` | Weight of a no vote in remove.chimeras command. [8.0] |
| Float | `dn_chimeras` | pseudo-count prior on number of no votes in remove.chimeras command. [1.4] |
| String | `classify_method` | Sequence classification method in classify.seqs command. ('wang' or 'knn') ['wang'] |
| String | `classify_output_format` | Format of *.tax.summary file in classify.seqs command. ('simple', 'detail') ['simple'] |
| Int | `classify_seqs_cutoff` | Bootstrap value below which to remove sequences in classify.seqs command. [80] |
| Boolean | `bootstrap_probabilities` | Parameter to include/exclude bootstrap probabilities from output in classify.seqs command. [false] |
| Boolean | `relative_abundance` | Value to indicate if summary file values should be displayed as relative abundances rather than raw abundances in classify.seqs command. [false] |
| Int | `classify_seqs_iters` | Number of iterations to perform when calculating the bootstrap confidence score for the taxonomy in classify.seqs command. [100] |
| Int | `tax_level` | Integer parameter representing the level of taxonomic classification to classify up to; a value of -1 will print the max level within the .taxonomy file. [-1] |
| String | `taxon_remove` | Taxa to remove from sequence data in remove.lineage command. ['Chloroplast-Mitochondria-unknown-Archaea-Eukaryota'] |
| Int | `classify_cutoff` | Consensus confidence threshold for taxonomy output. ['51'] |
String | `asv_classify_basis` | The basis parameter indicates what the summary file should represent in classify.otus command. ('sequence', 'otu') ['sequence'] |
| Int | `classify_printlevel` | Taxlevel of the *tax.summary file to print to in classify.otus command. The options are 1 to the max level in the final .taxonomy file; a value of -1 will print the max level within the .taxonomy file. [-1] |
| Boolean | `classify_probabilities` | Parameter that shuts off the outputting of the consensus confidence results in classify.otus command. If true, the confidence will be shown. [true] |
| Boolean | `classify_persample` | Parameter that indicates the assignment of a consensus taxonomy for each group in classify.otus command. [false] |
| Int | `classify_threshold` | Cutoff for the taxonomy input file in classify.otus command. It’s a way to after the fact “adjust” the cutoff used in classify.seqs command without having to reclassify. [90] |
| Float | `cluster_cutoff` | Similarity threshold used for clustering sequences in cluster command. [0.03] |
| String | `cluster_metric` | Metric in the opticluster method in cluster command. ('mcc', 'sens', 'spec', 'tptn', 'fpfn', 'tp', 'tn', 'fp', 'fn', 'f1score', 'accuracy', 'ppv', 'npv', 'fdr') ['mcc'] |
| String | `cluster_initialize` | Initial randomization for the opticluster method in cluster command. ('singleton', 'oneotu') ['singleton'] |
| Float | `cluster_delta` | Stable value for the metric in the opticluster method in cluster command. [0.0001] |
| Int | `cluster_iters` | Maxiters for the opticluster method in cluster command. [100] |
| String | `otu_cluster_method` | Method for clustering sequences based on similarity in cluster command. ('opti', 'average', 'furthest', 'nearest', 'agc' 'dgc', 'unique') ['opti'] |
| String | `otu_classify_basis` | The basis parameter indicates what the summary file should represent in classify.otus command. ('sequence', 'otu') ['otu'] |
| Float | `dist_cutoff` | Cutoff in dist.seqs command; if you know that you are not going to form OTUs with distances larger than 0.10, Mothur can be instructed to not save any distances larger than 0.10. This will significantly cut down on the amount of hard drive space required to store the matrix. [0.03] |
| String | `dist_calc` | Distance calculation method in dist.seqs command. ('onegap', 'nogaps', 'eachgap') ['onegap'] |
| Boolean | `count_ends` | Parameter to penalize gaps that occur at end of sequences in dist.seqs command. [false] |
| Int | `processors` | Number of processors to use in various tasks. [32] |
| String | `container_registry` | Registry that hosts workflow containers ["dnastack"] |

## Workflow outputs

| Type | Output | Description |
| :- | :- | :- |
| File | `make_contigs_fasta` | .fasta file with sequences from make_contigs_scrap_fasta removed during make_contigs task |
| File | `make_contigs_scrap_fasta` | Sequences that were removed from raw reads in make_contigs task |
| File | `align_summary` | Summary of sequence alignment output by align_seqs_to_silva |
| File | `classify_seqs_taxonomy` | Sequence and corresponding taxonomic classification prior to the removal of unwanted lineages |
| File | `remove_lineage_fasta` | .fasta file to be have singletons and doubletons removed and then to be used in subsequent ASV/OTU generation |
| File | `remove_lineage_taxonomy` | Sequence and corresponding taxonomic classification with lineages removed as per the taxon_remove input parameter |
| File | `fasta_singletons_and_doubletons_removed` | Final .fasta file to be used in subsequent ASV/OTU generation |
| File | `count_table_singletons_and_doubletons_removed` | Final .count_table file to be used in subsequent ASV/OTU generation |
| Array[File] | `ASVs` | Array of output files that correspond to ASV data to be used for downstream processing/analysis/visualization |
| Array[File] | `OTUs` | Array of output files that correspond to OTU data to be used for downstream processing/analysis/visualization |

## References

[1] Schloss, P. D., Westcott, S. L., Ryabin, T., Hall, J. R., Hartmann, M., Hollister, E. B., Lesniewski, R. A., Oakley, B. B., Parks, D. H., Robinson, C. J., Sahl, J. W., Stres, B., Thallinger, G. G., Van Horn, D. J., & Weber, C. F. (2009). Introducing mothur: Open-Source, Platform-Independent, Community-Supported Software for Describing and Comparing Microbial Communities. Applied and Environmental Microbiology, 75(23), 7537–7541. https://doi.org/10.1128/AEM.01541-09

[2] Wallen, Z.D., Appah, M., Dean, M.N. et al. Characterizing dysbiosis of gut microbiome in PD: evidence for overabundance of opportunistic pathogens. npj Parkinsons Dis. 6, 11 (2020). https://doi.org/10.1038/s41531-020-0112-6
