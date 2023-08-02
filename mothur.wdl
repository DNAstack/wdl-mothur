version 1.0

workflow mothur {
    input {
        File test_fastq_tar
        File silva_ref_fasta
        File silva_ref_tax
        File oligos
        String file_type
        String prefix

        # make_contigs
        Int max_ambig = 0
        Int max_homop = 8
        Int max_length = 300
        Int primer_differences = 4
        Boolean check_orient = true
        String sequencing_format = "illumina1.8+"
        Int quality_score_threshold = 20

        # pcr_seqs
        Int ref_fasta_start
        Int ref_fasta_end
        Boolean keepdots = false

        # align_seqs_to_silva
        String align_method = "needleman"
        Int match = 1
        Int mismatch = -1
        Int gap_open = -2
        Int gap_extend = -1
        Int ksize = 8
        String search_method = "kmer"
        Boolean flip = true
        Float threshold = 0.5

        # filter_screened_aligned_seqs
        String trump_character = "."
        Boolean vertical = true

        # pre_cluster
        Int sequence_differences = 2
        String precluster_method = "simple"

        # chimera_vsearch
        Boolean remove_chimeras = true
        Boolean dereplicate_chimeras = true
        Float minh_chimeras = 0.3
        Float mindiv_chimeras = 0.5
        Int mindiffs_chimeras = 3
        Float xn_chimeras = 8.0
        Float dn_chimeras = 1.4

        # classify_seqs
        String classify_method = "wang"
        String classify_output_format = "simple"
        Int classify_seqs_cutoff = 80
        Boolean bootstrap_probabilities = false
        Boolean relative_abundance = false
        Int classify_seqs_iters = 100
        Int tax_level = -1

        # remove_lineage
        String taxon_remove = "Chloroplast-Mitochondria-unknown-Archaea-Eukaryota"

        # make_ASVs
        Int classify_cutoff = 51
        String asv_classify_basis = "sequence"
        Int classify_printlevel = -1
        Boolean classify_probabilities = true
        Boolean classify_persample = false
        Int classify_threshold = 90

        # make_OTUs
        Float cluster_cutoff = 0.03
        String cluster_metric = "mcc"
        String cluster_initialize = "singleton"
        Float cluster_delta = 0.0001
        Int cluster_iters = 100
        String otu_cluster_method = "opti"
        String otu_classify_basis = "otu"
        Float dist_cutoff = 0.03
        String dist_calc = "onegap"
        Boolean count_ends = false

        # Runtime configuration
        Int processors = 32
        String container_registry = "dnastack"
    }

    call make_contigs {
        input:
            prefix = prefix,
            test_fastq_tar = test_fastq_tar,
            file_type = file_type,
            oligos = oligos,
            max_ambig = max_ambig,
            max_homop = max_homop,
            max_length = max_length,
            primer_differences = primer_differences,
            check_orient = check_orient,
            sequencing_format = sequencing_format,
            quality_score_threshold = quality_score_threshold,
            align_method = align_method,
            match = match,
            mismatch = mismatch,
            gap_open = gap_open,
            gap_extend = gap_extend,
            processors = processors,
            container_registry = container_registry
    }

    call unique_seqs {
        input:
            prefix = prefix,
            fasta = make_contigs.make_contigs_fasta,
            count_table = make_contigs.make_contigs_count_table,
            container_registry = container_registry
    }

    call pcr_seqs {
        input:
            silva_ref_fasta = silva_ref_fasta,
            ref_fasta_start = ref_fasta_start,
            ref_fasta_end = ref_fasta_end,
            silva_ref_tax = silva_ref_tax,
            keepdots = keepdots,
            processors = processors,
            container_registry = container_registry
    }

    call align_seqs_to_silva {
        input:
            prefix = prefix,
            unique_seqs_fasta = unique_seqs.unique_seqs_fasta,
            silva_v4_fasta = pcr_seqs.silva_v4_fasta,
            align_method = align_method,
            match = match,
            mismatch = mismatch,
            gap_open = gap_open,
            gap_extend = gap_extend,
            ksize = ksize,
            search_method = search_method,
            flip = flip,
            threshold = threshold,
            processors = processors,
            container_registry = container_registry
    }

    call screen_aligned_seqs {
        input:
            prefix = prefix,
            aligned_unique_seqs_fasta = align_seqs_to_silva.aligned_unique_seqs_fasta,
            container_registry = container_registry
    }

    call filter_screened_aligned_seqs {
        input:
            prefix = prefix,
            aligned_screened_unique_seqs_fasta = screen_aligned_seqs.aligned_screened_unique_seqs_fasta,
            count_table = unique_seqs.unique_seqs_count_table,
            trump_character = trump_character,
            vertical = vertical,
            processors = processors,
            container_registry = container_registry
    }

    call pre_cluster {
        input:
            prefix = prefix,
            unique_filtered_screened_aligned_seqs_fasta = filter_screened_aligned_seqs.unique_filtered_screened_aligned_seqs_fasta,
            unique_filtered_screened_aligned_seqs_count_table = filter_screened_aligned_seqs.unique_filtered_screened_aligned_seqs_count_table,
            align_method = align_method,
            match = match,
            mismatch = mismatch,
            gap_open = gap_open,
            gap_extend = gap_extend,
            sequence_differences = sequence_differences,
            precluster_method = precluster_method,
            processors = processors,
            container_registry = container_registry
    }

    call chimera_vsearch {
        input:
            pre_cluster_fasta = pre_cluster.pre_cluster_fasta,
            count_table = pre_cluster.pre_cluster_count_table,
            remove_chimeras = remove_chimeras,
            dereplicate_chimeras = dereplicate_chimeras,
            minh_chimeras = minh_chimeras,
            mindiv_chimeras = mindiv_chimeras,
            mindiffs_chimeras = mindiffs_chimeras,
            xn_chimeras = xn_chimeras,
            dn_chimeras = dn_chimeras,
            container_registry = container_registry
    }

    call classify_seqs {
        input:
            prefix = prefix,
            chimera_vsearch_fasta = chimera_vsearch.chimera_vsearch_fasta,
            count_table = chimera_vsearch.chimera_vsearch_count_table,
            silva_v4_fasta = pcr_seqs.silva_v4_fasta,
            silva_ref_tax = silva_ref_tax,
            ksize = ksize,
            search_method = search_method,
            classify_method = classify_method,
            classify_output_format = classify_output_format,
            classify_seqs_cutoff = classify_seqs_cutoff,
            bootstrap_probabilities = bootstrap_probabilities,
            relative_abundance = relative_abundance,
            classify_seqs_iters = classify_seqs_iters,
            tax_level = tax_level,
            container_registry = container_registry
    }

    call remove_lineage {
        input:
            chimera_vsearch_fasta = chimera_vsearch.chimera_vsearch_fasta,
            count_table = chimera_vsearch.chimera_vsearch_count_table,
            taxonomy = classify_seqs.taxonomy,
            taxon_remove = taxon_remove,
            container_registry = container_registry
    }

    call remove_singletons_and_doubletons {
        input:
            prefix = prefix,
            fasta = remove_lineage.remove_lineage_fasta,
            count_table = remove_lineage.remove_lineage_count_table,
            container_registry = container_registry
    }

    call make_ASVs {
        input:
            count_table = remove_singletons_and_doubletons.count_table_singletons_and_doubletons_removed,
            taxonomy = remove_lineage.remove_lineage_taxonomy,
            classify_output_format = classify_output_format,
            relative_abundance = relative_abundance,
            classify_cutoff = classify_cutoff,
            asv_classify_basis = asv_classify_basis,
            classify_printlevel = classify_printlevel,
            classify_probabilities = classify_probabilities,
            classify_persample = classify_persample,
            classify_threshold = classify_threshold,
            container_registry = container_registry
    }

    call make_OTUs {
        input:
            fasta = remove_singletons_and_doubletons.fasta_singletons_and_doubletons_removed,
            count_table = remove_singletons_and_doubletons.count_table_singletons_and_doubletons_removed,
            taxonomy = remove_lineage.remove_lineage_taxonomy,
            classify_output_format = classify_output_format,
            relative_abundance = relative_abundance,
            classify_cutoff = classify_cutoff,
            classify_printlevel = classify_printlevel,
            classify_probabilities = classify_probabilities,
            classify_persample = classify_persample,
            classify_threshold = classify_threshold,
            cluster_cutoff = cluster_cutoff,
            cluster_metric = cluster_metric,
            cluster_initialize = cluster_initialize,
            cluster_delta = cluster_delta,
            cluster_iters = cluster_iters,
            otu_cluster_method = otu_cluster_method,
            otu_classify_basis = otu_classify_basis,
            dist_cutoff = dist_cutoff,
            dist_calc = dist_calc,
            count_ends = count_ends,
            processors = processors,
            container_registry = container_registry
    }

    output {
        # make_contigs
        File make_contigs_fasta = make_contigs.make_contigs_fasta
        File make_contigs_scrap_fasta = make_contigs.make_contigs_scrap_fasta

        # align_seqs_to_silva
        File align_summary = align_seqs_to_silva.align_summary

        # classify_seqs
        File classify_seqs_taxonomy = classify_seqs.taxonomy

        # remove_lineage
        File remove_lineage_fasta = remove_lineage.remove_lineage_fasta
        File remove_lineage_taxonomy = remove_lineage.remove_lineage_taxonomy

        # remove_singletons_and_doubletons
        File fasta_singletons_and_doubletons_removed = remove_singletons_and_doubletons.fasta_singletons_and_doubletons_removed
        File count_table_singletons_and_doubletons_removed = remove_singletons_and_doubletons.count_table_singletons_and_doubletons_removed

        # make_ASVs
        Array[File] ASVs = make_ASVs.ASVs

        # make_OTUs
        Array[File] OTUs = make_OTUs.OTUs
    }

    meta {
        author: "Jesse Wolf"
        email: "jesse@dnastack.com"
        description: "Use raw 16S rRNA reads in Mothur to characterize gut microbial communities"
    }

    parameter_meta {
        # Inputs
        test_fastq_tar: {help: "Tarball of all paired end fastq files to analyzed"}
        silva_ref_fasta: {help: "Reference sequence database from SILVA"}
        silva_ref_tax: {help: "Reference taxonomy database from SILVA "}
        oligos: {help: "File with primers used in sequencing"}
        file_type: {help: "File type is fastq"}
        prefix: {help: "Prefix for output files"}
        ref_fasta_start: {help: "Start position to trim reference fasta"}
        ref_fasta_end: {help: "End position to trim reference fasta"}

        # Outputs
        make_contigs_fasta: {help:".fasta file with sequences from make_contigs_scrap_fasta removed during make_contigs task"}
        make_contigs_scrap_fasta: {help:"Sequences that were removed from raw reads in make_contigs task"}
        align_summary: {help: "Summary of sequence alignment output by align_seqs_to_silva"}
        classify_seqs_taxonomy : {help: "Sequence and corresponding taxonomic classification prior to the removal of unwanted lineages output by classify_seqs"}
        remove_lineage_fasta: {help: "Final .fasta file to be used in subsequent ASV/OTU generation output by remove_lineage"}
        remove_lineage_taxonomy: {help: "Sequence and corresponding taxonomic classification with lineages removed output by remove_lineage"}
        fasta_singletons_and_doubletons_removed: {help: "Final .fasta file to be used in subsequent ASV/OTU generation"}
        count_table_singletons_and_doubletons_removed:{help: "Final .count_table file to be used in subsequent ASV/OTU generation"}
        ASVs: {help: "Array of output files that correspond to ASV data to be used for downstream processing/analysis/visualization"}
        OTUs: {help: "Array of output files that correspond to OTU data to be used for downstream processing/analysis/visualization"}
    }
}

task make_contigs {
    input {
        String prefix
        File test_fastq_tar
        String file_type

        File oligos

        Int max_ambig
        Int max_homop
        Int max_length
        Int primer_differences
        Boolean check_orient
        String sequencing_format
        Int quality_score_threshold

        String align_method
        Int match
        Int mismatch
        Int gap_open
        Int gap_extend

        Int processors
        String container_registry
    }

    Int disk_size_make_contigs = ceil(size(test_fastq_tar, "GB") * 20 + 50)

    command <<<
        set -euo pipefail

        mkdir fastq_files
        tar -xf ~{test_fastq_tar} -C fastq_files

        for file in fastq_files/*; do
            gunzip "${file}";
        done

        mothur "#set.dir(output=output_makeContigs);
        make.file(inputdir=fastq_files, type=~{file_type}, prefix=~{prefix});
        make.contigs(file=~{prefix}.files, processors=~{processors}, oligos=~{oligos}, maxambig=~{max_ambig}, pdiffs=~{primer_differences}, maxlength=~{max_length}, maxhomop=~{max_homop}, align=~{align_method}, checkorient=~{check_orient}, format=~{sequencing_format}, match=~{match}, mismatch=~{mismatch}, gapopen=~{gap_open}, gapextend=~{gap_extend}, insert=~{quality_score_threshold});
        summary.seqs()"
    >>>

    output {
        File stability = "output_makeContigs/~{prefix}.files"
        File make_contigs_fasta = "output_makeContigs/~{prefix}.trim.contigs.fasta"
        File make_contigs_count_table = "output_makeContigs/~{prefix}.contigs.count_table"
        File make_contigs_report = "output_makeContigs/~{prefix}.contigs_report"
        File make_contigs_scrap_fasta = "output_makeContigs/~{prefix}.scrap.contigs.fasta"
    }

    runtime {
        docker: "~{container_registry}/mothur:1.48.0"
        cpu: processors
        memory: "64 GB"
        disks: "local-disk " + disk_size_make_contigs + " HDD"
        preemptible: 3
    }
}

task unique_seqs {
    input {
        String prefix
        File fasta
        File count_table

        String container_registry
    }

    Int disk_size_unique_seqs = ceil(size(fasta, "GB") * 3 + 30)

    command <<<
        set -euo pipefail

        mothur "#set.dir(output=output_uniqueSeqs);
        unique.seqs(fasta=~{fasta}, count=~{count_table});
        summary.seqs()"
    >>>

    output {
        File unique_seqs_fasta = "output_uniqueSeqs/~{prefix}.trim.contigs.unique.fasta"
        File unique_seqs_count_table = "output_uniqueSeqs/~{prefix}.trim.contigs.count_table"
    }

    runtime {
        docker: "~{container_registry}/mothur:1.48.0"
        cpu: 4
        memory: "64 GB"
        disks: "local-disk " + disk_size_unique_seqs + " HDD"
        preemptible: 3
    }
}

task pcr_seqs {
    input {
        File silva_ref_fasta
        Int ref_fasta_start
        Int ref_fasta_end
        File silva_ref_tax

        Boolean keepdots

        Int processors
        String container_registry
    }

    String silva_ref_fasta_basename = basename(silva_ref_fasta, ".align")
    Int disk_size_pcr_seqs = ceil(size(silva_ref_fasta, "GB") * 3 + 30)

    command <<<
        set -euo pipefail

        mothur "#set.dir(output=output_pcr_seqs);
        pcr.seqs(fasta=~{silva_ref_fasta}, start=~{ref_fasta_start}, end=~{ref_fasta_end}, keepdots=~{keepdots}, processors=~{processors}, taxonomy=~{silva_ref_tax});
        summary.seqs()"
    >>>

    output {
        File silva_v4_fasta = "output_pcr_seqs/~{silva_ref_fasta_basename}.pcr.align"
    }

    runtime {
        docker: "~{container_registry}/mothur:1.48.0"
        cpu: processors
        memory: "2 GB"
        disks: "local-disk " + disk_size_pcr_seqs + " HDD"
        preemptible: 3
    }
}

task align_seqs_to_silva {
    input {
        String prefix

        File unique_seqs_fasta
        File silva_v4_fasta

        String align_method
        Int match
        Int mismatch
        Int gap_open
        Int gap_extend
        Int ksize
        String search_method
        Boolean flip
        Float threshold

        Int processors
        String container_registry
    }

    Int disk_size_align_seqs = ceil(size(silva_v4_fasta, "GB") * 3 + 30)

    command <<<
        set -euo pipefail

        mothur "#set.dir(output=output_align_seqs_to_silva);
        align.seqs(fasta=~{unique_seqs_fasta}, reference=~{silva_v4_fasta}, search=~{search_method}, align=~{align_method}, ksize=~{ksize}, match=~{match}, mismatch=~{mismatch}, gapopen=~{gap_open}, gapextend=~{gap_extend}, flip=~{flip}, threshold=~{threshold}, processors=~{processors});
        summary.seqs()"
    >>>

    output {
        File aligned_unique_seqs_fasta = "output_align_seqs_to_silva/~{prefix}.trim.contigs.unique.align"
        File align_summary = "output_align_seqs_to_silva/~{prefix}.trim.contigs.unique.summary"
    }

    runtime {
        docker: "~{container_registry}/mothur:1.48.0"
        cpu: processors
        memory: "16 GB"
        disks: "local-disk " + disk_size_align_seqs + " HDD"
        preemptible: 3
    }
}

task screen_aligned_seqs {
    input {
        String prefix
        File aligned_unique_seqs_fasta

        String container_registry
    }

    Int disk_size_screen_aligned_seqs = ceil(size(aligned_unique_seqs_fasta, "GB") * 3 + 30)

    command <<<
        set -euo pipefail

        mothur "#set.dir(output=output_screen_aligned_seqs);
        screen.seqs(fasta=~{aligned_unique_seqs_fasta}, optimize=start-end-minlength-maxlength, criteria=90);
        summary.seqs()"
    >>>

    output {
        File aligned_screened_unique_seqs_fasta = "output_screen_aligned_seqs/~{prefix}.trim.contigs.unique.good.align"
    }

    runtime {
        docker: "~{container_registry}/mothur:1.48.0"
        cpu: 32
        memory: "4 GB"
        disks: "local-disk " + disk_size_screen_aligned_seqs + " HDD"
        preemptible: 3
    }
}

task filter_screened_aligned_seqs {
    input {
        String prefix
        File aligned_screened_unique_seqs_fasta
        File count_table

        String trump_character
        Boolean vertical

        Int processors
        String container_registry
    }

    Int disk_size_filter_screened_aligned_seqs = ceil(size(aligned_screened_unique_seqs_fasta, "GB") * 3 + 30)

    command <<<
        set -euo pipefail

        mothur "#set.dir(output=output_filter_screened_aligned_seqs);
        filter.seqs(fasta=~{aligned_screened_unique_seqs_fasta}, vertical=~{vertical}, trump=~{trump_character}, processors=~{processors});
        unique.seqs(count=~{count_table});
        summary.seqs()"
    >>>

    output {
        File unique_filtered_screened_aligned_seqs_fasta = "output_filter_screened_aligned_seqs/~{prefix}.trim.contigs.unique.good.filter.unique.fasta"
        File unique_filtered_screened_aligned_seqs_count_table = "output_filter_screened_aligned_seqs/~{prefix}.trim.contigs.unique.good.filter.count_table"
    }

    runtime {
        docker: "~{container_registry}/mothur:1.48.0"
        cpu: processors
        memory: "16 GB"
        disks: "local-disk " + disk_size_filter_screened_aligned_seqs + " HDD"
        preemptible: 3
    }
}

task pre_cluster {
    input {
        String prefix
        File unique_filtered_screened_aligned_seqs_fasta
        File unique_filtered_screened_aligned_seqs_count_table

        String align_method
        Int match
        Int mismatch
        Int gap_open
        Int gap_extend

        Int sequence_differences
        String precluster_method

        Int processors
        String container_registry
    }

    Int disk_size_pre_cluster = ceil(size(unique_filtered_screened_aligned_seqs_fasta, "GB") * 3 + 30)

    command <<<
        set -euo pipefail

        mothur "#set.dir(output=output_pre_cluster);
        pre.cluster(fasta=~{unique_filtered_screened_aligned_seqs_fasta}, count=~{unique_filtered_screened_aligned_seqs_count_table}, diffs=~{sequence_differences}, align=~{align_method}, match=~{match}, mismatch=~{mismatch}, gapopen=~{gap_open}, gapextend=~{gap_extend}, method=~{precluster_method}, processors=~{processors})"
    >>>

    output {
        File pre_cluster_count_table = "output_pre_cluster/~{prefix}.trim.contigs.unique.good.filter.unique.precluster.count_table"
        File pre_cluster_fasta = "output_pre_cluster/~{prefix}.trim.contigs.unique.good.filter.unique.precluster.fasta"
    }

    runtime {
        docker: "~{container_registry}/mothur:1.48.0"
        cpu: processors
        memory: "16 GB"
        disks: "local-disk " + disk_size_pre_cluster + " HDD"
        preemptible: 3
    }
}

task chimera_vsearch {
    input {
        File pre_cluster_fasta
        File count_table

        Boolean remove_chimeras
        Boolean dereplicate_chimeras

        Float minh_chimeras
        Float mindiv_chimeras
        Int mindiffs_chimeras

        Float xn_chimeras
        Float dn_chimeras

        String container_registry
    }

    Int processors = 16
    Int disk_size_chimera_vsearch = ceil(size(pre_cluster_fasta, "GB") * 3 + 30)

    command <<<
        set -euo pipefail

        # Symbolic linking is necessary here; mothur will not accept '-' in the filenames, which is introduced when referencing paths in Cromwell
        ln -s ~{pre_cluster_fasta} precluster.fasta
        ln -s ~{count_table} precluster.count_table

        mothur "#set.dir(output=output_chimera_vsearch);
        chimera.vsearch(fasta=precluster.fasta, count=precluster.count_table, dereplicate=~{dereplicate_chimeras}, removechimeras=~{remove_chimeras}, processors=~{processors}, minh=~{minh_chimeras}, mindiv=~{mindiv_chimeras}, mindiffs=~{mindiffs_chimeras}, xn=~{xn_chimeras}, dn=~{dn_chimeras});
        summary.seqs()"
    >>>

    output {
        File chimera_vsearch_chimeras = "output_chimera_vsearch/precluster.denovo.vsearch.chimeras"
        File chimera_vsearch_fasta = "output_chimera_vsearch/precluster.denovo.vsearch.fasta"
        File chimera_vsearch_count_table = "output_chimera_vsearch/precluster.denovo.vsearch.count_table"
    }

    runtime {
        docker: "~{container_registry}/mothur:1.48.0"
        cpu: processors
        memory: "8 GB"
        disks: "local-disk " + disk_size_chimera_vsearch + " HDD"
        preemptible: 3
    }
}

task classify_seqs {
    input {
        String prefix
        File chimera_vsearch_fasta
        File count_table

        File silva_v4_fasta
        File silva_ref_tax

        Int ksize
        String search_method

        String classify_method
        String classify_output_format
        Int classify_seqs_cutoff
        Boolean bootstrap_probabilities
        Boolean relative_abundance
        Int classify_seqs_iters
        Int tax_level

        String container_registry
    }

    Int processors = 16
    Int disk_classify_seqs = ceil(size(chimera_vsearch_fasta, "GB") + size(silva_v4_fasta, "GB") + size(silva_ref_tax, "GB") * 3 + 30)

    command <<<
        set -euo pipefail

        mothur "#set.dir(output=output_classify_seqs);
        classify.seqs(fasta=~{chimera_vsearch_fasta}, count=~{count_table}, reference=~{silva_v4_fasta}, taxonomy=~{silva_ref_tax}, method=~{classify_method}, cutoff=~{classify_seqs_cutoff}, probs=~{bootstrap_probabilities}, processors=~{processors}, output=~{classify_output_format}, printlevel=~{tax_level}, iters=~{classify_seqs_iters}, search=~{search_method}, ksize=~{ksize}, relabund=~{relative_abundance});
        rename.file(taxonomy=current, prefix=~{prefix})"
    >>>

    output {
        File taxonomy = "output_classify_seqs/~{prefix}.taxonomy"
    }

    runtime {
        docker: "~{container_registry}/mothur:1.48.0"
        cpu: processors
        memory: "8 GB"
        disks: "local-disk " + disk_classify_seqs + " HDD"
        preemptible: 3
    }
}

task remove_lineage {
    input {
        File chimera_vsearch_fasta
        File count_table

        File taxonomy
        String taxon_remove

        String container_registry
    }

    Int disk_size_remove_lineage = ceil(size(chimera_vsearch_fasta, "GB")+ size(taxonomy, "GB") * 3 + 30)

    command <<<
        set -euo pipefail

        # Symbolic linking is necessary here; mothur will not accept '-' in the filenames, which is introduced when referencing paths in Cromwell
        ln -s ~{chimera_vsearch_fasta} chimera_vsearch.fasta
        ln -s ~{count_table} chimera_vsearch.count_table
        ln -s ~{taxonomy} classify_seqs.taxonomy

        mothur "#set.dir(output=output_remove_lineage);
        remove.lineage(fasta=chimera_vsearch.fasta, count=chimera_vsearch.count_table, taxonomy=classify_seqs.taxonomy, taxon=~{taxon_remove});
        summary.tax();
        rename.file(fasta=current, count=current, taxonomy=current, prefix=final)"
    >>>

    output {
        File remove_lineage_taxonomy = "output_remove_lineage/final.taxonomy"
        File remove_lineage_fasta = "output_remove_lineage/final.fasta"
        File remove_lineage_count_table = "output_remove_lineage/final.count_table"
    }

    runtime {
        docker: "~{container_registry}/mothur:1.48.0"
        cpu: 4
        memory: "4 GB"
        disks: "local-disk " + disk_size_remove_lineage + " HDD"
        preemptible: 3
    }
}

task remove_singletons_and_doubletons {
    input {
        String prefix
        File fasta
        File count_table

        String container_registry
    }

    command <<<
        set -euo pipefail

        awk -F'\t' 'NR > 2 && $2 == 1 { print $1 }' ~{count_table} > ~{prefix}_singles.accnos
        awk -F'\t' 'NR > 2 && $2 == 2 { print $1 }' ~{count_table} > ~{prefix}_doubles.accnos
        cat ~{prefix}_doubles.accnos ~{prefix}_singles.accnos > ~{prefix}_combined.accnos


        # Symbolic linking is necessary here; mothur will not accept '-' in the filenames, which is introduced when referencing paths in Cromwell
        ln -s ~{fasta} remove_lineage.fasta
        ln -s ~{count_table} remove_lineage.count_table

        mothur "#set.dir(output=output_remove_singletons_and_doubletons);
        remove.seqs(fasta=remove_lineage.fasta, count=remove_lineage.count_table, accnos=~{prefix}_combined.accnos);
        rename.file(fasta=current, count=current, prefix=final_singletons_doubletons_removed)"
    >>>

    output {
        File fasta_singletons_and_doubletons_removed = "output_remove_singletons_and_doubletons/final_singletons_doubletons_removed.fasta"
        File count_table_singletons_and_doubletons_removed = "output_remove_singletons_and_doubletons/final_singletons_doubletons_removed.count_table"
    }

    runtime {
        docker: "~{container_registry}/mothur:1.48.0"
        cpu: 4
        memory: "8 GB"
        disks: "local-disk 250 HDD"
        preemptible: 3
    }
}

task make_ASVs {
    input {
        File count_table

        File taxonomy

        String classify_output_format
        Boolean relative_abundance

        Int classify_cutoff
        String asv_classify_basis
        Int classify_printlevel
        Boolean classify_probabilities
        Boolean classify_persample
        Int classify_threshold

        String container_registry
    }

    command <<<
        set -euo pipefail

        mothur "#set.dir(output=output_ASVs);
        set.seed(seed=444);
        make.shared(count=~{count_table});
        classify.otu(list=current, count=~{count_table}, taxonomy=~{taxonomy}, label=ASV, cutoff=~{classify_cutoff}, basis=~{asv_classify_basis}, relabund=~{relative_abundance}, output=~{classify_output_format}, printlevel=~{classify_printlevel}, probs=~{classify_probabilities}, persample=~{classify_persample}, threshold=~{classify_threshold})"
    >>>

    output {
        Array[File] ASVs = glob("output_ASVs/*")
    }

    runtime {
        docker: "~{container_registry}/mothur:1.48.0"
        cpu: 16
        memory: "16 GB"
        disks: "local-disk 250 HDD"
        preemptible: 3
    }
}

task make_OTUs {
    input {
        File fasta
        File count_table

        File taxonomy

        String classify_output_format
        Boolean relative_abundance

        Int classify_cutoff
        Int classify_printlevel
        Boolean classify_probabilities
        Boolean classify_persample
        Int classify_threshold

        Float cluster_cutoff
        String cluster_metric
        String cluster_initialize
        Float cluster_delta
        Int cluster_iters
        String otu_cluster_method
        String otu_classify_basis
        Float dist_cutoff
        String dist_calc
        Boolean count_ends

        Int processors
        String container_registry
    }

    command <<<
        set -euo pipefail

        mothur "#set.dir(output=output_OTUs_dist);
        set.seed(seed=444);
        dist.seqs(fasta=~{fasta}, processors=~{processors}, cutoff=~{dist_cutoff}, countends=~{count_ends}, calc=~{dist_calc});
        cluster(column=current, count=~{count_table}, method=~{otu_cluster_method}, cutoff=~{cluster_cutoff},initialize=~{cluster_initialize}, iters=~{cluster_iters}, delta=~{cluster_delta}, metric=~{cluster_metric});
        make.shared(list=current, count=~{count_table}, label=~{dist_cutoff});
        classify.otu(list=current, count=~{count_table}, taxonomy=~{taxonomy}, label=~{dist_cutoff}, cutoff=~{classify_cutoff}, basis=~{otu_classify_basis}, relabund=~{relative_abundance}, output=~{classify_output_format}, printlevel=~{classify_printlevel}, probs=~{classify_probabilities}, persample=~{classify_persample}, threshold=~{classify_threshold})"
    >>>

    output {
        Array[File] OTUs = glob("output_OTUs_dist/*")
    }

    runtime {
        docker: "~{container_registry}/mothur:1.48.0"
        cpu: processors
        memory: "32 GB"
        disks: "local-disk 500 HDD"
        preemptible: 3
    }
}
