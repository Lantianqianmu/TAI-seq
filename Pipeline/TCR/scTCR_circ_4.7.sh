#!/bin/bash

THREAD=24
workingDIR="/data/xrz/TAI-demo"
scriptDIR="/data/xrz/Pipeline"
mixcrDIR="/home/xrz/mixcr-4.7.0"
SPLIT_PATTERN="6 6"
SPECIES=hsa # or mmu


for i in test_circ;
do  
    echo "---------------------------------------------------------${i}------------------------------------------------------------------"
    cd ${workingDIR}/${i}

    rm *mixcr*
    if [[ ! -e ${i}_barcode_TCR_R1.fq.gz ]]; then
        # for cRT45 and Read45
        python ${scriptDIR}/process_fastq_phase.py \
            -1 ${i}_R1.fq.gz \
            -2 ${i}_R2.fq.gz \
            -q ${i}_barcode_TCR_R1.fq.gz \
            -e ${i}_barcode_TCR_R2.fq.gz \
            -c ${THREAD} \
            -d TCR_v2 -m 0 -s 0
    fi

    if [[ ! -e ${i}_barcode_cutadapt_TCR_R1.fq.gz || ! -e ${i}_barcode_cutadapt_TCR_R2.fq.gz ]]; then
        # Read3
        cutadapt \
            -j ${THREAD} -m 2:2 \
            -a "AGATGTGTATAAGAGACAG;max_error_rate=0.2...ACTCTGCGTTGATACCACTGCTT;min_overlap=1" \
            --pair-filter=any \
            -o ${i}_barcode_cutadapt_TCR_R1.fq.gz \
            -p ${i}_barcode_cutadapt_TCR_R2.fq.gz \
            ${i}_barcode_TCR_R1.fq.gz \
            ${i}_barcode_TCR_R2.fq.gz
    fi

    echo `${mixcrDIR}/mixcr -Xmx32g --version`

    # with umi
    ${mixcrDIR}/mixcr align -Xmx256g \
        --preset generic-ht-single-cell-fragmented-with-umi \
        --tag-pattern '^(R1:*)\^(CELL1:N{8})(CELL2:N{8})(CELL3:N{5})(UMI:N{8})' \
        --threads ${THREAD} \
        --rna \
        --species ${SPECIES} \
        --write-all \
        --report ${i}_mixcr_align.report \
        --keep-non-CDR3-alignments \
        --reset-export-clone-grouping \
        --assemble-contigs-by-cells \
        --trimming-quality-threshold 2 \
        --force-overwrite \
        ${i}_barcode_cutadapt_TCR_R1.fq.gz \
        ${i}_barcode_cutadapt_TCR_R2.fq.gz \
        ${i}_mixcr_alignment_umi.vdjca

    ${mixcrDIR}/mixcr refineTagsAndSort -Xmx400g \
        --use-local-temp \
        --report ${i}_mixcr_barcode_correction_report.txt \
        --json-report ${i}_mixcr_barcode_correction_report.json \
        ${i}_mixcr_alignment_umi.vdjca \
        ${i}_mixcr_alignment_corrected_umi.vdjca

    ${mixcrDIR}/mixcr assemblePartial -Xmx256g \
        --report ${i}_mixcr_assemblePartial_report.txt \
        --json-report ${i}_mixcr.assemblePartial_report.json \
        --cell-level \
        ${i}_mixcr_alignment_corrected_umi.vdjca \
        ${i}_mixcr_passembled_umi.vdjca

    ${mixcrDIR}/mixcr extend -Xmx256g \
        --report ${i}_mixcr_extend_report.txt \
        --json-report ${i}_mixcr_extend_report.json \
        ${i}_mixcr_passembled_umi.vdjca \
        ${i}_mixcr_extended_umi.vdjca

    ${mixcrDIR}/mixcr assemble -Xmx400g \
        --write-alignments \
        --dont-infer-threshold \
        --cell-level \
        --report ${i}_mixcr_assemble_report.txt \
        --json-report ${i}_mixcr_assemble_report.json \
        ${i}_mixcr_extended_umi.vdjca \
        ${i}_mixcr_assembled_umi.clna    

    ${mixcrDIR}/mixcr exportClones -Xmx256g -topChains \
        --impute-germline-on-export \
        -cellId dash \
        ${i}_mixcr_assembled_umi.clna \
        ${i}_mixcr_clones_umi.tsv

    ${mixcrDIR}/mixcr assembleContigs -Xmx256g \
        --threads ${THREAD} \
        --report ${i}_mixcr_assembleContigs_report.txt \
        --json-report ${i}_mixcr_assembleContigs_report.json \
        ${i}_mixcr_assembled_umi.clna \
        ${i}_mixcr_contigs_umi.clns

    ${mixcrDIR}/mixcr qc -Xmx64g ${i}_mixcr_contigs_umi.clns

    ${mixcrDIR}/mixcr exportClones -Xmx64g -topChains \
        -f \
        --impute-germline-on-export \
        -cellId dash \
        ${i}_mixcr_contigs_umi.clns \
        ${i}_mixcr_clones_umi.tsv

    ${mixcrDIR}/mixcr exportClonesPretty -Xmx64g \
        -f \
        ${i}_mixcr_contigs_umi.clns \
        ${i}_mixcr_clones_pretty_umi.txt

    # export QC
    ${mixcrDIR}/mixcr exportQc align -Xmx32g \
        --force-overwrite \
        ${i}_mixcr_alignment_corrected_umi.vdjca \
        ${i}_mixcr_alignQc_umi.pdf
    ${mixcrDIR}/mixcr exportQc coverage -Xmx32g \
        --force-overwrite \
        ${i}_mixcr_alignment_corrected_umi.vdjca \
        ${i}_mixcr_coverage_umi.pdf
    ${mixcrDIR}/mixcr exportQc chainUsage -Xmx32g \
        --force-overwrite \
        --hide-non-functional \
        ${i}_mixcr_assembled_umi.clna \
        ${i}_mixcr_chainUsage_umi.pdf 
    ${mixcrDIR}/mixcr exportQc tags -Xmx32g \
        --force-overwrite \
        ${i}_mixcr_assembled_umi.clna \
        ${i}_mixcr_barcodesFiltering_umi.pdf 



    # without umi
    ${mixcrDIR}/mixcr align -Xmx256g \
        --preset generic-ht-single-cell-fragmented \
        --tag-pattern '^(R1:*)\^(CELL1:N{8})(CELL2:N{8})(CELL3:N{5})' \
        --threads ${THREAD} \
        --rna \
        --species ${SPECIES} \
        --write-all \
        --report ${i}_mixcr_align.report \
        --keep-non-CDR3-alignments \
        --reset-export-clone-grouping \
        --assemble-contigs-by-cells \
        --trimming-quality-threshold 2 \
        --force-overwrite \
        ${i}_barcode_cutadapt_TCR_R1.fq.gz \
        ${i}_barcode_cutadapt_TCR_R2.fq.gz \
        ${i}_mixcr_alignment_noumi.vdjca

    ${mixcrDIR}/mixcr refineTagsAndSort -Xmx400g \
        --use-local-temp \
        --report ${i}_mixcr_barcode_correction_report.txt \
        --json-report ${i}_mixcr_barcode_correction_report.json \
        ${i}_mixcr_alignment_noumi.vdjca \
        ${i}_mixcr_alignment_corrected_noumi.vdjca

    ${mixcrDIR}/mixcr assemblePartial -Xmx256g \
        --report ${i}_mixcr_assemblePartial_report.txt \
        --json-report ${i}_mixcr.assemblePartial_report.json \
        --cell-level \
        ${i}_mixcr_alignment_corrected_noumi.vdjca \
        ${i}_mixcr_passembled_noumi.vdjca

    ${mixcrDIR}/mixcr extend -Xmx256g \
        --report ${i}_mixcr_extend_report.txt \
        --json-report ${i}_mixcr_extend_report.json \
        ${i}_mixcr_passembled_noumi.vdjca \
        ${i}_mixcr_extended_noumi.vdjca

    ${mixcrDIR}/mixcr assemble -Xmx400g \
        --dont-infer-threshold \
        --cell-level \
        --write-alignments \
        --report ${i}_mixcr_assemble_report.txt \
        --json-report ${i}_mixcr_assemble_report.json \
        ${i}_mixcr_extended_noumi.vdjca \
        ${i}_mixcr_assembled_noumi.clna    

    ${mixcrDIR}/mixcr exportClones -Xmx256g -topChains \
        --impute-germline-on-export \
        -cellId dash \
        ${i}_mixcr_assembled_noumi.clna \
        ${i}_mixcr_clones_noumi.tsv

    ${mixcrDIR}/mixcr assembleContigs -Xmx256g \
        --threads ${THREAD} \
        --report ${i}_mixcr_assembleContigs_report.txt \
        --json-report ${i}_mixcr_assembleContigs_report.json \
        ${i}_mixcr_assembled_noumi.clna \
        ${i}_mixcr_contigs_noumi.clns

    ${mixcrDIR}/mixcr qc -Xmx64g ${i}_mixcr_contigs_noumi.clns

    ${mixcrDIR}/mixcr exportClones -Xmx64g -topChains \
        -f \
        --impute-germline-on-export \
        -cellId dash \
        ${i}_mixcr_contigs_noumi.clns \
        ${i}_mixcr_clones_noumi.tsv

    ${mixcrDIR}/mixcr exportClonesPretty -Xmx64g \
        -f \
        ${i}_mixcr_contigs_noumi.clns \
        ${i}_mixcr_clones_pretty_noumi.txt

    # export QC
    ${mixcrDIR}/mixcr exportQc align -Xmx32g \
        --force-overwrite \
        ${i}_mixcr_alignment_corrected_noumi.vdjca \
        ${i}_mixcr_alignQc_noumi.pdf
    ${mixcrDIR}/mixcr exportQc coverage -Xmx32g \
        --force-overwrite \
        ${i}_mixcr_alignment_corrected_noumi.vdjca \
        ${i}_mixcr_coverage_noumi.pdf
    ${mixcrDIR}/mixcr exportQc chainUsage -Xmx32g \
        --force-overwrite \
        --hide-non-functional \
        ${i}_mixcr_assembled_noumi.clna \
        ${i}_mixcr_chainUsage_noumi.pdf 
    ${mixcrDIR}/mixcr exportQc tags -Xmx32g \
        --force-overwrite \
        ${i}_mixcr_assembled_noumi.clna \
        ${i}_mixcr_barcodesFiltering_noumi.pdf 



    ###################### Split Round 1 barcode tsv ######################
    echo "####################################################################################################################################"
    echo "spliting round 1 barcode..."
    echo "####################################################################################################################################"
 
    python ${scriptDIR}/split_round1.py \
        -i ${i}_mixcr_clones_umi.tsv \
        -u ${SPLIT_PATTERN} \
        -p ${i}_mixcr_clones_umi \
        -t tsv -l TCR  

    python ${scriptDIR}/split_round1.py \
        -i ${i}_mixcr_clones_noumi.tsv \
        -u ${SPLIT_PATTERN} \
        -p ${i}_mixcr_clones_noumi \
        -t tsv -l TCR   


    echo "[`date`] ${i} finished.\n"
done



