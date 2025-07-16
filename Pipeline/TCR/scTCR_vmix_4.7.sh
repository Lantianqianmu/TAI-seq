#!/bin/bash
THREAD=24
workingDIR="/data/xrz/TAI-demo"
scriptDIR="/data/xrz/Pipeline"
mixcrDIR="/home/xrz/mixcr-4.7.0"
SPLIT_PATTERN="6 6"
SPECIES=hsa # or mmu

for i in test_vmix;
do  
    echo "---------------------------------------------------------${i}------------------------------------------------------------------"
    cd ${workingDIR}/${i}


    # rm *mixcr*
    echo "[`date`] Begin to process ${i}"
    # Allow 0 shift, 0 mismatch, v2 = Read45
    if [[ ! -e ${i}_barcode_TCR_R1.fq.gz ]]; then
        python ${scriptDIR}/process_fastq_phase.py \
            -1 ${i}_R1.fq.gz \
            -2 ${i}_R2.fq.gz \
            -q ${i}_barcode_TCR_R1.fq.gz \
            -e ${i}_barcode_TCR_R2.fq.gz \
            -c ${THREAD} \
            -d Vmix_v2 -s 0 -m 0 -v 2
    fi

    echo `${mixcrDIR}/mixcr -Xmx32g --version`

    # umi
    ${mixcrDIR}/mixcr align -Xmx256g \
        --preset generic-ht-single-cell-amplicon-with-umi \
        --tag-pattern '^(R1:*)\^(CELL1:N{8})(CELL2:N{8})(CELL3:N{5})(UMI:N{8})' \
        --threads ${THREAD} \
        --write-all \
        --trimming-quality-threshold 2 \
        --assemble-clonotypes-by CDR3 \
        --rna \
        --species ${SPECIES} \
        --keep-non-CDR3-alignments \
        --reset-export-clone-grouping \
        --floating-left-alignment-boundary \
        --floating-right-alignment-boundary C \
        --report ${i}_mixcr_align.report \
        ${i}_barcode_TCR_R1.fq.gz \
        ${i}_barcode_TCR_R2.fq.gz \
        ${i}_mixcr_alignment_umi.vdjca

    ${mixcrDIR}/mixcr refineTagsAndSort -Xmx400g \
        -f \
        --use-local-temp \
        --report ${i}_mixcr_barcode_correction_report.txt \
        --json-report ${i}_mixcr_barcode_correction_report.json \
        ${i}_mixcr_alignment_umi.vdjca \
        ${i}_mixcr_alignment_corrected_umi.vdjca

    ${mixcrDIR}/mixcr assemble -Xmx400g \
        -f \
        --write-alignments \
        --dont-infer-threshold \
        --report ${i}_mixcr_assemble_report.txt \
        --json-report ${i}_mixcr_assemble_report.json \
        ${i}_mixcr_alignment_corrected_umi.vdjca \
        ${i}_mixcr_assembled_umi.clna    

    ${mixcrDIR}/mixcr qc -Xmx64g ${i}_mixcr_assembled_umi.clna

    ${mixcrDIR}/mixcr exportClones -Xmx64g -topChains \
        -f \
        --impute-germline-on-export \
        -cellId dash \
        ${i}_mixcr_assembled_umi.clna \
        ${i}_mixcr_clones_umi.tsv


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


    # noumi
    ${mixcrDIR}/mixcr align -Xmx256g \
        --preset generic-ht-single-cell-amplicon \
        --tag-pattern '^(R1:*)\^(CELL1:N{8})(CELL2:N{8})(CELL3:N{5})' \
        --trimming-quality-threshold 2 \
        --assemble-clonotypes-by CDR3 \
        --threads ${THREAD} \
        --write-all \
        --rna \
        --species ${SPECIES} \
        --keep-non-CDR3-alignments \
        --reset-export-clone-grouping \
        --floating-left-alignment-boundary \
        --floating-right-alignment-boundary C \
        --report ${i}_mixcr_align.report \
        ${i}_barcode_TCR_R1.fq.gz \
        ${i}_barcode_TCR_R2.fq.gz \
        ${i}_mixcr_alignment_noumi.vdjca


    ${mixcrDIR}/mixcr refineTagsAndSort -Xmx400g \
        -f \
        --use-local-temp \
        --report ${i}_mixcr_barcode_correction_report.txt \
        --json-report ${i}_mixcr_barcode_correction_report.json \
        ${i}_mixcr_alignment_noumi.vdjca \
        ${i}_mixcr_alignment_corrected_noumi.vdjca

    ${mixcrDIR}/mixcr assemble -Xmx400g \
        -f \
        --write-alignments \
        --dont-infer-threshold \
        --report ${i}_mixcr_assemble_report.txt \
        --json-report ${i}_mixcr_assemble_report.json \
        ${i}_mixcr_alignment_corrected_noumi.vdjca \
        ${i}_mixcr_assembled_noumi.clna    

    ${mixcrDIR}/mixcr qc -Xmx256g ${i}_mixcr_assembled_noumi.clna

    ${mixcrDIR}/mixcr exportClones -Xmx256g -topChains \
        -f \
        --impute-germline-on-export \
        -cellId dash \
        ${i}_mixcr_assembled_noumi.clna \
        ${i}_mixcr_clones_noumi.tsv

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
   
    echo "[`date`] ${i} finished."
    echo ""

done


