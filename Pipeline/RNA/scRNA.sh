#!/bin/bash
# # # # # CRITICAL # # # # # 

genome="mm10"
workingDIR="/data/xrz/TAI-demo"
scriptDIR="/data/xrz/Pipeline"
gtfDir="/data/xrz/ref/${genome}/${genome}_annotated.gtf"
genomeDir="/data/xrz/ref/${genome}/${genome}-STAR/"
SPLIT_NUM=2
SPLIT_PATTERN="6 6"
minUMI=0
THREAD=24
sample="test_rna"



keepOneReadUMI="T"
keepMultiMapping="F"
allowMultiAssignment="F"
useUMITOOLS="F"
calculateAlignmentDistribution="F"
N6="F"

echo "SPLIT_NUM=${SPLIT_NUM}"
echo "minUMI=${minUMI}"
echo "keepMultiMapping=${keepMultiMapping}"
echo "allowMultiAssignment=${allowMultiAssignment}"
echo "useUMITOOLS=${useUMITOOLS}"
echo "THREAD=${THREAD}"
echo "calculateAlignmentDistribution=${calculateAlignmentDistribution}"
echo "N6=${N6}"


for i in ${sample}; do
    echo "---------------------------------------------------------${i}------------------------------------------------------------------"   
    # extract barcodes from fastq files
    cd ${workingDIR}/${i}
    
    echo "Begin prcessing ${i}: `date`"


    if [[ ! -e ${i}_barcode_RNA_R1.fq.gz || ! -e ${i}_barcode_RNA_R2.fq.gz ]]; then
        # Allow 0 shift, 0 mismatch, v2 = Read45
        # RNA_n6_test: Allow N6 only in 1-48 but not 49-96
        python ${scriptDIR}/process_fastq_phase.py \
            -1 ${i}_R1.fq.gz \
            -2 ${i}_R2.fq.gz \
            -q ${i}_barcode_RNA_R1.fq.gz \
            -e ${i}_barcode_RNA_R2.fq.gz \
            -c ${THREAD} \
            -d Vmix_v2 -s 0 -m 0 -v 2
            # -d RNA_nolibbc -s 0 -m 0 -v 2 -p 138
    fi

    # # Remove adaptor sequences. CRITICAL: cutadapt >= 4.1
    # if [[ ! -e ${i}_barcode_cutadapt_RNA_R1.fq.gz ]]; then
    #     # Read1
    #     cutadapt -j ${THREAD} -m 2 -g "GCAGCGTCAGATGTGTATAAGAGACAG...CCATCTGAGCCACCAGCG;optional;min_overlap=1;max_error_rate=0.2" ${i}_barcode_RNA_R1.fq.gz | \
    #         cutadapt -j ${THREAD} -m 2 -g "AAGCAGTGGTATCAACGCAGAGT" - | \
    #         cutadapt -j ${THREAD} -u 7 -m 2 -a "A{15};min_overlap=1" -o ${i}_barcode_cutadapt_RNA_R1.fq.gz -
    # fi

    # # Remove adaptor sequences. CRITICAL: cutadapt >= 4.1
    # if [[ ! -e ${i}_barcode_cutadapt_RNA_R1.fq.gz ]]; then
    #     # Read3
    #     cutadapt -j ${THREAD} -m 2 -g "ACGAGACCGGAAAGATGTGTATAAGAGACAG...CCATCTGAGCCACCAGCG;optional;min_overlap=1;max_error_rate=0.2" ${i}_barcode_RNA_R1.fq.gz | \
    #         cutadapt -j ${THREAD} -m 2 -g "AAGCAGTGGTATCAACGCAGAGT" - | \
    #         cutadapt -j ${THREAD} -u 7 -m 2 -a "A{15};min_overlap=1" -o ${i}_barcode_cutadapt_RNA_R1.fq.gz -
    # fi

    # Remove adaptor sequences. CRITICAL: cutadapt >= 4.1
    if [[ ! -e ${i}_barcode_cutadapt_RNA_R1.fq.gz ]]; then
        # Read3
        cutadapt -j ${THREAD} -m 2 -g "AGATGTGTATAAGAGACAG...CCATCTGAGCCACCAGCG;optional;min_overlap=1;max_error_rate=0.2" ${i}_barcode_RNA_R1.fq.gz | \
            cutadapt -j ${THREAD} -m 2 -a "A{15};min_overlap=1" -o ${i}_barcode_cutadapt_RNA_R1.fq.gz -
    fi

    # Align with star.
    if [[ ! -e ${i}_Aligned.bam ]]; then
        STAR \
            --runMode alignReads \
            --runThreadN ${THREAD} \
            --readFilesCommand zcat \
            --outFilterMultimapNmax 20 \
            --outFilterScoreMinOverLread 0.3 \
            --outFilterMatchNminOverLread 0.3 \
            --outSAMattributes NH HI AS nM MD \
            --limitOutSJcollapsed 5000000 \
            --outReadsUnmapped None \
            --chimOutType WithinBAM \
            --readFilesIn ${i}_barcode_cutadapt_RNA_R1.fq.gz \
            --genomeDir ${genomeDir} \
            --outFileNamePrefix ${i}_ \
            --outSAMtype BAM SortedByCoordinate \
            --outBAMsortingThreadN ${THREAD}
        samtools view -@ ${THREAD} -h ${i}_Aligned.sortedByCoord.out.bam | \
            perl -F'\t' -alne 'BEGIN{$|++;}if(/^@/){print;}else{chomp; if($F[0] =~ /([A|G|C|T]{8}-[A|G|C|T]{8}-[A|G|C|T]{5}-{0,1}[A|G|C|T]{0,3}):([A|G|C|T|N]{8,10})/){print "$_\tBC:Z:$1\tUB:Z:$2";}}' | \
            samtools view -@ ${THREAD} -h -b > ${i}_Aligned.bam
        rm ${i}_Aligned.sortedByCoord.out.bam
        samtools index -@ ${THREAD} ${i}_Aligned.bam
    else
        echo "${i}_Aligned.bam already exists."
    fi
  

    if [ ${keepMultiMapping} == "F" ]; then
        echo "Multi-mapping reads are filtered."
        samtools view -q 30 -b ${i}_Aligned.bam > ${i}_Aligned_filtered.bam
        featureCounts -s 1 -Q 30 -T ${THREAD} -t gene -g gene_name -a ${gtfDir} -o ${i}_featureCounts_report.txt -R BAM ${i}_Aligned_filtered.bam
        samtools sort -@ ${THREAD} ${i}_Aligned_filtered.bam.featureCounts.bam -o ${i}_featureCounts_sorted.bam
        samtools index -@ ${THREAD} ${i}_featureCounts_sorted.bam
        rm ${i}_Aligned_filtered.bam.featureCounts.bam

        if [ ${calculateAlignmentDistribution} == "T" ]; then
            featureCounts -s 1 -Q 30 -T ${THREAD} -t exon -g exon_type -a ${gtfDir} -o ${i}_exon.txt -R BAM ${i}_Aligned_filtered.bam
            samtools sort -n -@ ${THREAD} ${i}_Aligned_filtered.bam.featureCounts.bam -o ${i}_exon.bam
            samtools sort -n -@ ${THREAD} ${i}_featureCounts_sorted.bam -o ${i}_gene.bam
            perl ${scriptDIR}/RNA/alignment_distribution.pl -e ${i}_exon.bam -g ${i}_gene.bam -o ${i}_alignment_distribution.txt
        fi
    else
        echo "Multi-mapping reads are retained."
        featureCounts -s 1 -Q 0 -M -T ${THREAD} -t gene -g gene_name -a ${gtfDir} -o ${i}_featureCounts_report.txt -R BAM ${i}_Aligned.bam
        samtools sort -@ ${THREAD} ${i}_Aligned.bam.featureCounts.bam -o ${i}_featureCounts_sorted.bam
        samtools index -@ ${THREAD} ${i}_featureCounts_sorted.bam
        rm ${i}_Aligned.bam.featureCounts.bam
    fi

    if [ ${N6} == "F" ]; then
        if [ ${useUMITOOLS} == "T" ]; then
            umi_tools group \
                --verbose=0 \
                --extract-umi-method=read_id \
                --umi-separator=: \
                --per-cell --per-gene --gene-tag=XT \
                --method=cluster --edit-distance-threshold=1 \
                --skip-tags-regex=Unassigned \
                -I ${i}_featureCounts_sorted.bam \
                --output-bam -S ${i}_grouped.bam \
                --group-out=${i}_grouped.tsv
            cat ${i}_grouped.tsv | perl -F'\t' -alne 'next if($. == 1); if(s/:/\t/g){print;}' | \
                perl -alne 'BEGIN{$, = "\t";} $id = -1 if $.==1; if($id != $F[10]){$id = $F[10]; print($F[5], $F[1], $F[9]);}' | \
                perl -F'\t' -alne '
                    BEGIN{$, = "\t";} 
                    if($. == 1){$gene = $F[0]; $bc = $F[1]; $umis = 1; $reads = $F[2];}else{
                        if($F[0] eq $gene && $F[1] eq $bc){
                            $umis += 1; $reads += $F[2];
                        }else{
                            print $gene, $bc, $umis, $reads;
                            $gene = $F[0]; $bc = $F[1]; $umis = 1; $reads = $F[2];
                        }
                    }
                    END{print $gene, $bc, $umis, $reads;}
            ' > ${i}_counts.tsv
        else
            samtools view -@ ${THREAD} ${i}_featureCounts_sorted.bam | grep XT:Z: | perl -F'\t' -alne 'BEGIN{$,="\t";} $F[0] =~ s/:/\t/g; print @F;' | \
                perl -F'\t' -alne 'BEGIN{$,="\t";} print @F[0,1,2,4,5,$#F];'| sed s/XT:Z:// | sort -k6,6 -k4,4 -k5,5n > ${i}_assigned.tsv
            if [ ${allowMultiAssignment} == "F" ] && [ ${keepMultiMapping} == "T" ]; then
                python ${scriptDIR}/RNA/remove_multi_assignment.py -i ${i}_assigned.tsv -o ${i}_assigned_filtered.tsv
                python ${scriptDIR}/RNA/group_umi.py -i ${i}_assigned_filtered.tsv -o ${i}_grouped.tsv --max_dist 1
            else
                python ${scriptDIR}/RNA/group_umi.py -i ${i}_assigned.tsv -o ${i}_grouped.tsv --max_dist 1
            fi
            sort --parallel=${THREAD} -k1,1 -k2,2 ${i}_grouped.tsv | perl -F'\t' -alne'BEGIN{$, = "\t";}
                if($. == 1){$gene = $F[0]; $bc = $F[1]; $umis = 1; $reads = $F[3];}else{
                    if($F[0] eq $gene && $F[1] eq $bc){
                        $umis += 1; $reads += $F[3];
                    }else{
                        print $gene, $bc, $umis, $reads;
                        $gene = $F[0]; $bc = $F[1]; $umis = 1; $reads = $F[3];
                    }
                }
                END{print $gene, $bc, $umis, $reads;}
            ' > ${i}_counts.tsv
        fi
    else
        samtools view -@ ${THREAD} ${i}_featureCounts_sorted.bam | grep XT:Z: | perl -F'\t' -alne 'BEGIN{$,="\t";} $F[0] =~ s/:/\t/g; print @F;' | \
        perl -F'\t' -alne 'BEGIN{$,="\t";} print @F[0,1,2,4,5,$#F];'| sed s/XT:Z:// | sort -k6,6 -k4,4 -k5,5n > ${i}_assigned.tsv
        grep -v 'NNNNNNNNNN' ${i}_assigned.tsv > ${i}_assigned_dT.tsv
        grep 'NNNNNNNNNN' ${i}_assigned.tsv > ${i}_assigned_N6.tsv
        if [ ${allowMultiAssignment} == "F" ] && [ ${keepMultiMapping} == "T" ]; then
            python ${scriptDIR}/RNA/remove_multi_assignment.py -i ${i}_assigned_dT.tsv -o ${i}_assigned_filtered_dT.tsv
            mv ${i}_assigned_filtered_dT.tsv ${i}_assigned_dT.tsv
        fi
        # count dT umi by barcode-umi-cutting position
        python ${scriptDIR}/RNA/group_umi.py -i ${i}_assigned_dT.tsv -o ${i}_grouped_dT.tsv --max_dist 1
        # count N6 umi by barcode-cutting position. omit umi
        cat ${i}_assigned_N6.tsv | perl -alne 'BEGIN{$,="\t";} print @F[1,2,3,4,5];' | sort --parallel=${THREAD} -k3,3 -k4,4n -k1,1 -k5,5 | \
            uniq -c | perl -alne 'BEGIN{$,="\t";} print @F[5,1,2,0];' > ${i}_grouped_N6.tsv
        sort --parallel=${THREAD} -k1,1 -k2,2 ${i}_grouped_dT.tsv ${i}_grouped_N6.tsv | perl -F'\t' -alne'BEGIN{$, = "\t";}
            if($. == 1){$gene = $F[0]; $bc = $F[1]; $umis = 1; $reads = $F[3];}else{
                if($F[0] eq $gene && $F[1] eq $bc){
                    $umis += 1; $reads += $F[3];
                }else{
                    print $gene, $bc, $umis, $reads;
                    $gene = $F[0]; $bc = $F[1]; $umis = 1; $reads = $F[3];
                }
            }
            END{print $gene, $bc, $umis, $reads;}
        ' > ${i}_counts.tsv       

    fi



    # Output file: Gene Barcode Umis Reads
    if [ ${keepOneReadUMI} == "F" ]; then
        cat ${i}_counts.tsv | perl -alne 'print if $F[3] > 1;' > ${i}_rmoneread_counts.tsv
        mv ${i}_rmoneread_counts.tsv ${i}_counts.tsv
    fi

    # calculate PCR duplicate rate
    UNIQUE=`cat ${i}_counts.tsv | perl -alne '$total += $F[2]; END{print $total;}'`
    TOTAL=`cat ${i}_counts.tsv | perl -alne '$total += $F[3]; END{print $total;}'`
    perl -e '
        $unique = shift; $total = shift; $dup = $total - $unique;
        $prop = $dup/$total; $prop = sprintf("%.4f", $prop); $prop *= 100;
        print("Total reads after quntification: $total\n");
        print("Duplicate reads after quntification: $dup\n");
        print("Duplicate reads proportion: $prop%\n");
    ' ${UNIQUE} ${TOTAL}
    unset UNIQUE TOTAL

    # mito and ribo umi counts
    cat ${i}_counts.tsv | perl -alne '$umi += $F[2]; $mito += $F[2] if(/^mt-|^MT-/); $ribo += $F[2] if(/^RP[S|L]|^Rp[s|l]|rRNA/);
        END{$mito_prop = sprintf("%.4f", $mito/$umi); $mito_prop *= 100; $ribo_prop = sprintf("%.4f", $ribo/$umi); $ribo_prop *= 100;
            print "Total umis: $umi."; print "Mitochondria umis: $mito ($mito_prop%)."; print "Ribosome umis: $ribo ($ribo_prop%).";}'
    cat ${i}_counts.tsv | perl -alne '$read += $F[3]; $mito += $F[3] if(/^mt-|^MT-/); $ribo += $F[3] if(/^RP[S|L]|^Rp[s|l]|rRNA/);
        END{$mito_prop = sprintf("%.4f", $mito/$read); $mito_prop *= 100; $ribo_prop = sprintf("%.4f", $ribo/$read); $ribo_prop *= 100;
            print "Total reads: $read."; print "Mitochondria reads: $mito ($mito_prop%)."; print "Ribosome reads: $ribo ($ribo_prop%).";}'

    cat ${i}_counts.tsv | \
        perl -alne '
            BEGIN{$min = shift;}
            $umis{$F[1]} += $F[2];
            END{while(($k, $v) = each(%umis)){print "$k\t$v" if $v >= $min;}}' ${minUMI} | sort -k1,1 > ${i}_umis_count.txt
    cat ${i}_counts.tsv | \
        perl -alne '
            BEGIN{$min = shift;}
            $genes{$F[1]} += 1; $umis{$F[1]} += $F[2];
            END{while(($k, $v) = each(%genes)){print "$k\t$v" if $umis{$k} >= $min;}}' ${minUMI} | sort -k1,1 > ${i}_genes_count.txt 
    cat ${i}_counts.tsv | \
        perl -alne '
            BEGIN{$min = shift;}
            $reads{$F[1]} += $F[3]; $umis{$F[1]} += $F[2];
            END{while(($k, $v) = each(%reads)){print "$k\t$v" if $umis{$k} >= $min;}}' ${minUMI} | sort -k1,1 > ${i}_reads_count.txt

    python ${scriptDIR}/estimate_library_size.py \
        -r ${i}_reads_count.txt \
        -u ${i}_umis_count.txt \
        -m ${minUMI} \
        -o ${i}_estimate_libsize.txt

    python ${scriptDIR}/plot_distribution.py \
        -i ${i}_reads_count.txt \
        -o ${i}_reads_distribution.png \
        --field 1 --logx --logy \
        --start 1 --end 1000000 --otsu
    python ${scriptDIR}/plot_distribution.py \
        -i ${i}_umis_count.txt \
        -o ${i}_umis_distribution.png \
        --field 1 --logx --logy \
        --start 1 --end 1000000 --otsu
    python ${scriptDIR}/plot_distribution.py \
        -i ${i}_genes_count.txt \
        -o ${i}_genes_distribution.png \
        --field 1 --logx --logy \
        --start 1 --end 1000000 --otsu
    perl ${scriptDIR}/RNA/make_10X.pl \
        -g ${gtfDir} \
        -i ${i}_counts.tsv \
        -o ${i}_out.tar.gz \
        -m ${minUMI}



    ###################### Split Round 1 barcode tsv ######################
    echo "####################################################################################################################################"
    echo "spliting round 1 barcode..."
    echo "####################################################################################################################################"
 
     python ${scriptDIR}/split_round1.py \
        -i ${i}_counts.tsv \
        -u ${SPLIT_PATTERN} \
        -p ${i}_counts \
        -t tsv -l RNA

    for j in $(eval echo {1..${SPLIT_NUM}})
    do
        echo "Statistics for ${i}_${j}:"
        # calculate PCR duplicate rate

        UNIQUE=`cat ${i}_counts_part${j}.tsv | perl -alne '$total += $F[2]; END{print $total;}'`
        TOTAL=`cat ${i}_counts_part${j}.tsv | perl -alne '$total += $F[3]; END{print $total;}'`
        perl -e '
            $unique = shift; $total = shift; $dup = $total - $unique;
            $prop = $dup/$total; $prop = sprintf("%.4f", $prop); $prop *= 100;
            print("Total reads after quntification: $total\n");
            print("Duplicate reads after quntification: $dup\n");
            print("Duplicate reads proportion: $prop%\n");
        ' ${UNIQUE} ${TOTAL}
        unset UNIQUE TOTAL

        # mito and ribo umi counts
        cat ${i}_counts_part${j}.tsv | perl -alne '$umi += $F[2]; $mito += $F[2] if(/^mt-|^MT-/); $ribo += $F[2] if(/^RP[S|L]|^Rp[s|l]|rRNA/);
            END{$mito_prop = sprintf("%.4f", $mito/$umi); $mito_prop *= 100; $ribo_prop = sprintf("%.4f", $ribo/$umi); $ribo_prop *= 100;
                print "Total umis: $umi."; print "Mitochondria umis: $mito ($mito_prop%)."; print "Ribosome umis: $ribo ($ribo_prop%).";}'
        cat ${i}_counts_part${j}.tsv | perl -alne '$read += $F[3]; $mito += $F[3] if(/^mt-|^MT-/); $ribo += $F[3] if(/^RP[S|L]|^Rp[s|l]|rRNA/);
            END{$mito_prop = sprintf("%.4f", $mito/$read); $mito_prop *= 100; $ribo_prop = sprintf("%.4f", $ribo/$read); $ribo_prop *= 100;
                print "Total reads: $read."; print "Mitochondria reads: $mito ($mito_prop%)."; print "Ribosome reads: $ribo ($ribo_prop%).";}'

        cat ${i}_counts_part${j}.tsv | \
            perl -alne '
                BEGIN{$min = shift;}
                $umis{$F[1]} += $F[2];
                END{while(($k, $v) = each(%umis)){print "$k\t$v" if $v >= $min;}}' ${minUMI} | sort -k1,1 > ${i}_${j}_umis_count.txt        
        cat ${i}_counts_part${j}.tsv | \
            perl -alne '
                BEGIN{$min = shift;}
                $genes{$F[1]} += 1; $umis{$F[1]} += $F[2];
                END{while(($k, $v) = each(%genes)){print "$k\t$v" if $umis{$k} >= $min;}}' ${minUMI} | sort -k1,1 > ${i}_${j}_genes_count.txt   
        cat ${i}_counts_part${j}.tsv | \
            perl -alne '
                BEGIN{$min = shift;}
                $reads{$F[1]} += $F[3]; $umis{$F[1]} += $F[2];
                END{while(($k, $v) = each(%reads)){print "$k\t$v" if $umis{$k} >= $min;}}' ${minUMI} | sort -k1,1 > ${i}_${j}_reads_count.txt

        python ${scriptDIR}/estimate_library_size.py \
            -r ${i}_${j}_reads_count.txt \
            -u ${i}_${j}_umis_count.txt \
            -m ${minUMI} \
            -o ${i}_${j}_estimate_libsize.txt

        python ${scriptDIR}/plot_distribution.py \
            -i ${i}_${j}_reads_count.txt \
            -o ${i}_${j}_reads_distribution.png \
            --field 1 --logx --logy \
            --start 1 --end 1000000 --otsu
        python ${scriptDIR}/plot_distribution.py \
            -i ${i}_${j}_umis_count.txt \
            -o ${i}_${j}_umis_distribution.png \
            --field 1 --logx --logy \
            --start 1 --end 1000000 --otsu
        python ${scriptDIR}/plot_distribution.py \
            -i ${i}_${j}_genes_count.txt \
            -o ${i}_${j}_genes_distribution.png \
            --field 1 --logx --logy \
            --start 1 --end 1000000 --otsu
        perl ${scriptDIR}/RNA/make_10X.pl \
            -g ${gtfDir} \
            -i ${i}_counts_part${j}.tsv \
            -o ${i}_${j}_out.tar.gz \
            -m ${minUMI}
    done






