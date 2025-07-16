#!/bin/bash

genome="mm10"
workingDIR="/data/xrz/TAI-demo"
scriptDIR="/data/xrz/Pipeline"
genomeDir="/data/xrz/ref/${genome}/${genome}-bowtie2/${genome}"
THREAD=24
MINFRAGS=0
SPLIT_NUM=2
SPLIT_PATTERN="6 6"
sample="test_ata"



for i in ${sample}; do
    echo "---------------------------------------------------------${i}------------------------------------------------------------------"
    cd ${workingDir}/${i}

    echo "Begin prcessing ${i}: `date`"
    if [[ ! -e ${i}_barcode_DNA_R1.fq.gz || ! -e ${i}_barcode_DNA_R2.fq.gz ]]; then
        # Allow 0 shift, 0 mismatch, v2 = Read45
        python ${scriptDIR}/process_fastq_phase.py \
        -1 ${i}_R1.fq.gz \
        -2 ${i}_R2.fq.gz \
        -q ${i}_barcode_DNA_R1.fq.gz \
        -e ${i}_barcode_DNA_R2.fq.gz \
        -d ATAC -c ${THREAD} -m 0 -s 0 -v 2
    fi

    # remove adaptor sequences. CRITICAL: cutadapt >= 4.1
    if [[ ! -e ${i}_barcode_cutadapt_DNA_R1.fq.gz || ! -e ${i}_barcode_cutadapt_DNA_R2.fq.gz ]]; then
        # # V1 (Read12)
        # cutadapt \
        #     -j ${THREAD} -m 2:2 \
        #     -a "GCAGCGTCAGATGTGTATAAGAGACAG;max_error_rate=0.4...CTGTCTCTTATACACATCT;optional;min_overlap=1;max_error_rate=0.2" \
        #     -A "AGATGTGTATAAGAGACAG;max_error_rate=0.2...CTGTCTCTTATACACATCT;optional;min_overlap=1;max_error_rate=0.2" \
        #     --pair-filter=any \
        #     -o ${i}_barcode_cutadapt_DNA_R1.fq.gz \
        #     -p ${i}_barcode_cutadapt_DNA_R2.fq.gz \
        #     ${i}_barcode_DNA_R1.fq.gz \
        #     ${i}_barcode_DNA_R2.fq.gz
        
        # V2 (Read45)
        cutadapt \
            -j ${THREAD} -m 2:2 \
            -a "AGATGTGTATAAGAGACAG;max_error_rate=0.2...CTGTCTCTTATACACATCT;optional;min_overlap=1;max_error_rate=0.2" \
            -A "AGATGTGTATAAGAGACAG;max_error_rate=0.2...CTGTCTCTTATACACATCT;optional;min_overlap=1;max_error_rate=0.2" \
            --pair-filter=any \
            -o ${i}_barcode_cutadapt_DNA_R1.fq.gz \
            -p ${i}_barcode_cutadapt_DNA_R2.fq.gz \
            ${i}_barcode_DNA_R1.fq.gz \
            ${i}_barcode_DNA_R2.fq.gz
    fi

    # if [[ ! -e ${i}_reads_per_barcode.txt ]]; then
    #     samtools import -@ ${THREAD} ${i}_barcode_DNA_R2.fq.gz | perl -alne 'BEGIN{$,="\t";} print $1 if($F[0] =~ /([A|G|C|T]{8}-[A|G|C|T]{8}-[A|G|C|T]{5}-{0,1}[A|G|C|T]{0,3})/);' | \
    #         sort --parallel=${THREAD} | uniq -c | perl -alne 'print "$F[1]\t$F[0]";' | sort --parallel=${THREAD} -k2,2nr > ${i}_reads_per_barcode.txt
    # fi    
    
    # python ${scriptDIR}/plot_distribution.py \
    #     -i ${i}_reads_per_barcode.txt \
    #     -o ${i}_reads_per_barcode.png \
    #     --field 1 --logx --logy \
    #     --start 1 --end 1000000  
        
    # mapping
    echo "Begin running bowtie2: `date`"
    bowtie2 \
        -p ${THREAD} -t -X2000 \
        --no-mixed --no-discordant \
        -x ${genomeDir} \
        -1 ${i}_barcode_cutadapt_DNA_R1.fq.gz \
        -2 ${i}_barcode_cutadapt_DNA_R2.fq.gz | perl -F'\t' -alne 'BEGIN{$|++;} if(/^@/){print;}else{chomp; if($F[0] =~ /:([A|G|C|T]{8}-[A|G|C|T]{8}-[A|G|C|T]{5})/){print "$_\tBC:Z:$1";}}' | \
        samtools sort -n -@ ${THREAD} > ${i}_Aligned.bam
    echo "bowtie2 finished: `date`"


    # sort by ID, retain paired q>30 reads, extract barcode as readID
    samtools view -h -@ ${THREAD} -f 2 -q 30 ${i}_Aligned.bam | \
        perl -alne 'BEGIN{$,="\t";} if(/^@/){print;}else{if($F[0]=~/.*:(.*)/){print $1,@F[1..@F];}}' | \
        samtools view -h -@ ${THREAD} -b > ${i}_sortedByID.bam

    # make fragment files, remove cell barcodes with less than ${MINFRAGS} fragments
    bedtools bamtobed -mate1 -bedpe -i ${i}_sortedByID.bam | perl -alne 'BEGIN{open POS, ">", shift; open NEG, ">", shift;}
        $F[8] eq "+" ? print POS : print NEG ; END{close POS; close NEG}' ${i}_R1_positive_strand.bedpe ${i}_R1_negative_strand.bedpe 
    cat ${i}_R1_positive_strand.bedpe | perl -alne 'BEGIN{$,="\t";} print($F[0],$F[1]+4,$F[5]-5,$F[6]);' | \
        sort --parallel=${THREAD} -k1,1 -k2,2n -k3,3n -k4,4 | uniq -c | perl -alne 'BEGIN{$,="\t";} print(@F[1..4],$F[0]) if($F[3]-$F[2] > 0);' > ${i}_ArchR_raw_positive.tsv
    cat ${i}_R1_negative_strand.bedpe | perl -alne 'BEGIN{$,="\t";} print($F[0],$F[4]+4,$F[2]-5,$F[6]);' | \
        sort --parallel=${THREAD} -k1,1 -k2,2n -k3,3n -k4,4 | uniq -c | perl -alne 'BEGIN{$,="\t";} print(@F[1..4],$F[0]) if($F[3]-$F[2] > 0);' > ${i}_ArchR_raw_negative.tsv
    cat ${i}_ArchR_raw_positive.tsv ${i}_ArchR_raw_negative.tsv | sort --parallel=${THREAD} -k1,1 -k2,2n -k3,3n -k4,4 | perl -F'\t' -alne ' BEGIN{$,="\t";}
        if($. == 1){$fragment = join("\t", @F[0..3]); $count = $F[4]; next;}
        if(join("\t", @F[0..3]) ne $fragment){
            print $fragment,$count; $fragment = join("\t", @F[0..3]); $count = $F[4];
        }else{
            $count += $F[4];
        }
        END{print $fragment,$count;} ' | sort --parallel=${THREAD} -k1,1 -k2,2n -k3,3n -k4,4  > ${i}_MACS2_raw.tsv
    cat ${i}_MACS2_raw.tsv | perl -F'\t' -alne ' BEGIN{$, = "\t"; $min = shift;}
        push @tmp, $_; $h{$F[3]} += 1;
        END{foreach(@tmp){@G = split /\t/; print($_) if($h{$G[3]} >= $min);}}' ${MINFRAGS} | \
        sort --parallel=${THREAD} -k1,1 -k2,2n -k3,3n -k4,4  > ${i}_ArchR_raw.tsv
    
    # # peak calling
    # grep -v -P "^chrM|^GL|^GH|^KI|^JH|^chrY|^MT" ${i}_MACS2_raw.tsv > ${i}_MACS2.tsv  
    # if [ ${genome} == "hg38" ]; then
    #     macs2 callpeak --keep-dup all --nomodel --nolambda --call-summits -f BEDPE -t ${i}_MACS2.tsv -n ${i} -g hs -p 0.05   
    # else
    #     macs2 callpeak --keep-dup all --nomodel --nolambda --call-summits -f BEDPE -t ${i}_MACS2.tsv -n ${i} -g mm -p 0.05
    # fi

    # calculate PCR duplicate rate
    # calculate mitochondria reads and fragments proportion
    cat ${i}_ArchR_raw.tsv | perl -F'\t' -alne '
        $read += $F[4] and $frag += 1 if(/^chrM|^MT/);
        $total_read += $F[4]; $total_frag += 1; 
        END{
            $prop_read = sprintf("%.4f", $read/$total_read); $prop_read *= 100;
            print("Total read count: $total_read"); print("Mitochondria read count: $read"); print("Mitochondria read proportion: $prop_read%");

            $prop_frag = sprintf("%.4f", $frag/$total_frag); $prop_frag *= 100; 
            print("Total fragment count: $total_frag"); print("Mitochondria fragment count: $frag"); print("Mitochondria fragment proportion: $prop_frag%");

            $dup = $total_read - $total_frag;
            $prop_dup = sprintf("%.4f", $dup/$total_read); $prop_dup *= 100;
            print("Total reads with q>30: $total_read"); print("Duplicate reads with q>30: $dup"); print("Duplicate reads proportion: $prop_dup%");
        }
    '

    # calculate library unique fragment length distribution 
    cat ${i}_ArchR_raw.tsv | perl -alne '$f_length = $F[2]-$F[1]; print $f_length;' | sort -k1,1n | uniq -c | \
        perl -alne 'BEGIN{$,="\t"; print("Length\tCounts");} print($F[1],$F[0]);' > ${i}_fragment_distribution.txt  
 
    # count unique fragments
    grep -v -P "^chrM|^MT" ${i}_ArchR_raw.tsv | perl -F'\t' -alne '
        $h{$F[3]} += 1;
        END{while(($k, $v) = each(%h)){print("$k\t$v");}}' > ${i}_fragments_count.txt
    grep -v -P "^chrM|^MT" ${i}_ArchR_raw.tsv | perl -F'\t' -alne '
        $h{$F[3]} += $F[4];
        END{while(($k, $v) = each(%h)){print("$k\t$v");}}' > ${i}_reads_count.txt

    python ${scriptDIR}/estimate_library_size.py \
        -r ${i}_reads_count.txt \
        -u ${i}_fragments_count.txt \
        -m ${MINFRAGS} \
        -o ${i}_estimate_libsize.txt

    grep -v -P "^chrM|^GL|^GH|^KI|^JH|^chrY|^MT" ${i}_ArchR_raw.tsv > ${i}_ArchR.tsv
    cat ${i}_ArchR.tsv | perl -F'\t' -alne '$total += 1; END{print("Total unique fragments after chromosome filtering: $total");}'

    # compress to gz file
    bgzip -f ${i}_ArchR.tsv 

    python ${scriptDIR}/plot_distribution.py \
        -i ${i}_reads_count.txt \
        -o ${i}_reads_distribution.png \
        --field 1 --logx --logy \
        --start 1 --end 1000000 --otsu
    python ${scriptDIR}/plot_distribution.py \
        -i ${i}_fragments_count.txt \
        -o ${i}_fragments_distribution.png \
        --field 1 --logx --logy \
        --start 1 --end 1000000 --otsu


    





    ###################### Split Round 1 barcode tsv ######################
    echo "####################################################################################################################################"
    echo "spliting round 1 barcode..."
    echo "####################################################################################################################################"

    python ${scriptDIR}/split_round1.py \
        -i ${i}_ArchR_raw.tsv \
        -u ${SPLIT_PATTERN} \
        -p ${i}_ArchR_raw \
        -t tsv -l ATAC

    for j in $(eval echo {1..${SPLIT_NUM}})
    do
        echo "Statistics for ${i}_${j}:"
        # macs2 callpeak --keep-dup all --nomodel --nolambda --call-summits -f BEDPE -t ${i}_MACS2_part${j}.tsv -n ${i}_${j} -g hs -p 0.05 
        cat ${i}_ArchR_raw_part${j}.tsv | perl -F'\t' -alne '
            $read += $F[4] and $frag += 1 if(/^chrM|^MT/);
            $total_read += $F[4]; $total_frag += 1; 
            END{
                $prop_read = sprintf("%.4f", $read/$total_read); $prop_read *= 100;
                print("Total read count: $total_read"); print("Mitochondria read count: $read"); print("Mitochondria read proportion: $prop_read%");

                $prop_frag = sprintf("%.4f", $frag/$total_frag); $prop_frag *= 100; 
                print("Total fragment count: $total_frag"); print("Mitochondria fragment count: $frag"); print("Mitochondria fragment proportion: $prop_frag%");

                $dup = $total_read - $total_frag;
                $prop_dup = sprintf("%.4f", $dup/$total_read); $prop_dup *= 100;
                print("Total reads with q>30: $total_read"); print("Duplicate reads with q>30: $dup"); print("Duplicate reads proportion: $prop_dup%");
            }
        '

        # calculate library unique fragment length distribution 
        cat ${i}_ArchR_raw_part${j}.tsv | perl -alne '$f_length = $F[2]-$F[1]; print $f_length;' | sort -k1,1n | uniq -c | \
            perl -alne 'BEGIN{$,="\t"; print("Length\tCounts");} print($F[1],$F[0]);' > ${i}_${j}_fragment_distribution.txt  
    
        # count unique fragments
        grep -v -P "^chrM|^MT" ${i}_ArchR_raw_part${j}.tsv | perl -F'\t' -alne '
            $h{$F[3]} += 1;
            END{while(($k, $v) = each(%h)){print("$k\t$v");}}' > ${i}_${j}_fragments_count.txt
        grep -v -P "^chrM|^MT" ${i}_ArchR_raw_part${j}.tsv | perl -F'\t' -alne '
            $h{$F[3]} += $F[4];
            END{while(($k, $v) = each(%h)){print("$k\t$v");}}' > ${i}_${j}_reads_count.txt

        python ${scriptDIR}/estimate_library_size.py \
            -r ${i}_${j}_reads_count.txt \
            -u ${i}_${j}_fragments_count.txt \
            -m ${MINFRAGS} \
            -o ${i}_${j}_estimate_libsize.txt
            
        grep -v -P "^chrM|^GL|^GH|^KI|^JH|^chrY|^MT" ${i}_ArchR_raw_part${j}.tsv > ${i}_${j}_ArchR.tsv
        cat ${i}_${j}_ArchR.tsv | perl -F'\t' -alne '$total += 1; END{print("Total unique fragments after chromosome filtering: $total");}'

        # compress to gz file
        bgzip -f ${i}_${j}_ArchR.tsv 

        python ${scriptDIR}/plot_distribution.py \
            -i ${i}_${j}_reads_count.txt \
            -o ${i}_${j}_reads_distribution.png \
            --field 1 --logx --logy \
            --start 1 --end 1000000
        python ${scriptDIR}/plot_distribution.py \
            -i ${i}_${j}_fragments_count.txt \
            -o ${i}_${j}_fragments_distribution.png \
            --field 1 --logx --logy \
            --start 1 --end 1000000
    done


    zip ${i}_out.zip *.log *.txt *.png *.tsv.gz

    echo "Finish prcessing ${i}: `date`\n"

done 


