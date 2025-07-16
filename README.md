# TAI-seq #

## TAI-seq Pipeline ##
### System requirements ###
The pipeline relies on _Perl 5.32.1_, and therefore runs only on a UNIX-based system. We have tested the pipeline on _Ubuntu 22.04_ server in our lab.  
The ATAC pipeline requires >= 32 GB memory, RNA >= 64 GB, and TCR/BCR >= 256 GB. It is recommended to have at least 8 cores. We typically run the pipeline on a server with dual Epyc 7542 and 512 GB memory.

### Software dependencies ###
Dependencies  | Version
------------- | -------------
Perl | 5.32.1
python | 3.10.12
cutadapt | 4.5
bzip2 | 1.0.8
bedtools | 2.31.0
bowtie2 | 2.5.1
pigz | 2.6
samtools | 1.18
star | 2.7.11a
openjdk | 17.0.8
mixcr | 4.7.0

### Installation ###
Installation should finish within 10 minutes (Epyc 7542). 

(1) Download the *Pipeline* folder.  
(2) Setup the environment using conda:  
`conda create -n taiseq python=3.10.12 cutadapt=4.5 bzip2=1.0.8 bedtools=2.31.0 bowtie2=2.5.1 pigz=2.6 samtools=1.18 star=2.7.11a openjdk=17.0.8`  
Then  
`conda activate taiseq`  
(3) To run the TCR or BCR pipeline, you also need to install mixcr 4.7.0. Please install mixcr according to its instructions.

### Usage ###
#### RNA ####
(1) Setup parameters. Open `./Pipeline/RNA/scRNA.sh`, and set:  
`scriptDIR`: the _Pipeline_ folder  
`workingDIR`: data folder that contains fastq file folders (see demo).   
`gtfDir`: the genome gtf file   
`genomeDir`: STAR references  
`THREAD`: threads to run the pipeline  
`minUMI`: minimal umis for a barcode to be kept in the pipeline. To run the demo, set to 0. Otherwise, set to 50.  
`sample`: enter the sample names (separated by space if you have multiple samples). **Critical: the fastq files must be named as ${i}_R1.fq.gz and ${i}_R2.fq.gz, where i is your sample name. They must be put in a folder named ${i}, and the folder is placed in ${workingDIR}**  

(2) Additionally, if different samples are loaded to plate 1, you may set the sample demultiplexing parameters:  
`SPLIT_NUM`: how many samples are loaded to plate 1 (e.g., 2)  
`SPLIT_PATTERN`: the column counts of each sample loaded, separated by space (e.g., "4 8" indicates that the first sample is loaded to column 1-4, and the second sample is loaded to column 5-12)  

(3) Execute `bash /Pipeline/RNA/scRNA.sh` to run the TAI-RNA pipeline. You will get a ${i}_out.tar.gz file. Decompresse this file, and load into Seurat with _Read10X_.

#### ATAC ####
(1) Setup parameters. Open `./Pipeline/ATAC/scATAC.sh`, and set:  
`scriptDIR`: the _Pipeline_ folder  
`workingDIR`: data folder that contains fastq file folders (see demo).   
`genomeDir`: bowtie2 references  
`THREAD`: threads to run the pipeline     
`MINFRAGS`: minimal fragments for a barcode to be kept in the pipeline. To run the demo, set to 0. Otherwise, set to 50.  
`sample`: enter the sample names (separated by space if you have multiple samples). **Critical: the fastq files must be named as ${i}_R1.fq.gz and ${i}_R2.fq.gz, where i is your sample name. They must be put in a folder named ${i}, and the folder is placed in ${workingDIR}**  

(2) Additionally, if different samples are loaded to plate 1, you may set the sample demultiplexing parameters:  
`SPLIT_NUM`: how many samples are loaded to plate 1 (e.g., 2)  
`SPLIT_PATTERN`: the column counts of each sample loaded, separated by space (e.g., "4 8" indicates that the first sample is loaded to column 1-4, and the second sample is loaded to column 5-12)  

(3) Execute `bash /Pipeline/ATAC/scATAC.sh` to run the TAI-ATAC pipeline. You will get a ${i}_ArchR.tsv.gz file. This file can be read by ArchR.  

#### TCR ####
(1) Setup parameters. Open `./Pipeline/TCR/scTCR_circ_4.7.sh` or `scTCR_vmix_4.7.sh`, and set:  
`scriptDIR`: the _Pipeline_ folder  
`workingDIR`: data folder that contains fastq file folders (see demo).   
`SPECIES`: hsa for human, and mmu for mouse  
`THREAD`: threads to run the pipeline  
`sample`: enter the sample names (separated by space if you have multiple samples). **Critical: the fastq files must be named as ${i}_R1.fq.gz and ${i}_R2.fq.gz, where i is your sample name. They must be put in a folder named ${i}, and the folder is placed in ${workingDIR}**  

(2) Additionally, if different samples are loaded to plate 1, you may set the sample demultiplexing parameters:  
`SPLIT_PATTERN`: the column counts of each sample loaded, separated by space (e.g., "4 8" indicates that the first sample is loaded to column 1-4, and the second sample is loaded to column 5-12)  

(3) Execute `bash /Pipeline/TCR/scTCR_circ_4.7.sh` to run the TAI-TCR (full-length) pipeline. Execute `bash /Pipeline/TCR/scTCR_vmix_4.7.sh` to run the TAI-TCR (V primer mix) pipeline. The ${i}_mixcr_clones_umi.tsv and ${i}_mixcr_clones_noumi.tsv are the desired output files.  

### Run demo ###
Demos of all modelity should finish within 10 minutes.  
(1) Download the *Pipeline* folder and the *TAI-demo* folder.  
(2) Install the pipeline.  
(3) Run the pipeline according to the instructions above. Note: `workingDIR` should be set to the *TAI-demo* folder.  
You are expected to get test_rna_out.tar.gz for the RNA pipeline, test_ata_ArchR.tsv.gz for the ATAC pipeline, test_vmix_mixcr_clones_noumi.tsv + test_vmix_mixcr_clones_umi.tsv for the TAI-TCR-vmix pipeline, and test_circ_mixcr_clones_noumi.tsv + test_circ_mixcr_clones_umi.tsv for the TAI-TCR-circ pipeline.  





