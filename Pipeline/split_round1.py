import argparse
import re
import os


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in', dest='IN', type=str, action='store', help='input bed file to split')
parser.add_argument('-s', '--split', dest='SPLIT', type=int, default=None, action='store', help='number of files to split into')
parser.add_argument('-u', '--sub', dest='SUB', nargs='+', type=int, default=None, action='store', help='tuple of column numbers of each sublibrary')
parser.add_argument('-p', '--prefix', dest='PREFIX', type=str, action='store', help='prefix of output files')
parser.add_argument('-t', '--type', dest='FTYPE', default='tsv', type=str, action='store', help='input file type')
parser.add_argument('-l', '--library', dest='LIB', default='ATAC', type=str, action='store', help='library type')


args = parser.parse_args()
infile = args.IN
number_of_split = args.SPLIT
prefix = args.PREFIX
file_type = args.FTYPE
lib_type = args.LIB
sublib = args.SUB
if sublib is not None:
    number_of_split = len(sublib)



round1_bc = [
    "AACGTGAT", "AAACATCG", "ATGCCTAA", "AGTGGTCA", "ACCACTGT", "ACATTGGC", "CAGATCTG", "CATCAAGT", "CGCTGATC", "ACAAGCTA",
    "CTGTAGCC", "AGTACAAG", "AACAACCA", "AACCGAGA", "AACGCTTA", "AAGACGGA", "AAGGTACA", "ACACAGAA", "ACAGCAGA", "ACCTCCAA",
    "ACGCTCGA", "ACGTATCA", "ACTATGCA", "AGAGTCAA", "AGATCGCA", "AGCAGGAA", "AGTCACTA", "ATCCTGTA", "ATTGAGGA", "CAACCACA",
    "GACTAGTA", "CAATGGAA", "CACTTCGA", "CAGCGTTA", "CATACCAA", "CCAGTTCA", "CCGAAGTA", "CCGTGAGA", "CCTCCTGA", "CGAACTTA",
    "CGACTGGA", "CGCATACA", "CTCAATGA", "CTGAGCCA", "CTGGCATA", "GAATCTGA", "CAAGACTA", "GAGCTGAA", "GATAGACA", "GCCACATA",
    "GCGAGTAA", "GCTAACGA", "GCTCGGTA", "GGAGAACA", "GGTGCGAA", "GTACGCAA", "GTCGTAGA", "GTCTGTCA", "GTGTTCTA", "TAGGATGA",
    "TATCAGCA", "TCCGTCTA", "TCTTCACA", "TGAAGAGA", "TGGAACAA", "TGGCTTCA", "TGGTGGTA", "TTCACGCA", "AACTCACC", "AAGAGATC",
    "AAGGACAC", "AATCCGTC", "AATGTTGC", "ACACGACC", "ACAGATTC", "AGATGTAC", "AGCACCTC", "AGCCATGC", "AGGCTAAC", "ATAGCGAC",
    "ATCATTCC", "ATTGGCTC", "CAAGGAGC", "CACCTTAC", "CCATCCTC", "CCGACAAC", "CCTAATCC", "CCTCTATC", "CGACACAC", "CGGATTGC",
    "CTAAGGTC", "GAACAGGC", "GACAGTGC", "GAGTTAGC", "GATGAATC", "GCCAAGAC"
]

def split_reference_barcode(reference, number_of_split: int) -> list:
    number_of_elements = int(96/number_of_split)
    groups = []
    for i in range(number_of_split):
        groups.append(reference[i*number_of_elements:(i+1)*number_of_elements])
    return groups

def split_reference_barcode_provided(reference: list, sublib: list) -> list:
    ref_used = reference.copy()
    number_of_split = len(sublib)
    groups = []
    for i in range(number_of_split):
        groups.append(ref_used[:sublib[i]*8])
        del ref_used[:sublib[i]*8]
    return groups

def open_output_files_bed(number_of_split, prefix):
    out = []
    for i in range(number_of_split):
        f = open(prefix + '_part' + str(i+1) + '.' + file_type, 'w')
        out.append(f)
    return(out)

def open_output_files_sam(number_of_split, prefix):
    out = []
    for i in range(number_of_split):
        f = open(prefix + '_part' + str(i+1) + '.sam', 'w')
        out.append(f)
    return(out)

def close_output_files(files):
    for f in files:
        f.close()

def extract_round1_barcode_bed(cols, lib_type):
    if lib_type == "ATAC":
        barcode = cols[3] # CTGAGCCA-GGAGAACA-TCGAC
    elif lib_type == "RNA":
        barcode = cols[1] # CTGAGCCA-GGAGAACA-TCGAC-GCA
    elif lib_type in ["TCR", "BCR"]:
        barcode = cols[-1]
    barcode_to_match = barcode[0:8]
    return barcode_to_match

def extract_round1_barcode_bam(cols):
    pattern = re.compile('[A|G|C|T]{8}-[A|G|C|T]{8}-[A|G|C|T]{5}')
    barcode = re.search(pattern, cols[0])[0]
    barcode_to_match = barcode[0:8]
    return barcode_to_match

if __name__ == "__main__":
    if sublib is not None:
        groups = split_reference_barcode_provided(round1_bc, sublib)
    else:
        groups = split_reference_barcode(round1_bc, number_of_split)
    if file_type in ["bed", "tsv"]:
        out_files = open_output_files_bed(number_of_split, prefix)
        with open(infile, "r") as f:
            for line in f:
                cols = line.strip().split('\t')
                if cols[-1] == "cellId":
                    for o in out_files:
                        o.write(line)                
                barcode_to_match = extract_round1_barcode_bed(cols, lib_type)
                for i, group in enumerate(groups):
                    if barcode_to_match in group:
                        out_files[i].write(line)
                        break
        close_output_files(out_files)
    elif file_type == "bam":
        out_files = open_output_files_sam(number_of_split, prefix)
        f = os.popen("samtools view -h " + infile)
        for line in f:
            if line.startswith("@"):
                for out in out_files:
                    out.write(line)
            else:
                cols = line.strip().split('\t')
                barcode_to_match = extract_round1_barcode_bam(cols)
                for i, group in enumerate(groups):
                    if barcode_to_match in group:
                        out_files[i].write(line)
                        break  
        close_output_files(out_files)      
        f.close()

        

        
