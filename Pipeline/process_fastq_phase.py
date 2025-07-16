from multiprocessing import Queue, Process, Pipe, connection
import argparse
import Levenshtein
import xopen
import dnaio
import time
import traceback
import io
from typing import Optional, BinaryIO, Tuple

# Author: Xie Runze #
# 2023/05/10 #
# This script extracts cell barcodes of each read.
# This script is for JoINT-seq (RNA + ATAC + TCR + BCR). RNA library is created with template switching stragety.

parser = argparse.ArgumentParser()
parser.add_argument('-1', '--r1', dest='R1', type=str, action='store', help='R1.fq.gz')
parser.add_argument('-2', '--r2', dest='R2', type=str, action='store', help='R2.fq.gz')
parser.add_argument('-q', '--o1', dest='O1', default=None, type=str, action='store', help='O1.fq.gz')
parser.add_argument('-e', '--o2', dest='O2', default=None, type=str, action='store', help='O2.fq.gz')
parser.add_argument('-m', '--mismatch', dest='MISMATCH', default=0, type=int, action='store', help='maximum mismatches in oligo dT')
parser.add_argument('-d', '--mode', dest='MODE', default="ATAC", type=str, action='store', help='Processed library type')
parser.add_argument('-s', '--shift', dest='SHIFT', default=0, type=int, action='store', help='maximum shift allowed when matching barcodes')
parser.add_argument('-n', '--number_of_dT', dest='NUM', default=8, type=int, action='store', help='number of oligo dT to match')
parser.add_argument('-c', '--cores', dest='CORES', default=6, type=int, action='store', help='number of threads to use')
parser.add_argument('-p', '--start_pos_dT', dest='START_POS_DT', default=142, type=int, action='store', help='starting position of oligo dT') # 123 for old, 142 for V1, 138 for V2
parser.add_argument('-v', '--primer_version', dest='PRIMER_VERSION', default=1, type=int, action='store', help='primer version. 1 = read12; 2 = read45') # 123, 142


args = parser.parse_args()

read1 = args.R1
read2 = args.R2
out1 = args.O1
out2 = args.O2
mode = args.MODE
mismatch = args.MISMATCH
T_num = args.NUM
match = 'T'*T_num
shift = args.SHIFT
cores = args.CORES
dt_start = args.START_POS_DT
pversion = args.PRIMER_VERSION

if pversion == 1:
    anchor = 'GTCTCGTGGGCTCGG'
elif pversion == 2:
    anchor = 'TGAAGCGACTGCTGT'

tr = str.maketrans('ATCGN', 'TAGCN')
# # if old library structure without sample barcode
if mode == "RNA_old":
    dt_start = 123

whitelist_bc1 = [
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
whitelist_bc2 = [
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
] # ACAGCAGT CGCTGATC
whitelist_bc3 = [
    "AACCA", "AACGG", "AACTC", "AAGAG", "AAGCC", "AAGGA", "AATCG", "AATGC", "ACACA", "ACAGG",
    "ACATC", "ACCAA", "ACCTG", "ACGAT", "ACGTA", "ACTAC", "ACTCT", "ACTGA", "AGAAG", "AGACC",
    "AGAGA", "AGCAC", "AGCGT", "AGCTA", "AGGAA", "AGGCT", "AGGTG", "AGTCA", "AGTGG", "AGTTC",
    "ATACG", "ATAGC", "ATCAG", "ATCCT", "ATCGA", "CAACA", "CAAGG", "CACAT", "CACGA", "CACTG",
    "CAGAC", "CAGCT", "CAGTA", "CATAG", "CATCC", "CATGT", "CCAAC", "CCAGA", "CCATG", "CCGAA",
    "CCTAT", "CCTCA", "CCTTC", "CGAAT", "CGATA", "GAACG", "GAAGA", "GAATC", "GACAG", "GACCT",
    "GACTA", "GAGAA", "GAGGT", "GAGTG", "GATAC", "GATCA", "GATGG", "GCAAG", "GCACT", "GCATA",
    "GCCAT", "GCGTT", "GCTAA", "GCTGT", "GCTTG", "GGAAC", "TAAGC", "TACAC", "TACCG", "TACGT",
    "TAGCA", "TAGTC", "TCACG", "TCAGT", "TCCAG", "TCCTA", "TCGAC", "TCGCT", "TCGGA", "TCGTG",
    "TCTCC", "TCTGG", "TGACT", "TGAGG", "TGATC", "TGCAT"
]
whitelist_sample = [
    "ACT", "GCA", "TCC", "CAT"
]

class InputFiles:
    def __init__(self, file1: BinaryIO, file2: Optional[BinaryIO] = None, interleaved: bool = False):
        self.file1 = file1
        self.file2 = file2
        self.interleaved = interleaved

    def open(self):
        return dnaio.open(self.file1, file2 = self.file2, interleaved=self.interleaved, mode="r")

    def close(self) -> None:
        self.file1.close()
        if self.file2 is not None:
            self.file2.close()


class OutputFiles:
    def __init__(
        self,
        out1: Optional[BinaryIO] = None,
        out2: Optional[BinaryIO] = None,

    ):
        self.out1 = out1
        self.out2 = out2
    def __iter__(self):
        for f in [self.out1, self.out2]:
            if f is not None:
                yield f

    def as_bytesio(self) -> "OutputFiles":
        """
        Create a new OutputFiles instance that has BytesIO instances for each non-None output file
        """
        result = OutputFiles()
        for attr in ("out1", "out2"):
            if getattr(self, attr) is not None:
                setattr(result, attr, io.BytesIO())
        return result


def open_input_files(file1, file2 = None) -> Tuple:
    r1 = xopen.xopen(file1, "rb")
    if file2 is not None:
        r2 = xopen.xopen(file2, "rb")
    return r1, r2

def open_output_files(mode, out1 = None, out2 = None) -> OutputFiles:
    """
    Select processing mode and open output file handles. Valid mode includes ATAC, RNA, TCR and BCR.
    """
    o1, o2, = None, None
    if mode in ["TCR", "TCR_v2", "TCR_v2_filt_untag"]:
        print("Barcoding TCR library...")
        o1 = xopen.xopen(out1, 'wb', compresslevel=6)
        o2 = xopen.xopen(out2, 'wb', compresslevel=6)
    elif mode in ["BCR", "BCR_v2", "BCR_v2_filt_untag"]:
        print("Barcoding BCR library...")
        o1 = xopen.xopen(out1, 'wb', compresslevel=6)
        o2 = xopen.xopen(out2, 'wb', compresslevel=6)
    elif mode in ["Vmix", "Vmix_v2"]:
        print("Barcoding JoINT Vmix amplicon TCR library...")
        o1 = xopen.xopen(out1, 'wb', compresslevel=6)
        o2 = xopen.xopen(out2, 'wb', compresslevel=6)
    elif mode in ["RNA", "RNA_old", "RNA_mixcr", "RNA_nolibbc", "RNA_n6", "RNA_n6_test"]:
        print("Barcoding RNA library...")
        o1 = xopen.xopen(out1, 'wb', compresslevel=6)
        o2 = xopen.xopen(out2, 'wb', compresslevel=6)
    elif mode == "ATAC":
        print("Barcoding ATAC library...")
        o1 = xopen.xopen(out1, 'wb', compresslevel=6)
        o2 = xopen.xopen(out2, 'wb', compresslevel=6)
    else:
        print("Barcoding unknown library...")
        o1 = xopen.xopen(out1, 'wb', compresslevel=6)
        o2 = xopen.xopen(out2, 'wb', compresslevel=6)
    return OutputFiles(o1, o2)

def find_phase(seq, mismatch_anchor = 3, weights = (15,15,1), anchor = anchor):
    # anchor = "GTCTCGTGGGCTCGG" # 96 111
    # anchor = 'TGAAGCGACTGCTGT' # 96 111
    # anchor = "GTGGCCGATGTTTCG" # 5 20
    for phase in range(4):
        if Levenshtein.distance(anchor, seq[96+phase:111+phase], weights = weights) <= mismatch_anchor:
        # if Levenshtein.distance(anchor, seq[5+phase:20+phase], weights = weights) <= mismatch_anchor:
            return phase
    return -1


def find_phase_immune_repertoire(seq, mismatch_anchor = 3, weights = (15,15,1)):
    anchor_immune = "CCATCTGAGCCACCA"
    for phase in range(4):
        if Levenshtein.distance(anchor_immune, seq[0+phase:15+phase], weights = weights) <= mismatch_anchor:
            return phase
    return -1


def assign_barcode(seq, round, mismatch, whitelist = [whitelist_bc3, whitelist_bc2, whitelist_bc1, whitelist_sample], shift = shift) -> Tuple:
    """
    Assign the barcode of each sequence. Allow 0 mismatch and 0 shift compared to the barcode whitelist.
    Return: (position shift, valid barcode)
    Will return (0,'') if no barcode is assigned.
    """
    pos_init = [(0,5), (35,43), (73,81), (121,124)]
    pos_shift = range(-shift, shift+1, 1)
    for valid_bc in whitelist[round-1]:
        if round == 2:
            mismatch = 1
        for i, to_check in enumerate(generate_checklist(pos_init[round-1], shift, seq)):
            if Levenshtein.distance(valid_bc, to_check) <= mismatch:
                return pos_shift[i], valid_bc
            else:
                pass
    return 0, ''

def assign_barcode_immune_repertoire(seq, round, mismatch, whitelist = [whitelist_sample, whitelist_bc1, whitelist_bc2, whitelist_bc3], shift = shift) -> Tuple:
    """
    Assign the barcode of each sequence. Allow 0 mismatch and 0 shift compared to the barcode whitelist.
    Return: (position shift, valid barcode)
    Will return (0,'') if no barcode is assigned.
    """
    # sample barcode: [19,21] | UMI [22,31] | barcode1 [62,69] | barcode2 [100,107] | barcode3 [138,142]
    pos_init = [(18,21), (61,69), (99,107), (137,142)]
    pos_shift = range(-shift, shift+1, 1)
    for valid_bc in whitelist[round-1]:
        for i, to_check in enumerate(generate_checklist(pos_init[round-1], shift, seq)):
            if Levenshtein.distance(valid_bc, reverse_complement(to_check, tr)) <= mismatch:
                return pos_shift[i], valid_bc
            else:
                pass
    return 0, ''

def assign_barcode_immune_repertoire_v2(seq, round, mismatch, whitelist = [whitelist_bc1, whitelist_bc2, whitelist_bc3], shift = shift) -> Tuple:
    """
    Assign the barcode of each sequence. Allow 0 mismatch and 0 shift compared to the barcode whitelist.
    Return: (position shift, valid barcode)
    Will return (0,'') if no barcode is assigned.
    """
    # UMI [19,28] | barcode1 [59,66] | barcode2 [97,104] | barcode3 [135,139]
    pos_init = [(56,64), (94,102), (132,137)]
    pos_shift = range(-shift, shift+1, 1)
    for valid_bc in whitelist[round-1]:
        for i, to_check in enumerate(generate_checklist(pos_init[round-1], shift, seq)):
            if Levenshtein.distance(valid_bc, reverse_complement(to_check, tr)) <= mismatch:
                return pos_shift[i], valid_bc
            else:
                pass
    return 0, ''


def generate_checklist(pos, shift, seq) -> list:
    out = list()
    for i in range(-shift, shift+1, 1):
        out.append(seq[(pos[0]+i):(pos[1]+i)])
    return out

def reverse_complement(seq: str, tr) -> str:
    return seq.translate(tr)[::-1]

def calculate_statistics(c1, c2, total, elapsed_time, mode):
    c1 = c1 * 2
    c2 = c2 * 2
    total = total * 2
    formatted_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    print("Total time elapsed running {}: {}".format(__file__, formatted_time))
    print("Total processed reads: {}".format(total))
    if mode in ["BCR", "BCR_v2", "BCR_v2_filt_untag"]:
        print("BCR read counts: {}". format(c2))
        print("BCR read proportion: {:.2%}". format(c2/total))
        print("Dropped read counts: {}". format(total-c2))
        print("Dropped read proportion: {:.2%}". format((total-c2)/total))
    elif mode in ["TCR", "TCR_v2", "TCR_v2_filt_untag"]:
        print("TCR read counts: {}". format(c2))
        print("TCR read proportion: {:.2%}". format(c2/total))
        print("Dropped read counts: {}". format(total-c2))
        print("Dropped read proportion: {:.2%}". format((total-c2)/total))
    elif mode == "ATAC":
        print("ATAC read counts: {}". format(c1))
        print("ATAC read proportion: {:.2%}". format(c1/total))
        print("Dropped read counts: {}". format(total-c1))
        print("Dropped read proportion: {:.2%}". format((total-c1)/total))
    elif mode in ["RNA", "RNA_old", "RNA_mixcr", "RNA_nolibbc"]:
        print("RNA read counts: {}". format(c1))
        print("RNA read proportion: {:.2%}". format(c1/total))
        print("Dropped read counts: {}". format(total-c1))
        print("Dropped read proportion: {:.2%}". format((total-c1)/total))
    elif mode in ["Vmix", "Vmix_v2"]:
        print("TCR read counts: {}". format(c1))
        print("TCR read proportion: {:.2%}". format(c1/total))
        print("Dropped read counts: {}". format(total-c1))
        print("Dropped read proportion: {:.2%}". format((total-c1)/total))
    elif mode in ["RNA_n6", "RNA_n6_test"]:
        print("RNA oligo dT read counts: {}". format(c1))
        print("RNA oligo dT read proportion: {:.2%}". format(c1/total))
        print("RNA random hexamer read counts: {}". format(c2))
        print("RNA random hexamer read proportion: {:.2%}". format(c2/total))
        print("Dropped read counts: {}". format(total-c1-c2))
        print("Dropped read proportion: {:.2%}". format((total-c1-c2)/total))        
    else:
        print("Valid read counts: {}". format(c1))
        print("Valid read counts: {:.2%}". format(c1/total))
        print("Dropped read counts: {}". format(total-c1))
        print("Dropped read proportion: {:.2%}". format((total-c1)/total))
    return None


class FileReader(Process):
    def __init__(self, handles, queue, conn_list) -> None:
        super(FileReader, self).__init__()
        self._conn_list = conn_list
        self.r1, self.r2 = handles
        self._q = queue

    def acquire_some_lines(self):
        """
        Read fastq files and send 4 lines to DataProcesser. Whenever working_queue has content, 
        it means a DataProcesser is idel and we send data to this DataProcesser.
        """
        index = 0
        for idx, (chunk1, chunk2) in enumerate(dnaio.read_paired_chunks(self.r1, self.r2, buffer_size=4*1024**2)):
            self.send_to_worker(idx, chunk1, chunk2)
        return index

    def send_to_worker(self, chunk_index, chunk1, chunk2 = None):
        """
        Get and idel worker index and send data to the worker.
        """
        worker_index = self._q.get()
        connection = self._conn_list[worker_index]
        connection.send(chunk_index)
        connection.send_bytes(chunk1)
        if chunk2 is not None:
            connection.send_bytes(chunk2)
        return None

    def _shutdown(self):
        """
        If reader process finished, send poison pills to all workers.
        """
        # Send poison pills to all workers
        for _ in range(len(self._conn_list)):
            worker_index = self._q.get()
            self._conn_list[worker_index].send(-1)
        return None

    def run(self) -> int:
        index = self.acquire_some_lines()
        self._shutdown()
        return index


class DataProcesser(Process):
    def __init__(self, index, queue, mode, conn_in, conn_out, original_outfiles: OutputFiles = None) -> None:
        super(DataProcesser, self).__init__()
        self._q = queue
        self._conn_in = conn_in
        self._conn_out = conn_out
        self._id = index
        self._mode = mode
        self._original_outfiles = original_outfiles
        self._reader = None        
    
    def process_read(self, read1, read2):
        """
        Receive data: dnaio.Sequence (R1), dnaio.Sequence (R2)
        Return data: flag, dnaio.Sequence (R1), dnaio.Sequence (R2)
        Will return flag = 0 if the sequence is dropped.

        flag == 0 -> dropped
        flag == 1 -> RNA
        flag == 2 -> TCR/BCR
        flag == 3 -> ATAC
        """
        flag = 0
        if self._mode == "RNA_nolibbc":
            phase = find_phase(read2.sequence)
            if phase == -1:
                flag = 0
                return flag, None, None
            shifted_seq = read2.sequence[phase:] + 'T'*phase
            dT = shifted_seq[dt_start:dt_start+T_num]
            if Levenshtein.distance(dT, match) <= 2:
                # then this sequence is RNA
                flag = 1
                shift3, bc3 = assign_barcode(shifted_seq, round=1, mismatch=mismatch)
                shift2, bc2 = assign_barcode(shifted_seq, round=2, mismatch=mismatch)
                shift1, bc1 = assign_barcode(shifted_seq, round=3, mismatch=mismatch)
                # shift4, lib_bc = assign_barcode(shifted_seq, round=4, mismatch=mismatch)
                if pversion == 1:
                    umi = shifted_seq[111+shift3:121+shift3]
                elif pversion == 2:
                    umi = shifted_seq[111+shift3:119+shift3]
                header1 = read1.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/1'
                read1.name = header1
                header2 = read2.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/2'
                read2.name = header2    
                if bc1 and bc2 and bc3:
                    return flag, read1, read2
                else:    
                    flag = 0
                    return flag, None, None  
            else: 
                flag = 0
                return flag, None, None  
        elif self._mode == "RNA_n6":
            phase = find_phase(read2.sequence)
            if phase == -1:
                flag = 0
                return flag, None, None
            shifted_seq = read2.sequence[phase:] + 'T'*phase
            dT = shifted_seq[dt_start:dt_start+T_num]
            if Levenshtein.distance(dT, match) <= 2:
                # then this sequence is oligo dT RNA
                flag = 1
                shift3, bc3 = assign_barcode(shifted_seq, round=1, mismatch=mismatch)
                shift2, bc2 = assign_barcode(shifted_seq, round=2, mismatch=mismatch)
                shift1, bc1 = assign_barcode(shifted_seq, round=3, mismatch=mismatch)
                # shift4, lib_bc = assign_barcode(shifted_seq, round=4, mismatch=mismatch)
                umi = shifted_seq[111+shift3:121+shift3]
                header1 = read1.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/1'
                read1.name = header1
                header2 = read2.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/2'
                read2.name = header2    
                if bc1 and bc2 and bc3:
                    return flag, read1, read2
                else:    
                    flag = 0
                    return flag, None, None  
            else:
                # then this sequence is random N6 RNA
                flag = 2
                shift3, bc3 = assign_barcode(shifted_seq, round=1, mismatch=mismatch)
                shift2, bc2 = assign_barcode(shifted_seq, round=2, mismatch=mismatch)
                shift1, bc1 = assign_barcode(shifted_seq, round=3, mismatch=mismatch)
                # shift4, lib_bc = assign_barcode(shifted_seq, round=4, mismatch=mismatch)
                umi = 'NNNNNNNNNN'
                header1 = read1.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/1'
                read1.name = header1
                header2 = read2.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/2'
                read2.name = header2    
                if bc1 and bc2 and bc3:
                    return flag, read1, read2
                else:    
                    flag = 0
                    return flag, None, None  
        elif self._mode == "RNA_n6_test":
            phase = find_phase(read2.sequence)
            if phase == -1:
                flag = 0
                return flag, None, None
            shifted_seq = read2.sequence[phase:] + 'T'*phase
            dT = shifted_seq[dt_start:dt_start+T_num]
            if Levenshtein.distance(dT, match) <= 2:
                # then this sequence is oligo dT RNA
                flag = 1
                shift3, bc3 = assign_barcode(shifted_seq, round=1, mismatch=mismatch)
                shift2, bc2 = assign_barcode(shifted_seq, round=2, mismatch=mismatch)
                shift1, bc1 = assign_barcode(shifted_seq, round=3, mismatch=mismatch)
                # shift4, lib_bc = assign_barcode(shifted_seq, round=4, mismatch=mismatch)
                umi = shifted_seq[111+shift3:121+shift3]
                header1 = read1.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/1'
                read1.name = header1
                header2 = read2.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/2'
                read2.name = header2    
                if bc1 and bc2 and bc3:
                    return flag, read1, read2
                else:    
                    flag = 0
                    return flag, None, None  
            else:
                # then this sequence is random N6 RNA
                flag = 2
                shift3, bc3 = assign_barcode(shifted_seq, round=1, mismatch=mismatch)
                shift2, bc2 = assign_barcode(shifted_seq, round=2, mismatch=mismatch)
                shift1, bc1 = assign_barcode(shifted_seq, round=3, mismatch=mismatch)
                if bc1 in whitelist_bc1[49:96]:
                    return 0, None, None
                # shift4, lib_bc = assign_barcode(shifted_seq, round=4, mismatch=mismatch)
                umi = 'NNNNNNNNNN'
                header1 = read1.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/1'
                read1.name = header1
                header2 = read2.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/2'
                read2.name = header2    
                if bc1 and bc2 and bc3:
                    return flag, read1, read2
                else:    
                    flag = 0
                    return flag, None, None  
        elif self._mode == "RNA":
            phase = find_phase(read2.sequence)
            if phase == -1:
                flag = 0
                return flag, None, None
            shifted_seq = read2.sequence[phase:] + 'T'*phase
            dT = shifted_seq[dt_start:dt_start+T_num]
            if Levenshtein.distance(dT, match) <= 2:
                # then this sequence is RNA
                flag = 1
                shift3, bc3 = assign_barcode(shifted_seq, round=1, mismatch=mismatch)
                shift2, bc2 = assign_barcode(shifted_seq, round=2, mismatch=mismatch)
                shift1, bc1 = assign_barcode(shifted_seq, round=3, mismatch=mismatch)
                shift4, lib_bc = assign_barcode(shifted_seq, round=4, mismatch=mismatch)
                umi = shifted_seq[111+shift4:121+shift4]
                header1 = read1.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + '-' + lib_bc + ':' + umi + '/1'
                read1.name = header1
                header2 = read2.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + '-' + lib_bc + ':' + umi + '/2'
                read2.name = header2    
                if bc1 and bc2 and bc3 and lib_bc:
                    return flag, read1, read2
                else:    
                    flag = 0
                    return flag, None, None  
            else:
                flag = 0
                return flag, None, None  
        elif self._mode == "RNA_mixcr":
            phase = find_phase(read2.sequence)
            if phase == -1:
                flag = 0
                return flag, None, None
            shifted_seq = read2.sequence[phase:] + 'T'*phase
            shifted_qual = read2.qualities[phase:] + '?'*phase
            dT = shifted_seq[dt_start:dt_start+T_num]
            if Levenshtein.distance(dT, match) <= 2:
                # then this sequence is RNA
                flag = 1
                shift3, bc3 = assign_barcode(shifted_seq, round=1, mismatch=mismatch)
                shift2, bc2 = assign_barcode(shifted_seq, round=2, mismatch=mismatch)
                shift1, bc1 = assign_barcode(shifted_seq, round=3, mismatch=mismatch)
                shift4, lib_bc = assign_barcode(shifted_seq, round=4, mismatch=mismatch)
                umi = shifted_seq[111+shift4:121+shift4]
                header1 = read1.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + '-' + lib_bc + ':' + umi + '/1'
                read1.name = header1
                header2 = read2.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + '-' + lib_bc + ':' + umi + '/2'
                read2.name = header2    
                if bc1 and bc2 and bc3 and lib_bc:
                    seq_R2 = bc1 + bc2 + bc3 + lib_bc + umi
                    # dT [143,150] | sample barcode [122,124] | UMI [112,121] | barcode 1 [74,81] | barcode 2 [36,43] | barcode 3 [1,5]
                    qual_R2 = shifted_qual[73+shift1:81+shift1] + shifted_qual[35+shift2:43+shift2] + shifted_qual[0+shift3:5+shift3] + shifted_qual[121+shift4:+124+shift4] + shifted_qual[111+shift4:121+shift4]
                    read2.sequence = seq_R2
                    read2.qualities = qual_R2
                    return flag, read1, read2
                else:    
                    flag = 0
                    return flag, None, None  
            else:
                flag = 0
                return flag, None, None  
        elif self._mode in ["BCR", "TCR"]:
            flag = 2
            phase = find_phase_immune_repertoire(read2.sequence)
            if phase == -1:
                flag = 0
                return flag, None, None
            shifted_seq = read2.sequence[phase:] + 'N'*phase
            shifted_qual = read2.qualities[phase:] + '?'*phase
            # shift0, lib_bc = assign_barcode_immune_repertoire(shifted_seq, round=1, mismatch=mismatch)
            shift1, bc1 = assign_barcode_immune_repertoire(shifted_seq, round=2, mismatch=mismatch)
            shift2, bc2 = assign_barcode_immune_repertoire(shifted_seq, round=3, mismatch=mismatch)
            shift3, bc3 = assign_barcode_immune_repertoire(shifted_seq, round=4, mismatch=mismatch)
            umi = reverse_complement(shifted_seq[21:31], tr)
            header1 = read1.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/1'
            read1.name = header1
            header2 = read2.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/2'
            read2.name = header2   
            if bc1 and bc2 and bc3:
                seq_R2 = bc1 + bc2 + bc3 + umi
                # sample barcode: [19,21] | UMI [22,31] | barcode1 [62,69] | barcode2 [100,107] | barcode3 [138,142]
                qual_R2 = shifted_qual[61+shift1:69+shift1] + shifted_qual[99+shift2:107+shift2] + shifted_qual[137+shift3:142+shift3] + shifted_qual[21:31]
                read2.sequence = seq_R2
                read2.qualities = qual_R2
                return flag, read1, read2
            else:    
                flag = 0
                return flag, None, None
        elif self._mode in ["BCR_v2", "TCR_v2"]:
            flag = 2
            phase = find_phase_immune_repertoire(read2.sequence)
            if phase == -1:
                flag = 0
                return flag, None, None
            shifted_seq = read2.sequence[phase:] + 'N'*phase
            shifted_qual = read2.qualities[phase:] + '?'*phase
            shift1, bc1 = assign_barcode_immune_repertoire_v2(shifted_seq, round=1, mismatch=mismatch)
            shift2, bc2 = assign_barcode_immune_repertoire_v2(shifted_seq, round=2, mismatch=mismatch)
            shift3, bc3 = assign_barcode_immune_repertoire_v2(shifted_seq, round=3, mismatch=mismatch)
            umi = reverse_complement(shifted_seq[18:26], tr)
            header1 = read1.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/1'
            read1.name = header1
            header2 = read2.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/2'
            read2.name = header2   
            if bc1 and bc2 and bc3:
                seq_R2 = bc1 + bc2 + bc3 + umi
                # CCATCTGAGCCACCAGCGNNNNNNNNACAGCAGTCGCTTCATCGGACGATCATGGGNNNNNNNNCAAGTATGCAGCGCGCTCAAGCACGTGGATNNNNNNNNAGTCGTACGCCGATGCGAAACATCGGCCACNNNNNAAGTCGGAGGCCAAGCGGAAGCAGTGGTATCAACGCAGAGT
                # UMI [19,26] | barcode1 [59,66] | barcode2 [97,104] | barcode3 [135,139]
                qual_R2 = shifted_qual[56+shift1:64+shift1] + shifted_qual[94+shift2:102+shift2] + shifted_qual[132+shift3:137+shift3] + shifted_qual[18:26]
                read2.sequence = seq_R2
                read2.qualities = qual_R2
                return flag, read1, read2
            else:    
                flag = 0
                return flag, None, None
        elif self._mode in ["BCR_v2_filt_untag", "TCR_v2_filt_untag"]:
            flag = 2
            phase = find_phase_immune_repertoire(read2.sequence)
            if phase == -1:
                flag = 0
                return flag, None, None
            shifted_seq = read2.sequence[phase:] + 'N'*phase
            shifted_qual = read2.qualities[phase:] + '?'*phase
            shift1, bc1 = assign_barcode_immune_repertoire_v2(shifted_seq, round=1, mismatch=mismatch)
            shift2, bc2 = assign_barcode_immune_repertoire_v2(shifted_seq, round=2, mismatch=mismatch)
            shift3, bc3 = assign_barcode_immune_repertoire_v2(shifted_seq, round=3, mismatch=mismatch)
            umi = reverse_complement(shifted_seq[18:26], tr)
            header1 = read1.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/1'
            read1.name = header1
            header2 = read2.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/2'
            read2.name = header2   
            # ACGAGACCGGAAAGATGTGTATAAGAGACAG AGTCTCTCAGCTGGTACACG
            if bc1 and bc2 and bc3:
                seq_R2 = bc1 + bc2 + bc3 + umi
                # CCATCTGAGCCACCAGCGTTACTCGTACAGCAGTCGCTTCATCGGACGATCATGGG
                # CCATCTGAGCCACCAGCGNNNNNNNNACAGCAGTCGCTTCATCGGACGATCATGGGNNNNNNNNCAAGTATGCAGCGCGCTCAAGCACGTGGATNNNNNNNNAGTCGTACGCCGATGCGAAACATCGGCCACNNNNNAAGTCGGAGGCCAAGCGGAAGCAGTGGTATCAACGCAGAGT
                # UMI [19,26] | barcode1 [59,66] | barcode2 [97,104] | barcode3 [135,139]
                qual_R2 = shifted_qual[56+shift1:64+shift1] + shifted_qual[94+shift2:102+shift2] + shifted_qual[132+shift3:137+shift3] + shifted_qual[18:26]
                # ALl bases of barcodes and umis must have Q>30 (>'?')
                if not all(map(lambda x : ord(x) > 63, list(qual_R2))):
                    return 0, None, None
                read2.sequence = seq_R2
                read2.qualities = qual_R2
                return flag, read1, read2
            else:    
                flag = 0
                return flag, None, None
        elif self._mode == "Vmix":
            phase = find_phase(read2.sequence)
            if phase == -1:
                flag = 0
                return flag, None, None
            shifted_seq = read2.sequence[phase:] + 'T'*phase
            shifted_qual = read2.qualities[phase:] + '?'*phase
            # then this sequence is RNA
            flag = 1
            shift3, bc3 = assign_barcode(shifted_seq, round=1, mismatch=mismatch)
            shift2, bc2 = assign_barcode(shifted_seq, round=2, mismatch=mismatch)
            shift1, bc1 = assign_barcode(shifted_seq, round=3, mismatch=mismatch)
            # shift4, lib_bc = assign_barcode(shifted_seq, round=4)
            umi = shifted_seq[111+shift3:121+shift3]
            header1 = read1.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/1'
            read1.name = header1
            header2 = read2.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/2'
            read2.name = header2    
            if bc1 and bc2 and bc3:
                seq_R2 = bc1 + bc2 + bc3 + umi
                # dT [143,150] | sample barcode [122,124] | UMI [112,121] | barcode 1 [74,81] | barcode 2 [36,43] | barcode 3 [1,5]
                qual_R2 = shifted_qual[73+shift1:81+shift1] + shifted_qual[35+shift2:43+shift2] + shifted_qual[0+shift3:5+shift3] + shifted_qual[111+shift3:121+shift3]
                read2.sequence = seq_R2
                read2.qualities = qual_R2                
                return flag, read1, read2
            else:    
                flag = 0
                return flag, None, None 
        elif self._mode == "Vmix_v2":
            phase = find_phase(read2.sequence)
            if phase == -1:
                flag = 0
                return flag, None, None
            shifted_seq = read2.sequence[phase:] + 'T'*phase
            shifted_qual = read2.qualities[phase:] + '?'*phase
            # then this sequence is RNA
            flag = 1
            shift3, bc3 = assign_barcode(shifted_seq, round=1, mismatch=mismatch)
            shift2, bc2 = assign_barcode(shifted_seq, round=2, mismatch=mismatch)
            shift1, bc1 = assign_barcode(shifted_seq, round=3, mismatch=mismatch)
            # shift4, lib_bc = assign_barcode(shifted_seq, round=4)
            umi = shifted_seq[111+shift3:119+shift3]
            header1 = read1.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/1'
            read1.name = header1
            header2 = read2.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/2'
            read2.name = header2    
            if bc1 and bc2 and bc3:
                seq_R2 = bc1 + bc2 + bc3 + umi
                # dT [143,150] | sample barcode [122,124] | UMI [112,121] | barcode 1 [74,81] | barcode 2 [36,43] | barcode 3 [1,5]
                qual_R2 = shifted_qual[73+shift1:81+shift1] + shifted_qual[35+shift2:43+shift2] + shifted_qual[0+shift3:5+shift3] + shifted_qual[111+shift3:119+shift3]
                read2.sequence = seq_R2
                read2.qualities = qual_R2                
                return flag, read1, read2
            else:    
                flag = 0
                return flag, None, None     
        elif self._mode == "ATAC":
            flag = 3
            phase = find_phase(read2.sequence)
            if phase == -1:
                flag = 0
                return flag, None, None
            shifted_seq = read2.sequence[phase:] + 'N'*phase
            shift3, bc3 = assign_barcode(shifted_seq, round=1, mismatch=mismatch)
            shift2, bc2 = assign_barcode(shifted_seq, round=2, mismatch=mismatch)
            shift1, bc1 = assign_barcode(shifted_seq, round=3, mismatch=mismatch)
            header1 = read1.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + '/1'
            read1.name = header1
            header2 = read2.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + '/2'
            read2.name = header2    
            if bc1 and bc2 and bc3:
                return flag, read1, read2
            else:    
                flag = 0
                return flag, None, None  
        elif self._mode == "RNA_old":
            dT = read2.sequence[dt_start:dt_start+T_num]
            if Levenshtein.distance(dT, match) <= 2:
                # then this sequence is RNA, with no sample barcode
                flag = 1
                shift3, bc3 = assign_barcode(read2.sequence, round=1, mismatch=mismatch)
                shift2, bc2 = assign_barcode(read2.sequence, round=2, mismatch=mismatch)
                shift1, bc1 = assign_barcode(read2.sequence, round=3, mismatch=mismatch)
                umi = read2.sequence[111+shift1:121+shift1]
                header1 = read1.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/1'
                read1.name = header1
                header2 = read2.name[:-2] + ':' + bc1 + '-' + bc2 + '-' + bc3 + ':' + umi + '/2'
                read2.name = header2    
                if bc1 and bc2 and bc3:
                    return flag, read1, read2
                else:    
                    flag = 0
                    return flag, None, None  
            else:
                flag = 0
                return flag, None, None  
        else:
            raise ValueError

    def run(self):
        try:
            stats = 0
            while True:
                total = 0
                c1 = 0
                c2 = 0
                # Notify the reader that we need data
                self._q.put(self._id)
                chunk_index = self._conn_in.recv()
                if chunk_index == -1:
                    # reader is done
                    break
                infiles = self._make_input_files()
                self._reader = infiles.open()
                # We turn original output files into a temporary buffer (io.Bytesio). 
                # We open this temporary buffer with dnaio.open and write dnaio.Sequence object to it.
                outfiles = self._original_outfiles.as_bytesio()
                o = dnaio.open(file1=outfiles.out1, file2=outfiles.out2, mode='w', qualities=True)
                # read1 and read2 are dnaio._core.Sequence objects
                for read1, read2 in self._reader:     
                    flag, read1, read2 = self.process_read(read1, read2)
                    total += 1
                    c1_add, c2_add = self.valid_read_counts(flag)
                    c1 += c1_add
                    c2 += c2_add
                    self._deliver_processed_reads(o, read1, read2)
                # Send the processed chunk to writer. The writer will receive chunk index + data(bytes), or chunk index + b''
                self._send_outfiles(outfiles, chunk_index, total, c1, c2)  
                self._close(o, infiles, outfiles) 
            self._conn_out.send(-1)
            self._conn_out.send(stats)
        except Exception as e:
            self._conn_out.send(-2)
            self._conn_out.send((e, traceback.format_exc()))
        return None

    def _make_input_files(self) -> InputFiles:
        data = self._conn_in.recv_bytes()
        input1 = io.BytesIO(data)
        data = self._conn_in.recv_bytes()
        input2 = io.BytesIO(data)
        return InputFiles(input1, input2)

    def _send_outfiles(self, outfiles: OutputFiles, chunk_index: int, n_reads: int, c1: int, c2: int):
        self._conn_out.send(chunk_index)
        self._conn_out.send(n_reads)
        self._conn_out.send(c1)
        self._conn_out.send(c2)
        for f in outfiles:
            f.flush()
            assert isinstance(f, io.BytesIO)
            processed_chunk = f.getvalue() # if empty, getvalue() returns b''
            self._conn_out.send_bytes(processed_chunk)
    
    def _deliver_processed_reads(self, o, read1, read2):
        if read1 is not None and read2 is not None:
            try:
                o.write(read1, read2)
            except TypeError:
                pass
        return None

    def valid_read_counts(self, flag: int = 0) -> Tuple:
        c1, c2 = 0, 0
        if flag == 3:
            c1 += 1
        else:
            if flag == 1:
                c1 += 1
            elif flag == 2:
                c2 += 1
        return c1, c2

    def _close(self, o, infiles, outfiles) -> None:
        o.close()
        self._close_input(infiles)
        self._close_output(outfiles)
        return None

    def _close_input(self, infiles) -> None:
        # self._reader: dnaio.TwoFilePairedEndReader
        self._reader.close()
        if infiles is not None:
            infiles.close() 
        return None   

    def _close_output(self, outfiles) -> None:
        for f in outfiles:
            assert isinstance(f, io.BytesIO)
            f.close()
        return None
 

class OrderedChunkWriter:
    """
    We may receive chunks of processed data from worker processes
    in any order. This class writes them to an output file in
    the correct order.
    """
    def __init__(self, outfile):
        self._chunks = dict()
        self._current_index = 0
        self._outfile = outfile

    def write(self, data, index):
        self._chunks[index] = data
        while self._current_index in self._chunks:
            self._outfile.write(self._chunks[self._current_index])
            del self._chunks[self._current_index]
            self._current_index += 1
        return None

def try_receive(connection):
    """
    Try to receive data over self.connection and return it.
    If an exception was received, raise it.
    """
    result = connection.recv()
    if result == -2:
        # An exception has occurred on the other end
        e, tb_str = connection.recv()
        raise e
    return result


if __name__ == "__main__":
    # Open files
    start_time = time.time()
    r1, r2 = open_input_files(read1, read2)
    outfiles = open_output_files(mode, out1, out2)
    # DataProcesser put its id in this queue to notify FileReader that it need data to process
    working_queue = Queue()
    # Total Process list
    process_list = list()
    # pipe_list_in is used to send data from FileReader to DataProcesser
    # pipe_list_out is used to send processed data from DataProcesser to main process (DataWriter)
    pipe_list_in, pipe_list_out = list(), list()
    for i in range(cores):
        pipe_list_in.append(Pipe())
    for i in range(cores):
        pipe_list_out.append(Pipe())
        pipe_out_recv = list([pipe_list_out[i][1] for i in range(len(pipe_list_out))])
    
    reader = FileReader([r1, r2], working_queue, [pipe_list_in[i][0] for i in range(len(pipe_list_in))])
    reader.start()
    for i in range(cores):
        process_list.append(DataProcesser(i, working_queue, mode, pipe_list_in[i][1], pipe_list_out[i][0], outfiles))
    # Start process
    for p in process_list:
        p.start()
    # Create writer
    writers = list()
    for f in outfiles:
        writers.append(OrderedChunkWriter(f))
    # Fetch data and write files
    total = 0
    c1 = 0
    c2 = 0
    while pipe_out_recv:
        for conn in connection.wait(pipe_out_recv):
            chunk_index = try_receive(conn)
            if chunk_index == -1:
                # the worker is done
                cur_stats = try_receive(conn)
                pipe_out_recv.remove(conn)
                continue
            total += try_receive(conn)
            c1 += try_receive(conn)
            c2 += try_receive(conn)
            for writer in writers:
                data = conn.recv_bytes()
                # assert isinstance(writer, OrderedChunkWriter)
                writer.write(data, chunk_index)

    for p in process_list:
        p.join()
    reader.join()

    end_time = time.time()
    elapsed_time = end_time - start_time
    calculate_statistics(c1, c2, total, elapsed_time, mode)
    r1.close()
    r2.close()
    for f in outfiles:
        f.close()


