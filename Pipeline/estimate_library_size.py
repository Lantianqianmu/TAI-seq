import pandas as pd
import argparse
import time
import numpy as np
from scipy.optimize import fsolve

# python /data/xrz/JoINT-seq/Pipeline/estimate_library_size.py -r /data/xrz/JoINT-seq/Data/plate2-rna/plate2-rna_hg38_1_reads_count.txt -u /data/xrz/JoINT-seq/Data/plate2-rna/plate2-rna_hg38_1_umis_count.txt -m 100 -o /data/xrz/JoINT-seq/Data/plate2-rna/plate2-rna_hg38_1_estimate_libsize.txt

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--reads', dest='READ', type=str, action='store', help='Reads per cell')
parser.add_argument('-u', '--umis', dest='UMI', type=str, action='store', help='Umis per cell')
parser.add_argument('-m', '--min', dest='MIN', type=int, action='store', help='Min umis/fragments per cell')
parser.add_argument('-o', '--out', dest='OUT', default=None, type=str, action='store', help='Output tsv file')
# parser.add_argument('-e', '--estimate', dest='EST', default=0, type=int, action='store', help='Estimated library size. Typically 3000 for RNA and 6000 for ATAC.')
args = parser.parse_args()

reads_file = args.READ
umis_file = args.UMI
min_umi = args.MIN
out = args.OUT
# est = args.EST


def read_file(file):
    data = pd.read_table(file, sep = "\t", header = None, names = ["Barcodes", "Counts"], index_col = None)
    return data

def filter_barcodes(min_umi, df):
    return df.loc[df['Counts'] >= min_umi]

def estimate_library_size(df, estimate=0, paired=False):
    size = []
    for i in range(data.shape[0]):
        reads = data.iloc[i,1]
        umis = data.iloc[i,2]
        # x = Symbol('x')
        # a = solve(x * (1 - exp(-reads/x)) - umis)[0]
        # a = np.float(a)
        # size.append(np.around(a).astype(int))
        a = fsolve(lambda x: x * (1 - np.exp(-reads/x)) - umis, 0)[0]
        est_libsize = np.around(a).astype(int)
        size.append(est_libsize)
    df['Size'] = size
    return df
      
if __name__ == "__main__":
    start_time = time.time()

    df_reads = read_file(file = reads_file)
    df_umis = read_file(file = umis_file)
    df_umis = filter_barcodes(min_umi, df_umis)
    data = pd.merge(df_reads, df_umis, how = "right", on = "Barcodes")
    data.columns = ["Barcodes", "Reads", "Umis"]
    data = estimate_library_size(data)
    data.to_csv(out, sep = "\t", header = True, index = False)

    end_time = time.time()
    elapsed_time = end_time - start_time
    formatted_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    print("Total time elapsed running {}: {}".format(__file__, formatted_time))



