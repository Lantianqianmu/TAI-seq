import matplotlib.pyplot as plt
import argparse
import numpy as np
import math
import re

# python /home/zeemeeuw/YangLab/JoINT-seq/Pipeline/TCR/plot_distribution.py -i /home/zeemeeuw/YangLab/JoINT-seq/Data/JH8-rna/JH8-rna_mm39_RNA.count.txt -o /home/zeemeeuw/YangLab/JoINT-seq/Data/JH8-rna/JH8-rna_mm39_RNA.count.png -s 1 -e 3000 -f 1
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest='IN', type=str, action='store', help='Count distribution')
parser.add_argument('-o', '--output', dest='OUT', type=str, action='store', help='Output plot')
parser.add_argument('-s', '--start', dest='START', type=int, action='store', help='Barcode rank of the first point (min=1)')
parser.add_argument('-e', '--end', dest='END', type=int, action='store', help='Barcode rank of the last point')
parser.add_argument('--logx', dest='LOGX', action='store_true', help='Using log10 x coordinate')
parser.add_argument('--logy', dest='LOGY', action='store_true', help='Using log10 y coordinate')
parser.add_argument('-f', '--field', dest='FIELD', type=int, action='store', default = 1, help='Column index to use (0-based)')
parser.add_argument('-t', '--header', dest='HEADER', action='store_true', help='The file contains header')
parser.add_argument('--otsu', dest='OTSU', action='store_true', help='The file contains header')

args = parser.parse_args()

input = args.IN
output = args.OUT
start = args.START-1
end = args.END
use_logx = args.LOGX
use_logy = args.LOGY
field = args.FIELD
header = args.HEADER
plot_otsu = args.OTSU
cell_threshold = 0

def plot_sorted_list(data, threshold, file):
    if use_logy:  
        data_y = np.log10(data)       
    else:
        data_y = data
    if len(data) >= end:
        data_x = list(range(start, end))
    else:
        data_x = list(range(start, len(data)))
    if use_logx:  
        data_x = [i+1 for i in data_x]
        data_x = np.log10(data_x)
    else:
        data_x = data_x

    plt.scatter(data_x, data_y, color = '#055189', s = 1)
    if use_logx:  
        plt.xlabel("Log10 barcode rank")
        plt.xlim(0, 6)
    else:
        plt.xlabel("Barcode rank")
        plt.xlim(0, end)
    if use_logy:  
        plt.ylabel("Log10 read counts")
        if threshold > 0:
            plt.axhline(np.log10(threshold), color='green', alpha = 0.5, linestyle = "dashed")
    else:
        plt.ylabel("Read counts")
        if threshold > 0:
            plt.axhline(threshold, color='green', alpha = 0.5, linestyle = "dashed")       
    # plt.show()
    plt.savefig(file)

def otsu(data, width = 0.2, logbase = 10):
    qc_flag = qc_otsu(data)
    if not qc_flag:
        return 0
    array = np.log(data)/np.log(logbase)
    counts, bins = np.histogram(array, bins = np.arange(0, max(array) + width, width))
    bin_centers = bins[:-1]
    # class probabilities for all possible thresholds
    weight1 = np.cumsum(counts)
    weight2 = np.cumsum(counts[::-1])[::-1]
    # class means for all possible thresholds
    mean1 = np.cumsum(counts * bin_centers) / weight1
    mean2 = (np.cumsum((counts * bin_centers)[::-1]) / weight2[::-1])[::-1]
    # Clip ends to align class 1 and class 2 variables:
    # The last value of ``weight1``/``mean1`` should pair with zero values in
    # ``weight2``/``mean2``, which do not exist.
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2
    if len(variance12) == 0:
        return 0
    idx = np.nanargmax(variance12)
    threshold = bin_centers[idx]    
    return(math.ceil(logbase ** threshold))

def qc_otsu(data):
    if(len(data)) < 20:
        print("data length < 20. Skip otsu.")
        return False
    if not all(np.array(data) > 0):
        print("All data must be > 0. Skip otsu.")
        return False
    return True

if __name__ == "__main__":
    counts = []
    n = 1
    with open(input, "r") as f:
        for line in f.readlines():
            if header and n == 1:
                n += 1
                continue
            line = line.rstrip()
            cols = line.split("\t")
            counts.append(cols[field])
    data = list(map(int, counts))
    data.sort(reverse = True)
    threshold = otsu(data)
    print("Threshold determined: " + str(threshold))
    cell_threshold = len([x for x in data if x > threshold])
    print("Cell number determined: " + str(cell_threshold))
    if not plot_otsu:
        threshold = 0
    data = data[start:end]
    
    plot_sorted_list(data, threshold, output)
    
    







