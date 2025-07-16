import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest='INPUT', type=str, action='store', help='input tsv file')
parser.add_argument('-o', '--output', dest='OUTPUT', type=str, action='store', help='output tsv file')

args = parser.parse_args()

infile = args.INPUT
outfile = args.OUTPUT

def remove_multi_assignment(data: pd.DataFrame):
    print("Removing multimapping reads that are assigned to multiple genes...")
    data["mapped_gene_number"] = data.groupby(["read_id"])["gene"].transform("nunique")
    return(data[data["mapped_gene_number"] == 1].drop_duplicates('read_id'))


if __name__ == "__main__":
    data = pd.read_csv(infile, sep = "\t", header = None, names = ["read_id", "barcode", "umi", "chr", "start", "gene"])
    data_out = remove_multi_assignment(data)
    data_out[["read_id", "barcode", "umi", "chr", "start", "gene"]].to_csv(outfile, sep = "\t", header = False, index = False)


    


