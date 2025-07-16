import Levenshtein
import argparse
import xopen
import time

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest='INPUT', type=str, action='store', help='input tsv file')
parser.add_argument('-o', '--output', dest='OUTPUT', type=str, action='store', help='output tsv file')
parser.add_argument('-g', '--gene', dest='GENE', action='store_true', help='Group umi at gene level. Omit alignment position.')
parser.add_argument('-m', '--max_dist', dest='MAXDIST', default=1, type=int, action='store', help='maximum distance of umi to be considered as the same')


args = parser.parse_args()

infile = args.INPUT
outfile = args.OUTPUT
max_dist = args.MAXDIST
gene_level = args.GENE


def cluster_levenshtein(words: list, max_dist: int = 1) -> list:
    groups = []
    clustered = False
    for word_tocheck in words:
        for group in groups:
            for word in group:
                if Levenshtein.distance(word, word_tocheck) <= max_dist:
                    group.append(word_tocheck)
                    clustered = True
                    break
            if clustered:
                break
        if clustered:
            # this word_tocheck has been assigned a cluster
            clustered = False
        else: # new cluster starting point
            groups.append([word_tocheck])    
    return(groups)

def dominant_word(words: list) -> str:
    word_set = set(words)
    max_count, max_word = 0, ""
    for word in word_set:
        current_count = words.count(word)
        if current_count > max_count:
            max_count, max_word = current_count, word
    return(max_word)



if __name__ == "__main__":
    start_time = time.time()
    print("Grouping umis...")

    fin = xopen.xopen(infile, mode = "rt")
    fout = xopen.xopen(outfile, mode = "wt")
    for line in fin:
        line.strip("\t")

    cell_dict = {}
    fin = xopen.xopen(infile, mode = "rt")
    fout = xopen.xopen(outfile, mode = "wt")
    _, barcode_0, umi_0, chr_0, start_0, gene_0 = fin.readline().strip().split("\t")
    cell_dict[barcode_0] = [umi_0]
    if gene_level:
        print("Assign umi by barcode and umi within each gene...")
        for line in fin:
            _, barcode_1, umi_1, chr_1, start_1, gene_1 = line.strip().split("\t")
            if gene_1 == gene_0:
                # same gene, different reads
                if barcode_1 in cell_dict.keys():
                    cell_dict[barcode_1].append(umi_1)
                else:
                    cell_dict[barcode_1] = [umi_1]   
            else:
                for barcode, umi_list in cell_dict.items():
                    if len(umi_list) == 1:
                        fout.write(gene_0 + '\t' + barcode + '\t' + umi_list[0] + '\t' + str(1) + '\n')
                        pass
                    else:
                        clustered_umi = cluster_levenshtein(umi_list, max_dist)
                        for umi_cluster in clustered_umi:
                            final_umi = dominant_word(umi_cluster)
                            fout.write(gene_0 + '\t' + barcode + '\t' + final_umi + '\t' + str(len(umi_cluster)) + '\n')
                barcode_0, umi_0, chr_0, start_0, gene_0 = barcode_1, umi_1, chr_1, start_1, gene_1
                cell_dict = {}
                cell_dict[barcode_0] = [umi_0]   
    else:
        print("Assign umi by barcode, umi and alignment position...")
        for line in fin:
            _, barcode_1, umi_1, chr_1, start_1, gene_1 = line.strip().split("\t")
            if chr_1 == chr_0 and start_1 == start_0 and gene_1 == gene_0:
                # same position, different reads
                if barcode_1 in cell_dict.keys():
                    cell_dict[barcode_1].append(umi_1)
                else:
                    cell_dict[barcode_1] = [umi_1]   
            else:
                for barcode, umi_list in cell_dict.items():
                    if len(umi_list) == 1:
                        fout.write(gene_0 + '\t' + barcode + '\t' + umi_list[0] + '\t' + str(1) + '\n')
                        pass
                    else:
                        clustered_umi = cluster_levenshtein(umi_list, max_dist)
                        for umi_cluster in clustered_umi:
                            final_umi = dominant_word(umi_cluster)
                            fout.write(gene_0 + '\t' + barcode + '\t' + final_umi + '\t' + str(len(umi_cluster)) + '\n')
                barcode_0, umi_0, chr_0, start_0, gene_0 = barcode_1, umi_1, chr_1, start_1, gene_1
                cell_dict = {}
                cell_dict[barcode_0] = [umi_0]
       
    # After looping, process the final umi group or last line
    for barcode, umi_list in cell_dict.items():
        if len(umi_list) == 1:
            fout.write(gene_0 + '\t' + barcode + '\t' + umi_list[0] + '\t' + str(1) + '\n')
            pass
        else:
            clustered_umi = cluster_levenshtein(umi_list, max_dist)
            for umi_cluster in clustered_umi:
                final_umi = dominant_word(umi_cluster)
                fout.write(gene_0 + '\t' + barcode + '\t' + final_umi + '\t' + str(len(umi_cluster)) + '\n')

    fin.close()
    fout.close()
    end_time = time.time()
    formatted_time = time.strftime("%H:%M:%S", time.gmtime(end_time - start_time))
    print("Time used to group umis: {}".format(formatted_time))

