import argparse
import math as m
import feature as f
import read_file as r
import gen_tools as t

class Universal_error(Exception):
    pass

parser = argparse.ArgumentParser(description=
    'play around to make estimation on required computational power')
parser.add_argument('--fasta', required=True, type=str,
    metavar='<path>', help='path to fasta file')
parser.add_argument('--dict', required=False, type=str,default='Jeanne.js',
    metavar='<path>', help='path to dict(.js) file')
arg = parser.parse_args()

def compute_num_of_shared_gapped_k_mer(seq1, seq2, k):
    if len(seq1) != len(seq2):
        raise Universal_error("Check seqs's len to do gapped k-mer")
    mismatch = 0
    length = len(seq1)
    for i in range(length):
        if seq1[i]!=seq2[i]: mismatch+=1
    l_m = length - mismatch
    if l_m < k: return 0
    else:
        return int(m.factorial(l_m)/(m.factorial(k)*m.factorial(l_m-k)))


# Get fasta sequence and look-up dictionary
seq_file = arg.fasta
dict_file = arg.dict
fas_read = r.FASTA(seq_file)
dict = t.decode_json(dict_file)

# Let's play with start site first
init = dict["init"]
for entry in fas_read:
    # Generate sliding-window seqs
    vectors = []
    for i in range(len(entry.seq)-10):
        seq = entry.seq[i:i+10]
        # For 10bp start site, the parameters are set as 4,3,3 now
        site = f.Site(seq, preLen=4, keyLen=3, suffLen=3)
        if site.key in init["start"]:
        # If the key has never seen in training set, then discard it
            vector = [init["start"][site.key]]
            for ele in init["pre"]:
                v = compute_num_of_shared_gapped_k_mer(ele, site.pre, 2)
                vector.append(v)
            for ele in init["aa1"]:
                v = compute_num_of_shared_gapped_k_mer(ele, site.suff, 2)
                vector.append(v)
            vectors.append(vector)
    print(vectors)
