import os
import re
import argparse
from sklearn import svm
import read_file as read
# import feature as fea

parser = argparse.ArgumentParser(description='Just for test now')
parser.add_argument('--fasta', required=True, type=str,
    metavar='<path>', help='path to fasta file')
parser.add_argument('--gff', required=True, type=str,
    metavar='<path>', help='path to GFF3 file')
arg = parser.parse_args()
seq_file = arg.fasta
gff_file = arg.gff


# X = [[1,1,1], [2,2,2]]
# y = [0,1]
# clf = svm.SVC()
# clf.fit(X,y)
#
# print(clf.predict([[2,2,1]]))


fas_read = read.FASTA(seq_file)
for entry in fas_read:
    print(entry.id)
    print(entry.seq[0:10])
fas_read.close()
gff_read = read.GFF(gff_file)
for srcseq in gff_read:
    for entry in srcseq:
        if entry.type == 'gene' or entry.type == 'CDS':
            print(entry.seqid)
gff_read.close()
