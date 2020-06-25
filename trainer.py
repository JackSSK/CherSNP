import os
import re
from sklearn import svm
# from itertools import chain
import read_file as read
import seq_tools as seqt
import gffer
# import feature as fea
# Nothing here
# X = [[1,1,1], [2,2,2]]
# y = [0,1]
# clf = svm.SVC()
# clf.fit(X,y)
# print(clf.predict([[2,2,1]]))

class ID_error(Exception):
    pass

class Trainer:
    def __init__(self, seq_file, gff_file):
        self.gff = gffer.Process(gff_file).gff
        self._processFAS(seq_file)

    def _processFAS(self, seq_file):
        fas_read = read.FASTA(seq_file)
        for entry in fas_read:
            for trans in self.gff[entry.id]:
                if self.gff[entry.id][trans]['type'] == 'transcript':
                    beg = self.gff[entry.id][trans]['beg']
                    end = self.gff[entry.id][trans]['end']
                    strand = self.gff[entry.id][trans]['strand']
                    seq = entry.seq[beg-1:end]
                    coord = []
                    cds_contest = []
                    init = []
                    ter = []
                    ent_cds = []
                    ext_cds = []
                    for ele in self.gff[entry.id][trans]['cds']:
                        temp = [int(ele[0])-beg, int(ele[1])-beg]
                        coord.append(temp)
                    if strand == '+':
                        coord = sorted(coord, key=lambda x: x[0])
                        init.append(coord[0][0])
                        ter.append(coord[-1][1])
                        for ele in coord:
                            enter = ele[0]
                            exit = ele[1]
                            if enter != init[0]: ent_cds.append(enter)
                            if exit != ter[0]: ext_cds.append(exit)
                            inseq = seqt.complementary(seq[enter-4:enter+6])
                            outseq = seqt.complementary(seq[exit-5:exit+5])
                            cds_contest.append([inseq, outseq])
                    elif strand == '-':
                        coord = sorted(coord, key=lambda x: -x[0])
                        init.append(coord[0][1])
                        ter.append(coord[-1][0])
                        for ele in coord:
                            enter = ele[1]
                            exit = ele[0]
                            if enter != init[0]: ent_cds.append(enter)
                            if exit != ter[0]: ext_cds.append(exit)
                            outseq = seqt.complementary(seq[exit-4:exit+6])
                            inseq = seqt.complementary(seq[enter-5:enter+5])
                            cds_contest.append([inseq, outseq])

                    print(coord, cds_contest, init, ent_cds, ext_cds, ter)
        fas_read.close()
