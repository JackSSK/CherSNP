import os
import re
import argparse
from sklearn import svm
# from itertools import chain
import read_file as read
import seq_tools as seqt
# import feature as fea

# X = [[1,1,1], [2,2,2]]
# y = [0,1]
# clf = svm.SVC()
# clf.fit(X,y)
# print(clf.predict([[2,2,1]]))

class ID_error(Exception):
    pass

class Trainer:
    def __init__(self, seq_file, gff_file):
        self.gff = {}
        self._processGFF(gff_file)
        self._processFAS(seq_file)

    def _processGFF(self,gff_file):
        gff_read = read.GFF(gff_file)
        unlink_cds = []
        for srcseq in gff_read:
            for entry in srcseq:
                if entry.seqid not in self.gff:
                    self.gff[entry.seqid] = {}

                if entry.type == 'mRNA':
                    if re.search(r'ID=',entry.attr):
                        id = entry.attr.split('ID=')[1].split(';')[0]
                        if id not in self.gff[entry.seqid]:
                            self.gff[entry.seqid][id] = {
                                'beg':int(entry.beg),
                                'end':int(entry.end),
                                'strand':entry.strand,
                                'cds':[]
                            }
                        else:
                            raise ID_error('mRNA read in twice')
                    else:
                        raise ID_error('ID missing')
                        continue

                elif entry.type == 'CDS':
                    if re.search(r'Parent=',entry.attr):
                        pid = entry.attr.split('Parent=')[1].split(';')[0]
                        if pid not in self.gff[entry.seqid]:
                            unlink_cds.append([entry.seqid, pid, int(entry.beg), int(entry.end)])
                        else:
                            ele = [int(entry.beg), int(entry.end)]
                            if ele not in self.gff[entry.seqid][pid]['cds']:
                                self.gff[entry.seqid][pid]['cds'].append(ele)
                            else:
                                raise ID_error('CDS already exist')
                    else:
                        raise ID_error('Parent ID missing')
                        continue
        gff_read.close()
        for cds in unlink_cds:
            seqid = cds[0]
            pid = cds[1]
            ele = [cds[2], cds[3]]
            try:
                if ele not in self.gff[seqid][pid]['cds']:
                    self.gff[seqid][pid]['cds'].append(ele)
                else:
                    raise ID_error('CDS already exist')
            except:
                raise ID_error('Cannot load unlink CDS')

    def _processFAS(self, seq_file):
        fas_read = read.FASTA(seq_file)
        for entry in fas_read:
            for trans in self.gff[entry.id]:
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
