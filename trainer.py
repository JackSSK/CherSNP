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
class Read_error(Exception):
    pass

class ID_error(Exception):
    pass

class Trainer:
    def __init__(self, seq_file, gff_file):
        #Initiating dicts for stroing codon observations
        self.aaSta = seqt.generate(k=3, type='dict')
        self.aaEnd = seqt.generate(k=3, type='dict')
        self.aa1 = seqt.generate(k=3, type='dict')
        self.aaL = seqt.generate(k=3, type='dict')
        #Process in GFF file
        self.gff = gffer.Process(gff_file).gff
        self._processWiFAS(seq_file)



    def _processWiFAS(self, seq_file):
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
                    init = 0
                    ter = 0
                    ent_cds = []
                    ext_cds = []
                    for ele in self.gff[entry.id][trans]['cds']:
                        temp = [int(ele[0])-beg, int(ele[1])-beg]
                        coord.append(temp)

                    if strand == '+':
                        coord = sorted(coord, key=lambda x: x[0])
                        init = coord[0][0]
                        ter = coord[-1][1]
                        for ele in coord:
                            enter = ele[0]
                            exit = ele[1]

                            inseq = seq[enter-4:enter+6]
                            outseq = seq[exit-5:exit+5]
                            cds_contest.append([inseq, outseq])

                            if enter != init: ent_cds.append(enter)
                            else: self._update_aaSta(inseq)
                            if exit != ter: ext_cds.append(exit)
                            else: self._update_aaEnd(outseq)

                    elif strand == '-':
                        coord = sorted(coord, key=lambda x: -x[0])
                        init = coord[0][1]
                        ter = coord[-1][0]
                        for ele in coord:
                            enter = ele[1]
                            exit = ele[0]

                            outseq = seqt.complementary(seq[exit-4:exit+6])
                            inseq = seqt.complementary(seq[enter-5:enter+5])
                            cds_contest.append([inseq, outseq])

                            if enter != init: ent_cds.append(enter)
                            else: self._update_aaSta(inseq)
                            if exit != ter: ext_cds.append(exit)
                            else: self._update_aaEnd(outseq)

                    print(trans,'\n', coord, cds_contest, init, ent_cds, ext_cds, ter)
        fas_read.close()

    # This is to store observations of how transcript initiated
    # Stuffs and Start Codon
    def _update_aaSta(self, contest, prefix=True):
        if prefix:
            start_codon = contest[4:7]
            first = contest[7:]
        try:
            self.aaSta[start_codon]+=1
            self.aa1[first]+=1
        except Exception as e:
            raise Read_error('What the Hell is ' +
                start_codon + ' or ' + first)

    # This is to store observations of how transcript Ended
    # Stuffs and End Codon
    def _update_aaEnd(self, contest, sufix=True):
        end_codon = contest[3:6]
        last = contest[:3]
        try:
            self.aaEnd[end_codon]+=1
            self.aaL[last]+=1
        except Exception as e:
            raise Read_error('What the Hell is ' +
                end_codon + ' or ' + last)
