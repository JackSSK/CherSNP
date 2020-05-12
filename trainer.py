import os
import re
import argparse
from sklearn import svm
import read_file as read
import seq_tools as seqt
# import feature as fea

class ID_error(Exception):
    pass

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
                for ele in self.gff[entry.id][trans]['cds']:
                    temp = [int(ele[0])-beg, int(ele[1])-beg]
                    coord.append(temp)
                if strand == '+':
                    coord = sorted(coord, key=lambda x: x[0])
                elif strand == '-':
                    coord = sorted(coord, key=lambda x: -x[0])
                    for ele in coord:
                        outseq = seqt.complementary(seq[ele[0]-4:ele[0]+6])
                        inseq = seqt.complementary(seq[ele[1]-5:ele[1]+5])
                        cds_contest.append([inseq, outseq])
                print(cds_contest)
        fas_read.close()


Trainer(seq_file, gff_file)
