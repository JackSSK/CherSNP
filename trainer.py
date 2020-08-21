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

class Initiate_Site:
    def __init__(self,seq):
        self.pre = seq[:4]
        self.start = seq[4:7]
        self.first = seq[7:]
    def out(self):
        return [self.pre, self.start, self.first]

class Trainer:
    def __init__(self, seq_file, gff_file):
        # Process in GFF file
        self.gff = gffer.Process(gff_file).gff
        # Get observations and correspond dictionary
        self.dict, self.observ = self._observer(seq_file)
        # Process in observations and get classifiers

    def _observer(self, seq_file):
        #Initiating dicts for stroing correct observations
        dict = {
            "init":{
                "pre":seqt.generate(k=4, type='dict'),
                "start":seqt.generate(k=3, type='dict'),
                "aa1":seqt.generate(k=3, type='dict')
            },
            "term":{
                "last":seqt.generate(k=3, type='dict'),
                "stop":seqt.generate(k=3, type='dict'),
                "suf":seqt.generate(k=4, type='dict'),
            },
            "entCDS":{

            },
            "outCDS":{

            }
        }

        observ = {
            "init":{
                "correct":[],
                "wrong":[]
            },
            "term":{
                "correct":[],
                "wrong":[]
            },
            "entCDS":{
                "correct":[],
                "wrong":[]
            },
            "outCDS":{
                "correct":[],
                "wrong":[]
            }
        }

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
                            inseq = seq[enter-5:enter+7]
                            outseq = seq[exit-6:exit+6]
                            cds_contest.append([inseq, outseq])

                            if enter != init: ent_cds.append(enter)
                            else: self._update_init(dict,observ,inseq)
                            if exit != ter: ext_cds.append(exit)
                            # else: self._update_aaEnd(outseq)

                    elif strand == '-':
                        coord = sorted(coord, key=lambda x: -x[0])
                        init = coord[0][1]
                        ter = coord[-1][0]
                        for ele in coord:
                            enter = ele[1]
                            exit = ele[0]
                            outseq = seqt.complementary(seq[exit-5:exit+7])
                            inseq = seqt.complementary(seq[enter-6:enter+6])
                            cds_contest.append([inseq, outseq])

                            if enter != init: ent_cds.append(enter)
                            else: self._update_init(dict,observ,inseq)
                            if exit != ter: ext_cds.append(exit)
                            # else: self._update_aaEnd(outseq)

                    # print(trans,'\n', coord, cds_contest, init, ent_cds, ext_cds, ter)
        fas_read.close()
        return dict, observ

    # This is to store observations of how transcript initiated
    # Stuffs and Start Codon
    # Now wrong cases are gained by shift 1bp,
    # Maybe shift 3bp instead of 1?
    def _update_init(self, dict, observ, contest, prefix=True):
        if prefix:
            wrong1 = Initiate_Site(contest[:-2])
            wrong2 = Initiate_Site(contest[2:])
            correct = Initiate_Site(contest[1:-1])
            observ["init"]["correct"].append(correct.out())
            observ["init"]["wrong"].append(wrong1.out())
            observ["init"]["wrong"].append(wrong2.out())
            try:
                dict["init"]["start"][correct.start]+=1
                dict["init"]["aa1"][correct.first]+=1
                dict["init"]["pre"][correct.pre]+=1
            except Exception as e:
                raise Read_error('What the Hell is ' +
                    correct.start + ' or ' + correct.first)

    # This is to store observations of how transcript Ended
    # Stuffs and End Codon
    # Need to change to sth like _update_init
    # def _update_aaEnd(self, contest, sufix=True):
    #     wrong1 = contest[:-2]
    #     wrong2 = contest[2:]
    #     correct = contest[1:-1]
    #
    #     end_codon = contest[4:7]
    #     last = contest[1:4]
    #     try:
    #         self.aaEnd[end_codon]+=1
    #         self.aaL[last]+=1
    #     except Exception as e:
    #         raise Read_error('What the Hell is ' +
    #             end_codon + ' or ' + last)
