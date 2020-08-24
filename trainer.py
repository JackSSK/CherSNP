import os
import re
from sklearn import svm
# from itertools import chain
import read_file as read
import seq_tools as seqt
import gffer
import feature
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
        # Process in GFF file
        self.gff = gffer.Process(gff_file).gff
        # Get observations and correspond dictionary
        self.dict, self.observ = self._observer(seq_file)
        # seqt.encode_json([self.dict, self.observ])
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
                "end2":seqt.generate(k=2, type='dict'),
                "first1":seqt.generate(k=1, type='dict'),
            },
            "outCDS":{
                "pre2":seqt.generate(k=2, type='dict'),
                "first2":seqt.generate(k=2, type='dict'),
                "next4":seqt.generate(k=4, type='dict'),
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
                    init = 0
                    ter = 0
                    for ele in self.gff[entry.id][trans]['cds']:
                        temp = [int(ele[0])-beg, int(ele[1])-beg]
                        coord.append(temp)

                    if strand == '+':
                        coord = sorted(coord, key=lambda x: x[0])
                        init = coord[0][0]
                        ter = coord[-1][1]
                        for ele in coord:
                            inseq = ""
                            outseq = ""
                            if ele[0] != init:
                                inseq = seq[ele[0]-15:ele[0]+3]
                                self._update_entCDS(dict,observ,inseq)
                            else:
                                inseq = seq[ele[0]-5:ele[0]+7]
                                self._update_init(dict,observ,inseq)
                            if ele[1] != ter:
                                outseq = seq[ele[1]-2:ele[1]+8]
                                self._update_outCDS(dict,observ,outseq)
                            # else:
                            #     outseq = seq[ele[1]-6:ele[1]+6]
                            #     self._update_aaEnd(outseq)

                    elif strand == '-':
                        coord = sorted(coord, key=lambda x: -x[0])
                        init = coord[0][1]
                        ter = coord[-1][0]
                        for ele in coord:
                            inseq = ""
                            outseq = ""
                            if ele[1] != init:
                                inseq = seqt.complementary(seq[ele[1]-2:ele[1]+16])
                                self._update_entCDS(dict,observ,inseq)
                            else:
                                inseq = seqt.complementary(seq[ele[1]-6:ele[1]+6])
                                self._update_init(dict,observ,inseq)
                            if ele[0] != ter:
                                outseq = seqt.complementary(seq[ele[0]-7:ele[0]+3])
                                self._update_outCDS(dict,observ,outseq)
                            # else:
                            #     outseq = seqt.complementary(seq[ele[0]-5:ele[0]+7])
                            #     self._update_aaEnd(outseq)
        fas_read.close()
        return dict, observ

    # This is to store observations of how transcript initiated
    # Stuffs and Start Codon
    # Now wrong cases are gained by shift 1bp,
    # Maybe shift 3bp instead of 1?
    def _update_init(self, dict, observ, contest, prefix=True):
        if prefix:
            wrong1 = feature.Initiate_Site(contest[:-2])
            wrong2 = feature.Initiate_Site(contest[2:])
            correct = feature.Initiate_Site(contest[1:-1])
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

    def _update_entCDS(self, dict, observ, contest):
        wrong1 = feature.Enter_CDS(contest[:-2])
        wrong2 = feature.Enter_CDS(contest[2:])
        correct = feature.Enter_CDS(contest[1:-1])
        observ["entCDS"]["correct"].append(correct.out())
        observ["entCDS"]["wrong"].append(wrong1.out())
        observ["entCDS"]["wrong"].append(wrong2.out())
        try:
            dict["entCDS"]["end2"][correct.end2]+=1
            dict["entCDS"]["first1"][correct.first1]+=1
        except Exception as e:
            raise Read_error('What the Hell is ' +
                correct.end2 + ' or ' + correct.first1)

    def _update_outCDS(self, dict, observ, contest):
        wrong1 = feature.Out_CDS(contest[:-2])
        wrong2 = feature.Out_CDS(contest[2:])
        correct = feature.Out_CDS(contest[1:-1])
        observ["outCDS"]["correct"].append(correct.out())
        observ["outCDS"]["wrong"].append(wrong1.out())
        observ["outCDS"]["wrong"].append(wrong2.out())
        try:
            dict["outCDS"]["pre2"][correct.pre2]+=1
            dict["outCDS"]["first2"][correct.first2]+=1
            dict["outCDS"]["next4"][correct.next4]+=1
        except Exception as e:
            raise Read_error(correct.pre2, correct.first2, correct.next4)

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
