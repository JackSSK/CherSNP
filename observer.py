"""
Classes for preprocess data for learning part
Will have observations stored in a dict
And correspond dict to translate input seq to numbers
"""
import sys
import os
import re
import read_file as read
import gen_tools as t
import gffer
import feature
import math
import warnings as w

class Read_error(Exception):
    pass

class ID_error(Exception):
    pass

class Observer:
    def __init__(self, seq_file, gff_file):
        # Test to see how many init / term don't have contest
        self.init_count = 0
        self.term_count = 0
        # Process in GFF file
        self.gff = gffer.Process(gff_file).gff
        # Get observations and correspond dictionary
        self.dict, self.subj = self._observe(seq_file)
        self.dict = self._prepare_dict(self.dict)
        self.subj = self._obs_to_numb(self.subj)

        # print(self.init_count, self.term_count)


    def _observe(self, seq_file):
        #Initiating dicts for stroing correct observations
        dict = {
            "init":{
                "pre":t.generate(k=4, type='dict'),
                "start":t.generate(k=3, type='dict'),
                "aa1":t.generate(k=3, type='dict')
            },
            "term":{
                "last1":t.generate(k=3, type='dict'),
                "stop":t.generate(k=3, type='dict'),
                "next4":t.generate(k=4, type='dict'),
            },
            "entCDS":{
                "end2":t.generate(k=2, type='dict'),
                "first1":t.generate(k=1, type='dict'),
            },
            "outCDS":{
                "pre2":t.generate(k=2, type='dict'),
                "first2":t.generate(k=2, type='dict'),
                "next4":t.generate(k=4, type='dict'),
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
                if self.gff[entry.id][trans]["type"] == "mRNA":
                    if re.search(r'X',self.gff[entry.id][trans]["name"]):
                        continue
                    beg = self.gff[entry.id][trans]["beg"]
                    end = self.gff[entry.id][trans]["end"]
                    strand = self.gff[entry.id][trans]["strand"]
                    seq = entry.seq[beg-1:end]
                    coord = []
                    init = 0
                    ter = 0

                    qualify = True
                    for ele in self.gff[entry.id][trans]["cds"]:
                        if abs(ele[0] - ele[1]) < 3:
                            w.warn("Warning: "+ trans+ " has cds less than 3bp")
                            qualify = False
                            break
                        temp = [int(ele[0])-beg, int(ele[1])-beg]
                        coord.append(temp)

                    if not qualify: continue

                    if strand == "+":
                        coord = sorted(coord, key=lambda x: x[0])
                        init = coord[0][0]
                        ter = coord[-1][1]
                        for ele in coord:

                            # Is 5'UTR guranteed or not...

                            inseq = ""
                            outseq = ""
                            if ele[0] != init:
                                inseq = seq[ele[0]-18:ele[0]+4]
                                self._update_entCDS(dict,observ,inseq)
                            else:
                                inseq = seq[ele[0]-7:ele[0]+9]
                                if len(inseq) < 16:
                                    self.init_count += 1
                                    continue
                                self._update_init(dict,observ,inseq)
                            if ele[1] != ter:
                                outseq = seq[ele[1]-4:ele[1]+10]
                                self._update_outCDS(dict,observ,outseq)
                            else:
                                outseq = seq[ele[1]-8:ele[1]+8]
                                if len(outseq) < 16:
                                    self.term_count += 1
                                    continue
                                # Let's see who has the weird stop codon
                                # if outseq[3:-3][3:6] not in ['TAA','TGA','TAG']:
                                #     print(self.gff[entry.id][trans]["name"], coord)
                                self._update_term(dict, observ, outseq)

                    elif strand == '-':
                        coord = sorted(coord, key=lambda x: -x[0])
                        init = coord[0][1]
                        ter = coord[-1][0]
                        for ele in coord:
                            inseq = ""
                            outseq = ""
                            if ele[1] != init:
                                inseq = t.complementary(seq[ele[1]-3:ele[1]+19])
                                self._update_entCDS(dict,observ,inseq)
                            else:
                                inseq = t.complementary(seq[ele[1]-8:ele[1]+8])
                                if len(inseq) < 16:
                                    self.init_count += 1
                                    continue
                                # Let's see who has the weird start codon
                                # if inseq[7:10] != 'ATG':
                                #     print(self.gff[entry.id][trans]["name"], inseq[7:10], coord)
                                self._update_init(dict,observ,inseq)
                            if ele[0] != ter:
                                outseq = t.complementary(seq[ele[0]-9:ele[0]+5])
                                self._update_outCDS(dict,observ,outseq)
                            else:
                                outseq = t.complementary(seq[ele[0]-7:ele[0]+9])
                                if len(outseq) < 16:
                                    self.term_count += 1
                                    continue
                                # Let's see who has the weird stop codon
                                # if outseq[3:-3][3:6] not in ['TAA','TGA','TAG']:
                                #     print(outseq[3:-3][3:6])
                                self._update_term(dict, observ, outseq)
        fas_read.close()
        return dict, observ

    # This is to store observations of how transcript initiated
    # Stuffs and Start Codon
    # Now wrong cases are gained by shift 1bp,
    # Maybe shift 3bp instead of 1? Yes it is better wth 3bp
    def _update_init(self, dict, observ, contest, prefix=True):
        if prefix:
            wrong1 = feature.Initiate_Site(contest[:-6])
            wrong2 = feature.Initiate_Site(contest[6:])
            correct = feature.Initiate_Site(contest[3:-3])
            observ["init"]["correct"].append(correct.out())
            observ["init"]["wrong"].append(wrong1.out())
            observ["init"]["wrong"].append(wrong2.out())
            try:
                dict["init"]["start"][correct.start]+=1
                dict["init"]["aa1"][correct.first]+=1
                dict["init"]["pre"][correct.pre]+=1
            except Exception as e:
                print(contest)
                raise Read_error('What the Hell is ' +
                    correct.start + ' or ' + correct.first)

    def _update_entCDS(self, dict, observ, contest):
        wrong1 = feature.Enter_CDS(contest[:-6])
        wrong2 = feature.Enter_CDS(contest[6:])
        correct = feature.Enter_CDS(contest[3:-3])
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
        wrong1 = feature.Out_CDS(contest[:-6])
        wrong2 = feature.Out_CDS(contest[6:])
        correct = feature.Out_CDS(contest[3:-3])
        observ["outCDS"]["correct"].append(correct.out())
        observ["outCDS"]["wrong"].append(wrong1.out())
        observ["outCDS"]["wrong"].append(wrong2.out())
        try:
            dict["outCDS"]["pre2"][correct.pre2]+=1
            dict["outCDS"]["first2"][correct.first2]+=1
            dict["outCDS"]["next4"][correct.next4]+=1
        except Exception as e:
            raise Read_error(correct.pre2, correct.first2, correct.next4)

    def _update_term(self, dict, observ, contest, suffix=True):
        if suffix:
            wrong1 = feature.Term_Site(contest[:-6])
            wrong2 = feature.Term_Site(contest[6:])
            correct = feature.Term_Site(contest[3:-3])
            observ["term"]["correct"].append(correct.out())
            observ["term"]["wrong"].append(wrong1.out())
            observ["term"]["wrong"].append(wrong2.out())
            try:
                dict["term"]["last1"][correct.last1]+=1
                dict["term"]["stop"][correct.stop]+=1
                dict["term"]["next4"][correct.next4]+=1
            except Exception as e:
                raise Read_error(correct.last1, correct.stop, correct.next4)

    # This function will remove elements with 0 observations
    # Also, change # of observations to log(#observed / total observation)
    def _prepare_dict(self,dict):
        for type in dict:
            for sub in dict[type]:
                total = sum(dict[type][sub].values())
                del_list = []
                for ele in dict[type][sub]:
                    if dict[type][sub][ele] == 0: del_list.append(ele)
                    elif dict[type][sub][ele] < 0:
                        raise Read_error('Detected negative time observations')
                    else:
                        temp = dict[type][sub][ele]/total
                        dict[type][sub][ele] = math.log(temp)

                for ele in del_list:
                    del dict[type][sub][ele]
        return dict

    # Change observations from chars to doubles by using converted self.dict
    def _obs_to_numb(self,obs):
        for part in obs:
            if part == "init": names = ["pre", "start", "aa1"]
            elif part == "term": names = ["last1", "stop", "next4"]
            elif part == "outCDS": names = ["pre2", "first2", "next4"]
            elif part == "entCDS": names = ["DONE", "end2", "first1"]
            for sector in obs[part]:
                temp1 = []
                for term in obs[part][sector]:
                    if len(term) != len(names): print("error in _obs_to_numb")
                    temp2 = []
                    for i in range(len(term)):
                        ele = term[i]
                        if names[i] == "DONE":ele = ele
                        elif ele in self.dict[part][names[i]]:
                            ele = self.dict[part][names[i]][ele]
                        else:
                            ele = sys.float_info.min
                        temp2.append(ele)
                    temp1.append(temp2)
                obs[part][sector] = temp1
        return obs
