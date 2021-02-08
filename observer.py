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
    def __init__(self, seq_file, gff_file, mode):

        # if mode == "noRNA":
        #     sites = [["init", 4, 3, 3], ["term", 3, 3, 4]]
        # elif mode == "hasRNA":
        #     sites = [["init", 4, 3, 3], ["term", 3, 3, 4], ["outCDS", 2,2,4], ["entCDS"]]

        sites = [["init", 4, 3, 3], ["term", 3, 3, 4], ["outCDS", 2,2,4], ["entCDS"]]

        self.dict, self.subj = self._initJeanne(sites)
        # Test to see how many init / term don't have contest
        self.init_count = 0
        self.term_count = 0
        # Process in GFF file
        self.gff = gffer.Process(gff_file,mode).gff
        # Get observations and correspond dictionary
        self.dict, self.subj = self._observe(seq_file,self.dict, self.subj, mode)
        self.dict = self._prepare_dict(self.dict)
        self.subj = self._obs_to_numb(self.subj)

        # print(self.init_count, self.term_count)

    def _initJeanne(self, sites):
        dict = {}
        observ = {}
        for site in sites:
            name = site[0]
            observ[name] = {
                "correct":[],
                "wrong":[]
            }
            if name == "entCDS":
                dict[name]= {
                    "end2":t.generate(k=2, type='dict'),
                    "first1":t.generate(k=1, type='dict')
                }
            else:
                dict[name] = {
                    "pre":t.generate(k=site[1], type='dict'),
                    "key":t.generate(k=site[2], type='dict'),
                    "suff":t.generate(k=site[3], type='dict')
                }

        return dict, observ


    def _observe(self, seq_file, dict, observ, mode):
        #Initiating dicts for stroing correct observations
        fas_read = read.FASTA(seq_file)
        for entry in fas_read:
            for trans in self.gff[entry.id]:
                # For GRCh38 training:
                if mode == 'hasRNA': types = ["mRNA"]
                # For cov19:
                elif mode == 'noRNA': types = ["mRNA", "Gene"]
                if self.gff[entry.id][trans]["type"] in types:
                    if re.search(r'X',self.gff[entry.id][trans]["name"]):
                        continue
                    beg = self.gff[entry.id][trans]["beg"]
                    end = self.gff[entry.id][trans]["end"]
                    strand = self.gff[entry.id][trans]["strand"]
                    # does not allow no-utr trecords
                    if mode == 'hasRNA': seq = entry.seq[beg-1:end]
                    # allow no-utr records
                    if mode == 'noRNA': seq = entry.seq
                    coord = []
                    init = 0
                    ter = 0

                    qualify = True
                    for ele in self.gff[entry.id][trans]["cds"]:
                        if abs(ele[0] - ele[1]) < 3:
                            w.warn("Warning: "+ trans+ " has cds less than 3bp")
                            qualify = False
                            break
                        if mode == 'hasRNA':
                            temp = [int(ele[0])-beg, int(ele[1])-beg]
                        if mode == 'noRNA':
                            temp = [int(ele[0]-1), int(ele[1]-1)]
                        coord.append(temp)
                    if not qualify: continue

                    if mode == 'noRNA' and len(coord)>1:
                        for ele1 in coord:
                            for ele2 in coord:
                                if ele1[1] == ele2[0]:
                                    new = [[ele1[0],ele2[1]]]
                                    coord.remove(ele1)
                                    coord.remove(ele2)
                                    if len(coord) > 0:
                                        self._getUpdates(coord, strand, seq, dict, observ)
                                    self._getUpdates(new, strand, seq, dict, observ)
                    else:
                        self._getUpdates(coord, strand, seq, dict, observ)

        fas_read.close()
        return dict, observ

    def _getUpdates(self, coord, strand, seq, dict, observ):
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
                    self._update_site(dict,observ,inseq, "init", 4,3,3)
                if ele[1] != ter:
                    outseq = seq[ele[1]-4:ele[1]+10]
                    self._update_site(dict,observ,outseq, "outCDS", 2,2,4)
                else:
                    outseq = seq[ele[1]-8:ele[1]+8]
                    if len(outseq) < 16:
                        self.term_count += 1
                        continue
                    # Let's see who has the weird stop codon
                    # if outseq[3:-3][3:6] not in ['TAA','TGA','TAG']:
                    #     print(self.gff[entry.id][trans]["name"], coord)
                    self._update_site(dict, observ, outseq, "term", 3,3,4)

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
                    self._update_site(dict,observ,inseq, "init", 4,3,3)
                if ele[0] != ter:
                    outseq = t.complementary(seq[ele[0]-9:ele[0]+5])
                    self._update_site(dict,observ,outseq, "outCDS", 2,2,4)
                else:
                    outseq = t.complementary(seq[ele[0]-7:ele[0]+9])
                    if len(outseq) < 16:
                        self.term_count += 1
                        continue
                    # Let's see who has the weird stop codon
                    # if outseq[3:-3][3:6] not in ['TAA','TGA','TAG']:
                    #     print(outseq[3:-3][3:6])
                    self._update_site(dict, observ, outseq, "term", 3,3,4)

    # This is to store observations of how transcript initiated
    # Stuffs and Start Codon
    # Now wrong cases are gained by shift 1bp,
    # Maybe shift 3bp instead of 1? Yes it is better wth 3bp

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

    def _update_site(self, dict, observ, seq, name,
        preLen, keyLen, suffLen,
        context = True):
        if context:
            wrong1 = feature.Site(seq[:-6] ,preLen, keyLen, suffLen)
            wrong2 = feature.Site(seq[6:],preLen, keyLen, suffLen)
            correct = feature.Site(seq[3:-3],preLen, keyLen, suffLen)
            observ[name]["correct"].append(correct.out())
            observ[name]["wrong"].append(wrong1.out())
            observ[name]["wrong"].append(wrong2.out())
            try:
                dict[name]["pre"][correct.pre]+=1
                dict[name]["key"][correct.key]+=1
                dict[name]["suff"][correct.suff]+=1
            except Exception as e:
                raise Read_error(correct.pre, correct.key, correct.suff)

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
            if part == "entCDS": names = ["DONE", "end2", "first1"]
            else: names = ["pre", "key", "suff"]
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
