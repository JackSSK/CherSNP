"""
Classes to process input sequences with trained classifiers
"""

import feature as f
import read_file as read
import sys

class Processor:
    def __init__(self, seq_file, clfs):
        self.clfs = clfs
        # Read in sequences
        fas_read = read.FASTA(seq_file)
        # Process sequence 1 by 1
        for entry in fas_read:
            print("Processing " + entry.id)
            temp = self._doEntry(entry)
            # print(entry.seq[10:20], temp)
        fas_read.close()

    # Process each entry
    def _doEntry(self, entry):
        init = []
        pos = 0
        outCDS = []
        entCDS = []
        term = []
        while pos < len(entry.seq):
            # Get sequence in current posistion and processed
            seq_len10 = entry.seq[pos: pos+10]
            seq_len8 = entry.seq[pos: pos+8]
            seq_len16 = entry.seq[pos: pos+16]

            # stop when not enough length
            if len(seq_len10) < 10: break

            # Make predicions
            init_result = self._isInit(seq_len10)
            term_result = self._isTerm(seq_len10)
            outCDS_result = self._isOutCDS(seq_len8)
            entCDS_result = self._isEntCDS(seq_len16)
            # Check predictions
            if init_result[0] == 1:
                score = init_result[1]
                init.append([pos+4, score])
            if term_result[0] == 1:
                score = term_result[1]
                term.append([pos+6, score])
            if outCDS_result[0] == 1:
                outCDS.append([pos+2, score])
            if len(seq_len16) == 16 and entCDS_result[0] == 1:
                entCDS.append([pos+15, score])

            pos += 1
        # Sort positions based on scores in decreasing order
        init = sorted(init, key = lambda x: x[1],
            reverse=True)[:int(len(init)/4)]
        term = sorted(term, key = lambda x: x[1],
            reverse=True)[:int(len(term)/4)]
        entCDS = sorted(entCDS, key = lambda x: x[1],
            reverse=True)[:int(len(entCDS)/4)]
        outCDS = sorted(outCDS, key = lambda x: x[1],
            reverse=True)[:int(len(outCDS)/4)]
        print(len(init), "\n", init, "\n",
            len(term), "\n", term, "\n", len(entCDS), len(outCDS),)
        # return answers

    # Check whether a potential init site exist or not
    def _isInit(self, seq):
        fea = f.Initiate_Site(seq)

        # Get according numbs from dict
        if fea.pre in self.clfs.dict["init"]["pre"]:
            pre = self.clfs.dict["init"]["pre"][fea.pre]
        else: return [False]

        if fea.start in self.clfs.dict["init"]["start"]:
            start = self.clfs.dict["init"]["start"][fea.start]
        else: return [False]

        if fea.first in self.clfs.dict["init"]["aa1"]:
            first = self.clfs.dict["init"]["aa1"][fea.first]
        else: return [False]

        obs = [pre, start, first]

        # Do prediction only if observations are all valid
        if sys.float_info.min not in obs:
            pred = self.clfs.init.predict([obs])
        else:
            pred = [0]

        # Check predicions
        if pred[0] == 1: return [True,sum(obs)]
        else: return [False]

    # Check whether a potential term site exist or not
    def _isTerm(self, seq):
        fea = f.Term_Site(seq)

        # Get according numbs from dict
        if fea.last1 in self.clfs.dict["term"]["last1"]:
            last1 = self.clfs.dict["term"]["last1"][fea.last1]
        else: return [False]

        if fea.stop in self.clfs.dict["term"]["stop"]:
            stop = self.clfs.dict["term"]["stop"][fea.stop]
        else: return [False]

        if fea.next4 in self.clfs.dict["term"]["next4"]:
            next4 = self.clfs.dict["term"]["next4"][fea.next4]
        else: return [False]

        obs = [last1, stop, next4]

        # Do prediction only if observations are all valid
        if sys.float_info.min not in obs:
            pred = self.clfs.term.predict([obs])
        else:
            pred = [0]

        # Check predicions
        if pred[0] == 1: return [True,sum(obs)]
        else: return [False]

    # Check whether a potential donor site exist or not
    def _isOutCDS(self, seq):
        fea = f.Out_CDS(seq)

        # Get according numbs from dict
        if fea.pre2 in self.clfs.dict["outCDS"]["pre2"]:
            pre2 = self.clfs.dict["outCDS"]["pre2"][fea.pre2]
        else: return [False]

        if fea.first2 in self.clfs.dict["outCDS"]["first2"]:
            first2 = self.clfs.dict["outCDS"]["first2"][fea.first2]
        else: return [False]

        if fea.next4 in self.clfs.dict["outCDS"]["next4"]:
            next4 = self.clfs.dict["outCDS"]["next4"][fea.next4]
        else: return [False]

        obs = [pre2, first2, next4]

        # Do prediction only if observations are all valid
        if sys.float_info.min not in obs:
            pred = self.clfs.outCDS.predict([obs])
        else:
            pred = [0]

        # Check predicions
        if pred[0] == 1: return [True,sum(obs)]
        else: return [False]

    # Check whether a potential acceptor site exist or not
    def _isEntCDS(self, seq):
        fea = f.Enter_CDS(seq)

        # Get according numbs from dict
        if fea.end2 in self.clfs.dict["entCDS"]["end2"]:
            end2 = self.clfs.dict["entCDS"]["end2"][fea.end2]
        else: return [False]

        if fea.first1 in self.clfs.dict["entCDS"]["first1"]:
            first1 = self.clfs.dict["entCDS"]["first1"][fea.first1]
        else: return [False]

        obs = [fea.Yratio, end2, first1]

        # Do prediction only if observations are all valid
        if sys.float_info.min not in obs:
            pred = self.clfs.entCDS.predict([obs])
        else:
            pred = [0]

        # Check predicions
        if pred[0] == 1: return [True,sum(obs)]
        else: return [False]
