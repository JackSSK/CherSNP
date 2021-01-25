"""
Classes to find cds patterns in input sequences with trained classifiers
"""

import feature as f
import read_file as read
import translator as trans
import sys
import gen_tools as t
import itertools

class Finder:
    def __init__(self, seq_file, clfs):
        self.clfs = clfs
        self.patterns = []
        # Read in sequences
        fas_read = read.FASTA(seq_file)
        # Process sequence 1 by 1
        for entry in fas_read:
            print("Processing " + entry.id)
            pattern = self._doEntry(entry,10,1)
            self.patterns.append({
                'id': entry.id,
                'patterns':pattern
            })
        fas_read.close()


    # Process each entry
    def _doEntry(self, entry, max_intron, ratio):
        init = []
        pos = 0
        outCDS = []
        entCDS = []
        term = []

        # Preprocess seq to U->T lower to upper
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
                score = outCDS_result[1]
                outCDS.append([pos+2, score])
            if len(seq_len16) == 16 and entCDS_result[0] == 1:
                score = entCDS_result[1]
                entCDS.append([pos+15, score])

            pos += 1
        # Sort positions based on scores in decreasing order
        init = self._partial(init, ratio)
        term = self._partial(term, ratio)
        entCDS = self._partial(entCDS, ratio)
        outCDS = self._partial(outCDS, ratio)
        print(len(init), len(term), len(entCDS), len(outCDS))
        return self._getPatterns(entry, init, term, entCDS, outCDS)

    # Anze
    # Also get potential combinations of start and stop sites
    # and try to find possible CDS patterns with these locations

    def _getSet(self, start, end, outCDS, entCDS):
        entCDS = self._getSites(start[0], end[0], entCDS)
        outCDS = self._getSites(start[0], end[0], outCDS)
        patterns = []
        self._getSetDfs(0, 0, entCDS, outCDS, len(entCDS), len(outCDS), [], patterns)
        return patterns

    def _getSetDfs(self, entCDS_index, outCDS_index, entCDS, outCDS, entCDSLen, outCDSLen, ans, ansSet):
        if entCDS_index >= entCDSLen or outCDS_index >= outCDSLen:
            ansSet.append(ans)
            # print(len(ansSet), ans)
            return
        for i in range(outCDS_index, outCDSLen):
            if outCDS[i][0] > entCDS[entCDS_index][0]:
                nextOutCDSIndex = self._getNextA(entCDS, outCDS, entCDS_index, i)
                self._getSetDfs(nextOutCDSIndex, i+1, entCDS, outCDS, entCDSLen,
                    outCDSLen, ans + [[entCDS[entCDS_index][0],
                    outCDS[i][0]]], ansSet)
        self._getSetDfs(entCDS_index+1, outCDS_index, entCDS, outCDS, entCDSLen,
            outCDSLen, ans, ansSet)

    def _getNextA(self, listA, listB, indexA, indexB):
        while(listA[indexA] <= listB[indexB]):
            indexA = indexA + 1
            if(indexA >= len(listA)):
                break
        return indexA


    def _getPatterns(self, entry, init, term, entCDS, outCDS):
        patterns = []
        for start in init:
            for end in term:
                start = [161, -12312]
                end = [3912, -12312]
                if start[0] < end[0]:
                    score = start[1] + end[1]
                    # If it does not surely need an intron to be validated
                    stopers = ['TGA', 'TAA', 'TAG']
                    if trans.validate(entry.seq[start[0]:end[0]],stopers):
                        patterns.append({
                            'start':start[0],
                            'end':end[0],
                            'introns':[],
                            'score':score
                        })
                    print("Trace On!")
                    temp = self._getSet(start, end, outCDS, outCDS)
                    print("Unlimited Blade Works!")
                    for introns in temp:
                        patterns.append({
                          'start':start[0],
                          'end':end[0],
                          'introns':introns,
                          'score':score
                        })
                break
        for ele in patterns:
            if ele['start'] == 161:
                print(ele)
        return patterns

    # Get intron intervals which have containing introns overlap with each other
    def _getIntervals(self, introns):
        intervals = []
        for intron in introns:
            int_start = intron[0][0]
            int_end = intron[0][1]
            next = True
            for ele in intervals:
                mini = ele[0]
                maxi = ele[1]
                if mini <= int_start <= maxi or int_start <= mini <= int_end:
                    ele[0] = max(mini, int_start)
                    ele[1] = min(maxi, int_end)
                    ele[2].append(intron)
                    next = False
                    break
            if next:
                intervals.append([int_start, int_end, [intron]])
        answer = []
        for ele in intervals:
            answer.append(ele[2])

        return answer

    # Based on outCDS and entCDS, get introns that may exist
    # Collapese similar results
    def _getPotentialIntrons(self, entCDS, outCDS):
        pos_introns = []
        for donor in outCDS:
            for acceptor in entCDS:
                if donor[0] < acceptor[0]-18:
                    pair = [donor[0], acceptor[0]]
                    score = donor[1] + acceptor[1]
                    pos_introns.append([pair,score])
        # sort possible combinations by starting position
        pos_introns = sorted(pos_introns, key = lambda x: x[0][0],reverse=True)
        # collapese similar predictions which have boundary location differences
        # within 3bp
        i = 0
        answer = []
        while i < len(pos_introns)-1:
            ele1 = pos_introns[i]
            ele2 = pos_introns[i+1]
            if abs(ele1[0][0]-ele2[0][0]) + abs(ele1[0][1] - ele2[0][1]) <= 3:
                if ele1[1] > ele2[1]:
                    answer.append(ele1)
                    i += 2
                else:
                    answer.append(ele2)
                    i += 1
            else:
                i += 1
                answer.append(ele1)
        return answer

    # Get potential donor/acceptor sites in range of start to end:
    def _getSites(self, start, end, set):
        set = sorted(set, key=lambda x: x[0])
        answer = []
        for ele in set:
            if start < ele[0] < end:
                answer.append(ele)
            elif ele[0] >= end: break
        return answer

    # Partial a coord set by score and a desired ratio
    def _partial(self, set, ratio):
        return sorted(set, key = lambda x: x[1],
            reverse=True)[:int(len(set)*ratio)]


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
