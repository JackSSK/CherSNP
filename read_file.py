"""
Classes for read in files
"""

import os
import gzip
import re

class Read_error(Exception):
	pass

class BioSeq:
    def __init__(self, id, desc, seq):
        self.id = id
        self.desc = desc
        self.seq = seq
    def hgvsMod(self, hgvs=''):
        # Return a bioseq object after modification based on
        # hgvs annotation
        return BioSeq(id,desc,seq)

class FASTA:
    def __init__(self, file):
        self.filename = file
        if re.search(r'\.gz$', self.filename):
            self._fas = gzip.open(file)
        else:
            self._fas = open(self.filename, 'r')
        self.entries = []
        self.sta_coords = {}
        while (True):
            line = self._fas.readline()
            if line == '': break
            if line[0:1] == '>':
                m = re.search(r'>\s*(\S+)', line)
                id = m[1]
                if id in self.entries:
                    raise Read_error('duplicate id: ' + id)
                else: self.entries.append(id)
                self.sta_coords[id] = self._fas.tell() - len(line)
        self._fas.close()
        if len(self.entries) == 0:
            raise Read_error('You sure this is FASTA?')

    def test(self):
        return(self.entries)

    # def sequence(self, id):
    #     self._fas.seek(self.dict[id]['pos'], os.SEEK_SET)
    #     descibtion = self._fas.readline()
    #     seq = []
    #     while(True):
    #         pos = self._fas.tell()
    #         if pos == self.end:
    #             break
    #         else:
    #             line = self._fas.readline()
    #             content = line.split()
    #             if(len(content) == 0):
    #                 continue
    #             elif(len(content) == 1):
    #                 seq.append(content[0])
    #             elif(len(content) > 1):
    #                 break
    #     seq = ''.join(seq)
    #     return seq
    # def describe(self, id):
    #     self._fas.seek(self.dict[id]['pos'], os.SEEK_SET)
    #     descibtion = self._fas.readline()
    #     return descibtion
    # def num(self, id):
    #     return self.dict[id]['num']
    # def dict(self):
    #     return self.dict
    # def close(self):
    #     self._fas.close()

def gff_read(file):
    """Readin *.gff file entirely
    return a dictionary"""
    dict = {}
    for line in file:
        content = line.split("\t")
        seq_name = content[0]
        src = content[1]
        feature = content[2]
        start = int(content[3])
        end = int(content[4])
        score = content[5]
        strand = content[6]
        phase = content[7]
        atr = content[8]
        if src not in dict:
            dict[src] = {}
        if seq_name not in dict[src]:
            dict[src][seq_name] = {}
        if strand not in dict[src][seq_name]:
            dict[src][seq_name][strand] = {}
        if feature not in dict[src][seq_name][strand]:
            dict[src][seq_name][strand][feature] = []
        info = [start, end, score, phase, atr]
        dict[src][seq_name][strand][feature].append(info)
    return dict


print(FASTA('data/foo.fa').test())
