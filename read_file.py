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
        return self.seq

class GFF_entry:
    def __init__(self, line):
        content = line.split('\t')
        self.seqid = content[0]
        self.source = content[1]
        self.type = content[2]
        self.beg = int(content[3])
        self.end = int(content[4])
        self.score = content[5]
        self.strand = content[6]
        self.phase = content[7]
        self.id = ''
        self.pid = []
        self.cid = []
        if len(content) == 9:
            self.attr = content[8]
        else:
            self.attr = ''

class FASTA:
    def __init__(self, file):
        self.filename = file
        if re.search(r'\.gz$', self.filename):
            self._fas = gzip.open(file)
        else:
            self._fas = open(self.filename, 'r')
        self.entries = []
        self.itr_mark = 0
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
        if len(self.entries) == 0:
            raise Read_error('You sure this is FASTA?')

    def get(self, id):
        self._fas.seek(self.sta_coords[id])
        header = self._fas.readline()
        desc = header.split(id)[1].rstrip('\n')
        seq = []
        while(True):
            line = self._fas.readline()
            if line[0:1] == '>': break
            if line == '': break
            line = line.replace(' ', '')
            seq.append(line.strip())
        seq = ''.join(seq)
        return BioSeq(id, desc, seq)

    def __iter__(self):
        return self

    def __next__(self):
        if self.itr_mark == len(self.entries):
            raise StopIteration
        else:
            id = self.entries[self.itr_mark]
            self.itr_mark += 1
            return FASTA.get(self, id)

    def close(self):
        self._fas.close()

# May need to revise based on expect input format
class GFF:
    def __init__(self, file):
        self.filename = file
        if re.search(r'\.gz$', self.filename):
            self._gff = gzip.open(file)
        else:
            self._gff = open(self.filename, 'r')
        self.seqids = []
        self.itr_mark = 0
        self.sta_coords = {}
        while (True):
            line = self._gff.readline()
            if line[0:1]=='#': continue
            if line == '': break
            content = line.split("\t")
            if len(content) < 8: raise Read_error('Bad GFF Format')
            id = content[0]
            coord = self._gff.tell() - len(line) - 1
            if id not in self.seqids:
                self.seqids.append(id)
                self.sta_coords[id] = [coord]
            else:
                self.sta_coords[id].append(coord)
        if len(self.seqids) == 0:
            raise Read_error('You sure this is GFF?')

    def get(self, id):
        info = []
        for coord in self.sta_coords[id]:
            self._gff.seek(coord)
            entry = self._gff.readline()
            info.append(GFF_entry(entry))
        return info

    def __iter__(self):
        return self

    def __next__(self):
        if self.itr_mark == len(self.seqids):
            raise StopIteration
        else:
            id = self.seqids[self.itr_mark]
            self.itr_mark += 1
            return GFF.get(self, id)

    def close(self):
        self._gff.close()
