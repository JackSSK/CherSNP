"""
Classes for read in files
"""
import os
import gzip
import re

class Read_error(Exception):
    pass

class BioSeq:
    def __init__(self, id, description, seq):
        self.id = id
        self.description = description
        self.seq = seq
    def hgvsMod(self, hgvs=''):
        # Return a bioseq object after modification based on
        # hgvs annotation
        return self.seq

# Class to read in a FASTA sequence file
class FASTA:
    def __init__(self, file):
        self.filename = file
        # Open as .gz file
        if re.search(r'\.gz$', self.filename):
            self.file = gzip.open(file)
        # Open directly
        else:
            self.file = open(self.filename, 'r')
        self.entries = []
        self.iteration = 0
        self.entryStartCoordinates = {}
        # Read through lines
        while (True):
            line = self.file.readline()
            if line == '': break
            # Detect an entry
            if line[0:1] == '>':
                # Get ID
                id = re.search(r'>\s*(\S+)', line)[1]
                if id in self.entries:
                    raise Read_error('Duplicate id: ' + id)
                else: self.entries.append(id)
                self.entryStartCoordinates[id] = self.file.tell() - len(line)
        if len(self.entries) == 0:
            raise Read_error('You sure this is FASTA?')
    # Get specific sequence
    def get(self, id):
        self.file.seek(self.entryStartCoordinates[id])
        header = self.file.readline()
        description = header.split(id)[1].rstrip('\n')
        seq = []
        # Iterate through lines till next entry
        # or enconter an empty line
        while(True):
            line = self.file.readline()
            if line[0:1] == '>': break
            if line == '': break
            line = line.replace(' ', '')
            seq.append(line.strip())
        seq = ''.join(seq)
        return BioSeq(id, description, seq)

    def __iter__(self):
        return self

    def __next__(self):
        if self.iteration == len(self.entries):
            raise StopIteration
        else:
            id = self.entries[self.iteration]
            self.iteration += 1
            return FASTA.get(self, id)

    def close(self):
        self.file.close()

# May need to revise based on expect input format
# Class to read in a GFF file
class GFF:
    def __init__(self, file):
        self.filename = file
        if re.search(r'\.gz$', self.filename):
            self.file = gzip.open(file)
        else:
            self.file = open(self.filename, 'r')
        self.seqIDs = []
        self.iteration = 0
        self.entryStartCoordinates = {}
        while (True):
            line = self.file.readline()
            if line[0:1]=='#': continue
            if line == '': break
            content = line.split("\t")
            if len(content) < 8:
                if content == ["\n"]: continue
                raise Read_error('Bad GFF Format')
            id = content[0]

            # Something goes wrong here when test with TTC5.gff or Chr14.gff
            # -1 work with TTC5, and no - 1 wokrs with Chr14
            coordinate = self.file.tell() - len(line)

            if id not in self.seqIDs:
                self.seqIDs.append(id)
                self.entryStartCoordinates[id] = [coordinate]
            else:
                self.entryStartCoordinates[id].append(coordinate)
        if len(self.seqIDs) == 0:
            raise Read_error('You sure this is GFF?')

    def get(self, id):
        info = []
        for coordinate in self.entryStartCoordinates[id]:
            self.file.seek(coordinate)
            entry = self.file.readline()
            info.append(GFF_entry(entry))
        return info

    def __iter__(self):
        return self

    def __next__(self):
        if self.iteration == len(self.seqIDs):
            raise StopIteration
        else:
            id = self.seqIDs[self.iteration]
            self.iteration += 1
            return GFF.get(self, id)

    def close(self):
        self.file.close()

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
