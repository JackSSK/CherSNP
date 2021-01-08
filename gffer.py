import os
import re
import numpy as np
import read_file as read
import gen_tools as tool

class ID_error(Exception):
    pass

class Process:
    def __init__(self, file, mode):
        self.gff = {}
        self.mode = mode
        self._allway(file)
        # For test
        # tool.encode_json(self.gff)

    def _allway(self, file):
        gff_read = read.GFF(file)
        unlinkCDS = np.array([])
        for srcseq in gff_read:
            for entry in srcseq:
                if entry.seqid not in self.gff:
                    self.gff[entry.seqid] = {}

                if entry.type =='gene':
                    self._doGene(entry)

                elif entry.type == 'mRNA':
                    self._doRNA(entry)

                elif entry.type == 'CDS':
                    self._doCDS(entry, unlinkCDS)
        gff_read.close()
        self._doUnlinkCDS(unlinkCDS)


    def _doGene(self, entry):
        if entry.type =='gene':
            if re.search(r'ID=',entry.attr):
                id = entry.attr.split('ID=')[1].split(';')[0]
                hgvs = []

                #If HGVS is already annotated in gff
                if re.search(r'HGVS=',entry.attr):
                    var = entry.attr.split('HGVS=')[1].split(';')[0]
                    hgvs.append(var)

                if id not in self.gff[entry.seqid]:
                    self.gff[entry.seqid][id] = {
                        'type':'Gene',
                        'beg':int(entry.beg),
                        'end':int(entry.end),
                        'strand':entry.strand,
                        'name':id,
                        'cds':[],
                        'hgvs':hgvs
                    }
                else:
                    raise ID_error('mRNA read in twice')
            else:
                raise ID_error('ID missing')

    def _doRNA(self, entry):
        # Parent of mRNA is not mandatory
        if re.search(r'Parent=',entry.attr):
            parent = entry.attr.split('Parent=')[1].split(';')[0]
        else: parent = None

        if re.search(r'ID=',entry.attr):
            id = entry.attr.split('ID=')[1].split(';')[0]
            hgvs = []

            name = entry.attr.split('Name=')[1].split(';')[0]

            #If HGVS is already annotated in gff
            if re.search(r'HGVS=',entry.attr):
                var = entry.attr.split('HGVS=')[1].split(';')[0]
                hgvs.append(var)

            if id not in self.gff[entry.seqid]:
                self.gff[entry.seqid][id] = {
                    'type':'mRNA',
                    'parent':parent,
                    'beg':int(entry.beg),
                    'end':int(entry.end),
                    'strand':entry.strand,
                    'hgvs':hgvs,
                    'name':name,
                    'cds':[],
                }
            else:
                raise ID_error('mRNA read in twice')
        else:
            raise ID_error('ID missing')

    def _doCDS(self, entry, unlinkCDS):
        if re.search(r'Parent=',entry.attr):
            pid = entry.attr.split('Parent=')[1].split(';')[0]
            # Assuming CDS has to be a feature of a transcript
            # Not in Cov19 case
            if re.search('gene', pid) and self.mode == 'hasRNA':
                return
            # Ignore single nucleotide entry
            if int(entry.beg) == int(entry.end):
                return
            if pid not in self.gff[entry.seqid]:
                cds = [entry.seqid, pid, int(entry.beg),
                    int(entry.end)]
                unlinkCDS = np.append(unlinkCDS, cds)
            else:
                ele = [int(entry.beg), int(entry.end)]
                if ele not in self.gff[entry.seqid][pid]['cds']:
                    self.gff[entry.seqid][pid]['cds'].append(ele)
                else:
                    raise ID_error('CDS already exist')
        else:
            raise ID_error('Parent ID missing')

    def _doUnlinkCDS(self, unlinkCDS):
        for cds in unlinkCDS:
            seqid = cds[0]
            pid = cds[1]
            ele = [cds[2], cds[3]]
            try:
                if pid in self.gff[seqid]:
                    if ele not in self.gff[seqid][pid]['cds']:
                        self.gff[seqid][pid]['cds'].append(ele)
                    else:
                        raise ID_error('CDS already exist')
            except:
                raise ID_error('Cannot load unlink CDS')
