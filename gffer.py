import os
import re
import read_file as read
import gen_tools as tool

class ID_error(Exception):
    pass

class Process:
    def __init__(self, file):
        self.gff = {}
        self._allway(file)
        # For test
        # tool.encode_json(self.gff)

    def _allway(self, file):
        gff_read = read.GFF(file)
        unlink_cds = []
        for srcseq in gff_read:
            for entry in srcseq:
                if entry.seqid not in self.gff:
                    self.gff[entry.seqid] = {}

                if entry.type =='gene':
                    if re.search(r'ID=',entry.attr):
                        id = entry.attr.split('ID=')[1].split(';')[0]
                        hgvs = []

                        #If HGVS is already annotated in gff
                        if re.search(r'HGVS=',entry.attr):
                            hgvs.append(entry.attr.split('HGVS=')[1].split(';')[0])

                        if id not in self.gff[entry.seqid]:
                            self.gff[entry.seqid][id] = {
                                'type':'Gene',
                                'beg':int(entry.beg),
                                'end':int(entry.end),
                                'strand':entry.strand,
                                'hgvs':hgvs
                            }
                        else:
                            raise ID_error('mRNA read in twice')
                    else:
                        raise ID_error('ID missing')
                        continue

                elif entry.type == 'mRNA':
                    # Parent of mRNA is not mandatory
                    if re.search(r'Parent=',entry.attr):
                        parent = entry.attr.split('Parent=')[1].split(';')[0]
                    else: parent = None

                    if re.search(r'ID=',entry.attr):
                        id = entry.attr.split('ID=')[1].split(';')[0]
                        hgvs = []

                        #If HGVS is already annotated in gff
                        if re.search(r'HGVS=',entry.attr):
                            hgvs.append(entry.attr.split('HGVS=')[1].split(';')[0])

                        if id not in self.gff[entry.seqid]:
                            self.gff[entry.seqid][id] = {
                                'type':'transcript',
                                'parent':parent,
                                'beg':int(entry.beg),
                                'end':int(entry.end),
                                'strand':entry.strand,
                                'hgvs':hgvs,
                                'cds':[],
                            }
                        else:
                            raise ID_error('mRNA read in twice')
                    else:
                        raise ID_error('ID missing')
                        continue

                elif entry.type == 'CDS':
                    if re.search(r'Parent=',entry.attr):
                        pid = entry.attr.split('Parent=')[1].split(';')[0]
                        if pid not in self.gff[entry.seqid]:
                            unlink_cds.append([entry.seqid, pid, int(entry.beg), int(entry.end)])
                        else:
                            ele = [int(entry.beg), int(entry.end)]
                            if ele not in self.gff[entry.seqid][pid]['cds']:
                                self.gff[entry.seqid][pid]['cds'].append(ele)
                            else:
                                raise ID_error('CDS already exist')
                    else:
                        raise ID_error('Parent ID missing')
                        continue
        gff_read.close()
        for cds in unlink_cds:
            seqid = cds[0]
            pid = cds[1]
            ele = [cds[2], cds[3]]
            try:
                if ele not in self.gff[seqid][pid]['cds']:
                    self.gff[seqid][pid]['cds'].append(ele)
                else:
                    raise ID_error('CDS already exist')
            except:
                raise ID_error('Cannot load unlink CDS')
