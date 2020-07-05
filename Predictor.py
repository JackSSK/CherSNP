import os
import re
import read_file as read
import seq_tools as seqt
import gffer

class ID_error(Exception):
    pass

class Predicter:
    def __init__(self, seq_file, gff_file):
        self.gff = gffer.Process(gff_file).gff
        # Test with json sample
        # self.gff = seqt.decode_json('temp_json.js')

        # There should be something to make sure HGVS annotations are integrated

        self.results = {}
        self._processFAS(seq_file)

    def _processFAS(self, seq_file):
        fas_read = read.FASTA(seq_file)
        for entry in fas_read:
            for record in self.gff[entry.id]:
                beg = self.gff[entry.id][record]['beg']
                end = self.gff[entry.id][record]['end']
                strand = self.gff[entry.id][record]['strand']
                type = self.gff[entry.id][record]['type']
                hgvs = self.gff[entry.id][record]['hgvs']

                if type == 'Gene' and hgvs is not None:
                    seq = entry.seq[beg-1:end]


                elif type == 'transcript' and hgvs is not None:
                    seq = entry.seq[beg-1:end]

                    coord = []
                    cds_contest = []
                    init = []
                    ter = []
                    ent_cds = []
                    ext_cds = []
                    for ele in self.gff[entry.id][record]['cds']:
                        temp = [int(ele[0])-beg, int(ele[1])-beg]
                        coord.append(temp)
                    if strand == '+':
                        coord = sorted(coord, key=lambda x: x[0])
                        init.append(coord[0][0])
                        ter.append(coord[-1][1])
                        for ele in coord:
                            enter = ele[0]
                            exit = ele[1]
                            if enter != init[0]: ent_cds.append(enter)
                            if exit != ter[0]: ext_cds.append(exit)
                            inseq = seqt.complementary(seq[enter-4:enter+6])
                            outseq = seqt.complementary(seq[exit-5:exit+5])
                            cds_contest.append([inseq, outseq])
                    elif strand == '-':
                        coord = sorted(coord, key=lambda x: -x[0])
                        init.append(coord[0][1])
                        ter.append(coord[-1][0])
                        for ele in coord:
                            enter = ele[1]
                            exit = ele[0]
                            if enter != init[0]: ent_cds.append(enter)
                            if exit != ter[0]: ext_cds.append(exit)
                            outseq = seqt.complementary(seq[exit-4:exit+6])
                            inseq = seqt.complementary(seq[enter-5:enter+5])
                            cds_contest.append([inseq, outseq])

                    print(coord, cds_contest, init, ent_cds, ext_cds, ter)
        fas_read.close()
