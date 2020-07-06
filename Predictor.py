import os
import re
import read_file as read
import seq_tools as seqt
import gffer
import hgvser

class ID_error(Exception):
    pass

class Type_error(Exception):
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
                # strand = self.gff[entry.id][record]['strand']
                type = self.gff[entry.id][record]['type']
                hgvs = self.gff[entry.id][record]['hgvs']

                if type == 'Gene' and len(hgvs)>0:
                    seq = entry.seq[beg-1:end]


                elif type == 'transcript' and len(hgvs)>0:
                    seq = entry.seq[beg-1:end]

    # In this level check for whether the consequence would be:
    # exon_variant: UTR variant, coding_sequence_variant
    # intron_variant:
    def _trans_1st_classify(trans, hgvs):
        beg = trans['beg']
        end = trans['end']
        strand = trans['end']
        anno = hgvser.HGVS(hgvs)
        prefix = anno.prefix
        if anno.type != 'Substitution':
            raise Type_error('HGVS in 1st trans clssifier is not substitution')
        if prefix != 'RNA':
            raise Type_error('HGVS in 1st trans classifier is not r')
        pos = anno.info.pos
        ref = anno.info.ref
        alt = anno.info.alt
