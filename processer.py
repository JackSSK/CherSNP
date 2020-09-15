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
            self._doEntry(entry)
        fas_read.close()

    # Process each entry
    def _doEntry(self, entry):
            pos = 0
            cds_num = 1
            status = 'utr'
            mark = ''
            cds = []
            enter = 0
            while status != 'term':
                if status == 'utr':
                    enter, mark = self._inUTR(entry, pos, mark)
                    if enter == -1:
                        print('No CDS')
                    else:
                        status = 1
                elif isinstance(status, int):
                    # Test now
                    print(mark, enter, entry.seq)
                    status = 'term'

    # Process while not findin initial site
    def _inUTR(self, entry, pos, mark):
        while True:
            # What if no Init site found
            if len(mark) == len(entry.seq): return -1, mark

            # Get sequence in current posistion and processed
            seq = entry.seq[pos: pos+10]
            fea = f.Initiate_Site(seq)

            # Get according numbs from dict
            if fea.pre in self.clfs.dict["init"]["pre"]:
                pre = self.clfs.dict["init"]["pre"][fea.pre]
            else: pre =sys.float_info.min

            if fea.start in self.clfs.dict["init"]["start"]:
                start = self.clfs.dict["init"]["start"][fea.start]
            else: start =sys.float_info.min

            if fea.first in self.clfs.dict["init"]["aa1"]:
                first = self.clfs.dict["init"]["aa1"][fea.first]
            else: first =sys.float_info.min

            obs = [pre, start, first]
            pred = self.clfs.init.predict([obs])
            if pred[0] == 0:
                mark += 'u'
                pos += 1

            elif pred[0] == 1:
                pos += 4
                mark += 'uuuu'
                return pos, mark
