import gen_tools as t
STOP_CODONS = ['TGA', 'TAA', 'TAG']

# This function/class is for validating predicted cds pattern
def validate(seq, STOP_CODONS):
    if len(seq) % 3 > 0:
        return False
    else:
        i = 0
        while i < len(seq)-3:
            j= i+3
            read = seq[i:j]
            if read in STOP_CODONS:
                return False
            i+=3
        if seq[i:] in STOP_CODONS:
            return True


class Tranlate:
    def __init__(self, seq, patterns, snps):
        self.results = {}
        non_coding = []
        for ele in snps:
            non_coding.append(ele)
            self.results[ele.origin] = {
                'hgvs':ele.origin,
                'predicts':[]
            }
        for pattern in patterns:
            id = pattern["id"]
            regions = pattern["patterns"]
            tests = []
            for ele in regions:
                temp = [ele, []]
                for snp in snps:
                    if ele[0] < snp.info.pos < ele[1]:
                        temp[1].append(snp)
                        if snp in non_coding: non_coding.remove(snp)
                if len(temp[1]) > 0:
                    tests.append(temp)
            for ele in tests:
                region = ele[0]
                self._translate(seq, region, ele[1])

        for ele in non_coding:
            self.results[ele.origin]['predicts'] = ['None']
        t.encode_json(self.results)

    def _translate(self, seq, region, snps):
        start = region[0]
        end = region[1]
        alerts = {}
        for ele in snps:
            pos = ele.info.pos
            ref = ele.info.ref
            alt = ele.info.alt
            if pos not in alerts:
                alerts[pos] = [ele]
            else:
                alerts[pos].append(ele)
        table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }
        if len(seq[start:end])%3 == 0:
            for i in range(start, end, 3):
                for ele in alerts:
                    pos = ele-1
                    if pos in [i,i+1,i+2]:
                        codon = seq[i:i + 3]
                        for snp in alerts[ele]:
                            alt = snp.info.alt
                            if pos - i == 0:
                                new = alt + codon[1:]
                            elif pos - i == 1:
                                new = codon[0] + alt + codon[2]
                            else:
                                new = codon[:2] + alt
                            ori = table[codon]
                            new = table[new]
                            pred = "p." + ori + str(int(1+(i-start)/3)) + new
                            self.results[snp.origin]['predicts'].append(
                                [pred, [start,end]]
                            )
