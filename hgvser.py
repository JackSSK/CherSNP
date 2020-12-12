import re

# Now it is only for substitution of gene or transcript
# No deletion or insertion etc
# No allele annotation (changes seperate with ; )

class Prefix_error(Exception):
    pass

class Change_error(Exception):
    pass

class HGVS:
    def __init__(self, anno):
        prefix = anno.split('.')[0]
        change = anno.split('.')[1]
        self.origin = anno
        if prefix == 'g': self.prefix = 'Genomic'
        elif prefix == 'c': self.prefix = 'coding Gene'
        elif prefix == 'r': self.prefix = 'RNA'
        else: raise Prefix_error(prefix + 'is not an acceptable HGVS prefix')

        if re.search(r';',change):
            raise Change_error('No multiple changes acceptable now')

        if re.search(r'>',change):
            self.type = 'Substitution'
            self.info = substitution(change)

class substitution:
    def __init__(self, change):
        self.alt = change.split('>')[1]
        temp = re.findall(r'\d+', change)
        if len(temp)!=1:
            raise Change_error('No multiple position in substitution')
        self.pos = int(temp[0])
        self.ref = change.split('>')[0].split(temp[0])[1]
