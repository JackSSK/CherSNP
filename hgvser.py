import re

# Class for pattern HGVS annotation
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

        elif re.search(r'ins',change):
            self.type = 'Insertion'
            self.info = substitution(change)

        elif re.search(r'del', change):
            self.type = 'Deletion'
            self.info = substitution(change)

        elif re.search(r'dup', change):
            self.type = 'Duplicatin'
            self.info = substitution(change)

        elif re.search(r'delins', change):
            self.type = 'Delet_Insertion'
            self.info = substitution(change)

        elif re.search(r'[', change):
            self.type = 'Repeat'
            self.info = substitution(change)

class substitution:
    def __init__(self, change):
        self.alt = change.split('>')[1]
        temp = re.findall(r'\d+', change)
        if len(temp)!=1:
            raise Change_error('No multiple position in substitution')
        self.pos = int(temp[0])
        self.ref = change.split('>')[0].split(temp[0])[1]

class insertion:
    def __init__(self, change):
        temp = change.split('ins')
        self.alt = temp[1]
        poss = temp[0].split('_')
        pos1 = poss[0]
        pos2 = poss[1]
        if pos2 == pos1 += 1: self.pos=pos1
        else:
            raise Change_error('Unqualified insertion variant')

class deletion:
    def __init__(self, change):
        temp = change.split('del')
        poss = temp[0].split('_')
        self.start = poss[0]
        self.end = poss[1]

class duplication:
    def __init__(self, change):
        temp = change.split('dup')
        poss = temp[0].split('_')
        self.start = poss[0]
        self.end = poss[1]

class del_insert:
    def __init__(self, change):
        temp = change.split('delins')
        self.alt = temp[1]
        poss = temp[0].split('_')
        self.start = poss[0]
        self.end = poss[1]

class repeated:
    def __init__(self, change):
        temp = change.split('[')
        self.times = temp[1].split(']')[0]
        pos = re.findall(r'\d+', temp[0])[0]
        self.pos = int(pos)
        self.seq = temp[0].split(pos)[0]
