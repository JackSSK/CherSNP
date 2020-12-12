import observer as obs
import argparse
import hgvser
import trainer as train
import pattern
import process as p
import gffer
import translator as t
import read_file as read
import hgvser as h

parser = argparse.ArgumentParser(description='Just for test now')
parser.add_argument('--fasta', required=True, type=str,
    metavar='<path>', help='path to fasta file for query sequences')
parser.add_argument('--path', required=False, type=str,default='default/',
    metavar='<path>', help='path to classifier (.pkl) files and dictionary')
parser.add_argument('--dict', required=False, type=str,default='Jeanne.js',
    metavar='<path>', help='path to dict(.js) file')
parser.add_argument('--init', required=False, type=str,default='Katyusha.pkl',
    metavar='<path>', help='path to init classifier (.pkl) file')
parser.add_argument('--term', required=False, type=str,default='Erika.pkl',
    metavar='<path>', help='path to term classifier (.pkl) file')
parser.add_argument('--entCDS', required=False, type=str,default='Nadeshiko.pkl',
    metavar='<path>', help='path to entCDS classifier (.pkl) file')
parser.add_argument('--outCDS', required=False, type=str,default='Juliet.pkl',
    metavar='<path>', help='path to outCDS classifier (.pkl) file')

arg = parser.parse_args()
seq_file = arg.fasta
path = arg.path
dict_file = arg.dict
init_file = arg.init
term_file = arg.term
entCDS_file = arg.entCDS
outCDS_file = arg.outCDS

filenames = [
    path + init_file,
    path + term_file,
    path + entCDS_file,
    path + outCDS_file,
    path + dict_file
]

if "None" in filenames:
    print("Predict: One or more filename missing, using default")
    filenames = "Default"

clfs = train.Classifiers(filenames = filenames)
results = pattern.Finder(seq_file, clfs)
a = p.cov19_Processer(results.patterns)

for ele in a.patterns:
    id = ele["id"]
    pattern = ele["patterns"]
    print(pattern)
    gff = "data/Cov19/" + id +".gff3"
    gff = gffer.Process(gff, "Cov19").gff
    ref = []
    missing = []
    diff = 0
    for id in gff:
        for gene in gff[id]:
            record = gff[id][gene]["cds"]
            if len(record) != 1:
                for ele1 in record:
                    stop = False
                    for ele2 in record:
                        if ele1[1] == ele2[0]:
                            new = [ele1[0],ele2[1]]
                            record.remove(ele1)
                            record.remove(ele2)
                            ref.append(new)
                            if len(record) != 0:
                                ref.append(record[0])
            else: ref.append(record[0])
    for record in ref:
        adjust = [record[0]-1, record[1]]
        if adjust not in pattern:
            missing.append(record)

snps = []
targets = []
with open("data/Cov19/SNPs.txt", "r") as file:
    for line in file:
        line = line.strip()
        if len(line) == 0: break
        line = line.split()
        ann = h.HGVS(line[0])
        snps.append(ann)
        targets.append(line[1])

fas_read = read.FASTA(seq_file)
for entry in fas_read:
    preds = t.Tranlate(entry.seq, a.patterns, snps)



# with open("data/Cov19/geoAnno.txt", "r") as ann:
#     geos = []
#     acc = []
#     for line in ann:
#         line = line.strip()
#         if len(line) == 0: break
#         line = line.split(':')
#         geos.append(line[0])
#         acc.append(line[1])
#
#     print(geos, acc)



# test = hgvser.HGVS('c.123A>T')
# print(test.type, test.prefix, test.info.alt, test.info.ref, test.info.pos)
