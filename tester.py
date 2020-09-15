import observer as obs
import argparse
import hgvser
import trainer as train
import processer as process

parser = argparse.ArgumentParser(description='Just for test now')
parser.add_argument('--fasta', required=True, type=str,
    metavar='<path>', help='path to fasta file for query sequences')
parser.add_argument('--dict', required=False, type=str,default='None',
    metavar='<path>', help='path to dict(.js) file')
parser.add_argument('--init', required=False, type=str,default='None',
    metavar='<path>', help='path to init classifier (.pkl) file')
parser.add_argument('--term', required=False, type=str,default='None',
    metavar='<path>', help='path to term classifier (.pkl) file')
parser.add_argument('--entCDS', required=False, type=str,default='None',
    metavar='<path>', help='path to entCDS classifier (.pkl) file')
parser.add_argument('--outCDS', required=False, type=str,default='None',
    metavar='<path>', help='path to outCDS classifier (.pkl) file')

arg = parser.parse_args()
seq_file = arg.fasta
dict_file = arg.dict
init_file = arg.init
term_file = arg.term
entCDS_file = arg.entCDS
outCDS_file = arg.outCDS

filenames = [
    init_file,
    term_file,
    entCDS_file,
    outCDS_file,
    dict_file
]
if "None" in filenames:
    print("Predict: One or more filename missing, using default")
    filenames = "Default"

clfs = train.Classifiers(filenames = filenames)
process.Processor(seq_file, clfs)


# test = hgvser.HGVS('c.123A>T')
# print(test.type, test.prefix, test.info.alt, test.info.ref, test.info.pos)
