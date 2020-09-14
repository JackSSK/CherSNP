import observer as obs
import argparse
import hgvser
import trainer as train

parser = argparse.ArgumentParser(description=
    'Train Classifiers Based on Sequence and Associated GFF Annotation')
# Sequence and Annotation's path
parser.add_argument('--fasta', required=True, type=str,
    metavar='<path>', help='path to fasta file')
parser.add_argument('--gff', required=True, type=str,
    metavar='<path>', help='path to GFF3 file')
# Clf files and Dict file's path
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

# Need more parser for filenames

arg = parser.parse_args()
seq_file = arg.fasta
gff_file = arg.gff

train.Trainer(seq_file, gff_file)
