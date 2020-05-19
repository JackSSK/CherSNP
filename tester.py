import trainer as train
import argparse

parser = argparse.ArgumentParser(description='Just for test now')
parser.add_argument('--fasta', required=True, type=str,
    metavar='<path>', help='path to fasta file')
parser.add_argument('--gff', required=True, type=str,
    metavar='<path>', help='path to GFF3 file')
arg = parser.parse_args()
seq_file = arg.fasta
gff_file = arg.gff

train.Trainer(seq_file, gff_file)
