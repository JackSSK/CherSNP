import observer as obs
import argparse
import hgvser
import gen_tools as tool

parser = argparse.ArgumentParser(description='Just for test now')
parser.add_argument('--fasta', required=True, type=str,
    metavar='<path>', help='path to fasta file')
parser.add_argument('--gff', required=True, type=str,
    metavar='<path>', help='path to GFF3 file')
arg = parser.parse_args()
seq_file = arg.fasta
gff_file = arg.gff

obs.Observer(seq_file, gff_file)


# test = hgvser.HGVS('c.123A>T')
# print(test.type, test.prefix, test.info.alt, test.info.ref, test.info.pos)
