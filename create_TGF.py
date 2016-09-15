import sys, os
from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser(description='read a multi fasta file and create a dummy TGF file for Syzygy.') 
parser.add_argument("--in", dest="infile", help="fasta file to be read (required)", required=True)
parser.add_argument("--out", dest="outfile", help="tgf file to be written (required)", required=True)
parser.add_argument("--ids", dest="ids", help="fasta ids sequence in outputfile eg. X,2R,2L,3L,3R,4,U (def: all sequence IDs)", default=False)
args = parser.parse_args()
infile = vars(args)['infile']
outfile= vars(args)['outfile']
record_dict = SeqIO.index(infile,"fasta")
if (vars(args)['ids']):
    ids=vars(args)['ids'].split(",")
else:
    ids=list(record_dict)
no_match = [x for x in ids if x not in list(record_dict) ]
if (len(no_match) > 0):
    sys.exit("some IDs do not match: "+str(no_match))
output_handle = open(outfile, "w")
print >>output_handle, "\#FEATURE_NAME\tCHR\tSTART_POSITION\tEND_POSITION\tLENGTH\tGENOME_BUILD"
for i in record_dict.keys():
    length=len(record_dict[i].seq)
    print >>output_handle, "CHR_{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(i,i,1,length,length,999)
output_handle.close()

