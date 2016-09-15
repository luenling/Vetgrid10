import sys, os, re 
from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser(description='read a multi fasta file and only output sequences over a threshold length.') 
parser.add_argument("--in", "-i",dest="infile", help="fasta file to be read (required)", required=True)
parser.add_argument("--len","-l", dest="length", help="minimal length of sequences", default=True)
args = parser.parse_args()
infile = vars(args)['infile']
length=int(vars(args)['length'])
for sequence in SeqIO.parse(infile,"fasta"):
    if len(sequence) >= length:
        print sequence.format("fasta")
        #SeqIO.write(sequence, "fasta")
