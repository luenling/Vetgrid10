import sys, os, re 
from Bio import SeqIO
import argparse
from Bio import SearchIO
from collections import defaultdict

parser = argparse.ArgumentParser(description='read a multi fasta file to get the sequence IDs, and remove all other query entries from a blast xml file.') 
parser.add_argument("-f",dest="ffile", help="fasta file to be read (required)", required=True)
parser.add_argument("-b", dest="bfile", help="blast xml file",  required=True)
parser.add_argument("-o", dest="outfile", help="new, cleaned blast xml file",  required=True)

args = parser.parse_args()
ffile = vars(args)['ffile']
bfile = vars(args)['bfile']
outfile = vars(args)['outfile']
ids = set(SeqIO.index(ffile,"fasta"))

blast_handles = SearchIO.parse(bfile, 'blast-xml')
queries=0
matched=0
qresults=[]

for blast_handle in blast_handles:
    quid=blast_handle.id
    queries += 1
    if quid in ids:
        matched += 1
        qresults.append(blast_handle)    
blast_handles.close()
print >> sys.stderr, "matched queries: "+str(matched)+" not matched found in blast xml:"+str(queries - matched)+" left over in fasta:"+str(len(ids)-matched)

SearchIO.write(qresults, outfile, 'blast-xml')
