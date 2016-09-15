import sys, os, re 
from Bio import SeqIO
from collections import defaultdict
import argparse
parser = argparse.ArgumentParser(description='read a multi fasta file and only output sequences over a threshold length.') 
parser.add_argument("--in", "-i",dest="infile", help="fasta file to be read (required)", required=True)
parser.add_argument("--getorf", "-g",dest="getorf",  action="store_true" ,help="getorf created fasta (default:False)", default=False)

args = parser.parse_args()
infile = vars(args)['infile']
getorf = vars(args)['getorf']
outfile= re.split("\.\w+$",infile)[0]+"_longest_orf.fna"
handle = open(infile, "rU")
records = list(SeqIO.parse(handle, "fasta"))
handle.close()
seq_dict = {}
re_id=re.compile("(\w+)[|_]\d+")
if getorf:
    re_len=re.compile("\[(\d+)\s?-\s?(\d+)\]")
else:
    re_len=re.compile("len(?:gth)?[=:](\d+)")
# conid=re.search(re_id,line).group(1)
# length=int(re.search(re_len,records[0].description).group(1))
seq_dict=defaultdict(list)
contigs=[]
for i in range(0,len(records)):
    conid=re.search(re_id,records[i].id).group(1)
    if getorf:
        length=abs( int(re.search(re_len,records[i].description).group(1))-int(re.search(re_len,records[i].description).group(2))) 
    else:
        length=int(re.search(re_len,records[i].description).group(1))
    seq_dict[conid].append([i,length])
for conid in seq_dict.keys():
    seq_dict[conid].sort(key=lambda x: x[1],reverse=True)
contigs=seq_dict.keys()
contigs.sort(key=lambda x: seq_dict[x][0][0])
output_handle = open(outfile, "w")
for contig in contigs:
    idx=seq_dict[contig][0][0]
    records[idx].id=contig
    records[idx].name=None
    records[idx].description=""
    SeqIO.write(records[idx], output_handle, "fasta")
output_handle.close()


