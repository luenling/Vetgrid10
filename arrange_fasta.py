import sys, os, re 
from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser(description='read a multi fasta file and either give the contig or sequence IDs, if no output file is given, or write a file with sequences arranged according to a list of IDs given. the script can also create a seperate file containing each individual sequence.') 
parser.add_argument("--in", "-i",dest="infile", help="fasta file to be read (required)", required=True)
parser.add_argument("--out","-o", dest="outfile", help="fasta file to be written (optional)", default=False)
parser.add_argument("--ids", dest="ids", help="fasta ids sequence in outputfile eg. X,2R,2L,3L,3R,4,U (def: all sequence IDs)", default=False)
parser.add_argument("--idf", dest="idf", help="fasta ids in a file, takes only first field of each line (exclusive with --ids)", default=False)
parser.add_argument("--split", dest="splt", action="store_true", help="write each chrom or contig in IDs into different file (outfile(wo. fasta/fa) + ID + fasta) (def: false)", default=False)
parser.add_argument("-v", dest="inverse", action="store_true", help="only write contigs not in ID list or file (def: false)", default=False)
args = parser.parse_args()
splt = vars(args)['splt']
infile = vars(args)['infile']
outfile= vars(args)['outfile']
inverse= vars(args)['inverse']
record_dict = SeqIO.index(infile,"fasta")
if (vars(args)['ids']):
    ids=vars(args)['ids'].split(",")
else:
    ids=list(record_dict)
if (vars(args)['idf']):
    ids=[]
    with open(vars(args)['idf']) as f:
        pat=re.compile(r'\s*#')
        pat2=re.compile(r'\s*>')
        for line in f:
            # if line starts with a comment, continue
            if (re.search(pat,line)):
                continue
            elif (re.search(pat2,line)):
                line=re.sub(pat2,"",line)
            line=line.rstrip()
            ids.append(line.split()[0])
#print ids
no_match = [x for x in ids if x not in list(record_dict) ]
match_ids = [x for x in ids if x in list(record_dict) ]
if (len(no_match) > 0):
    print >>sys.stderr, "some "+str(len(no_match))+" IDs do not match, "+str(len(ids)-len(no_match))+" match"
if inverse:
    match_ids =  [x for x in ids if not ( x in list(record_dict) ) ]
if (outfile):
    output_handle=False
    for i in match_ids:
        outfile_temp=outfile
        if splt:
            outfile_temp = re.sub("\.fa(?:sta)?","",outfile)+"_"+i+".fasta"
            if output_handle:
                output_handle.close()
        if not output_handle or output_handle.closed:
            output_handle=open(outfile_temp, "w")
        SeqIO.write(record_dict[i], output_handle, "fasta")
    output_handle.close()
else:
    print ",".join(list(record_dict))
