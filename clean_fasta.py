import sys, os, re
import argparse
import gzip
parser = argparse.ArgumentParser(description="""
takes a file with contig IDs possibly in some blast formats and a multi fasta file and gives only contigs not matching to E coli or not matching anything with a threshold. For the minimal alignment length and minimal covered query sequence cutoffs, the shorter criterium counts for each contig.
eg.: blastn -db nt -query /Volumes/Temp/Lukas/Documents/Projects/Ivana_Bilic/Hmeleagridis/Ecoli_cleanup/HP_full.fasta -out HP_full_blast.txt -outfmt \'7 qseqid sseqid evalue qstart qend length sstart send slen pident stitle qcovs qlen\' -num_threads 15 -evalue 1e-5
"""
) 
parser.add_argument("--fasta","-f", dest="fastaf", help="Fasta file name", required=True)
parser.add_argument("--out","-o", dest="outf", help="cleaned Fasta file (default: False)", default=False)
parser.add_argument("--blastfile","-b", dest="blastf", help="blast file",required=True)
parser.add_argument("--eval","-e",  action="store_true", dest="eval", help="only use evalue cutoff",default=False)
parser.add_argument("--maxe", dest="maxe", help="max evalue cutoff (def: 1e-6)",default=1e-6)
parser.add_argument("--mincov", dest="mincov", help="minimal query sequence coverage (def: 0.1)",default=0.1)
parser.add_argument("--ident", dest="ident", help="minimal identity cutoff (def: 80)",default=80)
parser.add_argument("--minal", dest="minal", help="minimal aligned sequence length cutoff (def: 100)",default=100)
parser.add_argument("--any",  action="store_true", dest="any", help="remove all sequences having any match above criteria",default=False)
parser.add_argument("--good",  action="store_true", dest="good_conts", help="print list of good contigs and not bad ones (default: false)",default=False)

args = parser.parse_args()
fastaf = vars(args)['fastaf']
blastf = vars(args)['blastf']
outf = vars(args)['outf']
m_eval = vars(args)['eval']
good_conts = vars(args)['good_conts']
m_any = vars(args)['any']
# cutoff values
max_e=float(vars(args)['maxe'])
min_coverage=float(vars(args)['mincov'])
min_al=float(vars(args)['minal'])
ident=float(vars(args)['ident'])
if m_eval:
    min_coverage=0.0
    min_al=0
    ident=0.0
contig_pat=re.compile("\#\s+Query:\s+(\w+)") # using the qlen  else: \s+length\=(\d+)\s+")
tax_pat = re.compile("escherichia\s+coli",re.I)
# create contig list to ditch
contigs={}
bad_cont=set()
res_fields={}
bfh=open(blastf,"r")
for line in bfh:
    if re.match("\#",line):
        a=re.search(contig_pat,line)
        if a:
            contigs[a.group(1)]=1# using the qlen  float(a.group(2))
        elif re.match("\#\s+Fields:",line):
            line=line.rstrip()
            fields=re.split("[:,]",line)[1:]
            res_fields={ fields[x].lstrip() : x for x in range(0,len(fields))}
        continue
    if re.search(tax_pat,line) or m_any:
        # is e coli hit
        fields=line.split("\t")
        cont_id=fields[res_fields["query id"]]
        assert (cont_id in contigs.keys()), cont_id + " not dictionary of contigs"
        evalue=float(fields[res_fields["evalue"]])
        al_len=float(fields[res_fields["alignment length"]])
        perc_id=float(fields[res_fields['% identity']])
        q_len=float(fields[res_fields['query length']])
        if  evalue <= max_e and ( al_len/q_len >= min_coverage or al_len >= min_al ) and perc_id >= ident:
            bad_cont.add(cont_id)

bfh.close()
good_cont=set.difference(set(contigs.keys()),bad_cont)
if good_conts:
    for a in sorted(good_cont):
        print a
else:
    for a in sorted(bad_cont):
        print a
# go through fasta file and only output contigs that are not to be ditched
if outf == False:
    sys.exit()
print_entry=False
ofh=open(outf,"w")
ffh=open(fastaf,"r")
for line in ffh:
    if re.search(">",line):
        contig_name=line.split()[0][1:]
        if contig_name in bad_cont:
            print_entry=False
        else:
            print_entry=True
    if print_entry:
        print >>ofh, line.rstrip()
ffh.close()
ofh.close()
