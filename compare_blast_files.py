def check_hit(algn, qlen, evalue,ident,hitlen,hitfr):
    """
     gets hsp object and query length, min evalue,ident,hitlen,hitfr and returns True if good hsp, else False
    """
    #idcs = False
    if not (algn.seq_len >= hitlen or float(algn.seq_len)/qlen >= hitfr):
        # too short
        return False 
    for hsp in algn.hsps:
        hid=float(hsp.ident_num)/hsp.aln_span
        if (hsp.evalue <= evalue and hid >= ident ):
            #and ( hsp.aln_span >= hitlen or float(hsp.aln_span)/qlen >= hitfr  ) ):
            return True
    # nothing good matching found :(
    return False

def parse_bfile(bfile,gi_contigs,gi_defs, evalue,ident,hitlen,hitfr,max_q=False,verb=False):
    """
    get a blast xml-file name, go through it and fill two dicts, one mapping hot gis to a set of query ids, the other to definitions
    """
    #bfh = open(bfile)
    #blast_handles=NCBIXML.parse(bfh)
    blast_handles = SearchIO.parse(bfile, 'blast-xml')
    queries=0
    for blast_handle in blast_handles:
        if (max_q != False and max_q < queries):
            print "reached "+str(max_q)+" entries in "+bfile
            break 
        quid=blast_handle.id
        qlen=blast_handle.seq_len
        queries += 1
        if ( verb and queries%100 == 0):
            print queries
        for algn in blast_handle.hits:
            # get algn accession and length
            if ( check_hit(algn,qlen,evalue,ident,hitlen,hitfr) == True):
                gi_contigs[algn.id].add(quid)
                gi_defs[algn.id]=algn.description
            else:
                continue
    blast_handles.close()
    return True

def get_contig_dict(gi_contig,gis=False): 
    """
    turn a dict of gis to contigs into a contigs to gi dict
    """
    cont_gi=defaultdict(set)
    if not gis:
        gis=gi_contig.keys()
    for gi in gis:
        for contig in gi_contig[gi]:
            cont_gi[contig].add(gi)
    return cont_gi

# this script should get all positions aligned to a certain position of a subject sequence
#from Bio.Blast import NCBIXML
from Bio import SearchIO
import sys
import os 
from collections import defaultdict
import argparse
#from string import maketrans

parser = argparse.ArgumentParser(description="""
 this script should compare two blast xml files and get all entries with unique gi number or other identifier hits above a certain not really finished up to now ;)
needs Biopython > 1.6
""") 

parser.add_argument("--hpex", dest="hpex", help="file with contigs only in hp", default=False)
parser.add_argument("--lpex", dest="lpex", help="file with contigs only in lp", default=False)
parser.add_argument("--blast1","-b1", dest="bfile1", help="xml file with blast results 1", required=True)
parser.add_argument("--blast2","-b2", dest="bfile2", help="xml file with blast results 2", required=True)
parser.add_argument("-e", dest="evalue", help="max evalue threshold; default: 1e-6", default=1e-6)
parser.add_argument("-i", dest="ident", help="min identity; default: 0.2", default=0.2)
parser.add_argument("-l", dest="hitlen", help="min hit length; default: 60", default=60)
parser.add_argument("-f", dest="hitfr", help="min hit fraction; default: 0.1", default=0.1)
parser.add_argument("--out1", dest="out1", help="basename of file 1 output files", default="HP")
parser.add_argument("--out2", dest="out2", help="basename of file 2 output files", default="LP")
parser.add_argument("--max_q", dest="max_q", help="max number of query entries to consider, default: False", default=False)
parser.add_argument("--verb", dest="verb", action="store_true",help="max number of query entries to consider, default: False", default=False)


args = parser.parse_args()
bfile1 = vars(args)['bfile1']
bfile2 = vars(args)['bfile2']
hpex = vars(args)['hpex']
lpex = vars(args)['lpex']
out1 = vars(args)['out1']
out2 = vars(args)['out2']
evalue = float(vars(args)['evalue'])
ident = float(vars(args)['ident'])
hitlen = int(vars(args)['hitlen'])
hitfr = float(vars(args)['hitfr'])
verb = vars(args)['verb']
max_q = vars(args)['max_q']
if max_q != False:
    max_q=int(max_q)
# create dicts for all gi numbers of good hits holding a set of contig ids and definitions
gi1_contigs=defaultdict(set)
gi1_defs={}
parse_bfile(bfile1,gi1_contigs,gi1_defs, evalue,ident,hitlen,hitfr,max_q,verb)
gi2_contigs=defaultdict(set)
gi2_defs={}
parse_bfile(bfile2,gi2_contigs,gi2_defs, evalue,ident,hitlen,hitfr,max_q,verb)

gi1=set(gi1_contigs.keys())
gi2=set(gi2_contigs.keys())
# common gis
gic=gi1.intersection(gi2)

# only bfile1
gi1x=gi1.difference(gi2)
cont1x=get_contig_dict(gi1_contigs,gi1x)
# only bfile2
gi2x=gi2.difference(gi1)
cont2x=get_contig_dict(gi2_contigs,gi2x)

# print gi1x
outf=open(out1+"_gis.txt",'w')
for gi in sorted(gi1x):
    print >> outf, "{}\t{}\t{}".format(gi,gi1_defs[gi]," ".join(sorted(gi1_contigs[gi])))
outf.close()
outf=open(out1+"_conts.txt",'w')
for cont in sorted(cont1x.keys()):
    print >> outf, "{}\t{}".format(cont," ".join(sorted(cont1x[cont])))
outf.close()


# print gi2x
outf=open(out2+"_gis.txt",'w')
for gi in sorted(gi2x):
    print >> outf, "{}\t{}\t{}".format(gi,gi2_defs[gi]," ".join(sorted(gi2_contigs[gi])))
outf.close()
outf=open(out2+"_conts.txt","w")
for cont in sorted(cont2x.keys()):
    print >> outf, "{}\t{}".format(cont," ".join(sorted(cont2x[cont])))
outf.close()


if (hpex and lpex):
    cont1=get_contig_dict(gi1_contigs)
    cont2=get_contig_dict(gi2_contigs)
    gi1x=set()
    gi2x=set()
    with open(hpex) as f:
        hpex = set([x.rstrip() for x in f])  
    with open(lpex) as f:
        lpex = set([x.rstrip() for x in f])
    # write the contigs exclusive to condition 1 (HP)
    outf=open(out1+"exclusive_conts.txt",'w')
    for cont in sorted(hpex):
        # check whether gis not in other gi
        if cont1[cont].intersection(gi2):
            # if overlap between gis of contig and other blast gis
            continue
        gi1x.update(cont1[cont])
        print >> outf, "{}\t{}".format(cont," ".join(sorted(cont1[cont])))
    outf.close()
    # write the contigs exclusive to condition 2 (LP)
    outf=open(out2+"exclusive_conts.txt",'w')
    for cont in sorted(lpex):
        # check whether gis not in other gi
        if cont2[cont].intersection(gi1):
            # if overlap between gis of contig and other blast gis
            continue
        gi2x.update(cont2[cont])
        print >> outf, "{}\t{}".format(cont," ".join(sorted(cont2[cont])))
    outf.close()
    # write the gis and defs in these contigs
    outf=open(out1+"exclusive_gis.txt",'w')
    for gi in sorted(gi1x):
        print >> outf, "{}\t{}\t{}".format(gi,gi1_defs[gi]," ".join(sorted(gi1_contigs[gi])))
    outf.close()
    outf=open(out2+"exclusive_gis.txt",'w')
    for gi in sorted(gi2x):
        print >> outf, "{}\t{}\t{}".format(gi,gi2_defs[gi]," ".join(sorted(gi2_contigs[gi])))
    outf.close()
    
         
    

