#!/usr/bin/env python
import sys, os, glob, re
import numpy as np
from numpy import array, zeros
from collections import defaultdict
import argparse
parser = argparse.ArgumentParser(description='get groups of snps in which LL or DD are fully linked.') 
parser.add_argument("-i", dest="infile", help="file with SNPs and linkage", required=True)
parser.add_argument("-p", dest="pop", help="population to use", default=1)
parser.add_argument("-c", dest="min_cnt", help="minimal counts", default="100")
parser.add_argument("-m", dest="max_perc", help="maximal percent of LD linked", default="0.025")
parser.add_argument("-r", dest="min_r2", help="minimal r2 value", default="1.0")
args = parser.parse_args()
infile = vars(args)['infile']
pop = int(vars(args)['pop'])
min_cnt = int(vars(args)['min_cnt'])
min_r2 = float(vars(args)['min_r2'])
max_perc = float(vars(args)['max_perc'])
counts=(pop-1)*5+3
perc=(pop-1)*5+6
r2=(pop-1)*5+7
inf = open(infile,"r")
snp_dict=defaultdict(list)
for line in inf:
    entries = line.rstrip().split()
    # if entries[0] == "X" and entries[1] == "9121094":
    #    print "\t".join(entries[0:3])+"\t"+entries[perc]+"\t"+entries[counts]+"\t"+str(max_perc)
    #    print str(( entries[counts] == "nan" or int(entries[counts]) < min_cnt ) or ( float(entries[perc]) > max_perc and float(entries[r2]) < min_r2))
    if ( entries[counts] == "nan" or int(entries[counts]) < min_cnt ) or ( float(entries[perc]) > max_perc and float(entries[r2]) < min_r2) :
        continue
    # print "\t".join(entries[0:3])
    snp_dict[entries[0]].append(set([ int(x) for x in entries[1:3]]))
inf.close()
for chrom in snp_dict.keys():
    for i,s_a in enumerate(snp_dict[chrom]):
        for j in range(0,i):
            if snp_dict[chrom][j].intersection(s_a):
                snp_dict[chrom][j]=snp_dict[chrom][j].union(s_a)
                snp_dict[chrom][i]=set()
                continue
out_hap=open(infile+"_mc_"+str(min_cnt)+"_maxP_"+str(max_perc)+"_mr2_"+str(min_r2)+".hap_stretch","w")
out_sync=open(infile+"_mc_"+str(min_cnt)+"_maxP"+str(max_perc)+"_mr2_"+str(min_r2)+".hap_sync","w")

for chrom in snp_dict:
    hapblck=0
    for s_a in snp_dict[chrom]:
        if len(s_a) == 0:
            continue
        bpos="\t".join([str(x) for x in sorted(s_a)])
        print >>out_hap,chrom+"\t"+bpos
        bpos=[str(x) for x in sorted(s_a)]
        for bps in bpos:
            print >>out_sync,chrom+"\t"+bps+"\t"+chrom+"_H_"+str(hapblck)
        hapblck+=1
out_hap.close()
out_sync.close()





    
    
