# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 14:02:54 2016

@author: lukasendler

This file should just read a haploid vcf file and transform all Ns into . , change the GT of other alleles and clean them up
"""

import sys, re
import os 
#import numpy as np
#from scipy import stats
import argparse
import gzip


parser = argparse.ArgumentParser(description="""Go through a VCF file and add fisher strand bias (FS), ratio of fw to rev (SB), mean tail distance for the ref and alt allele (RTD and ATD),tail distance bias (TDB), and the absolute value of the RPB calculated by samtools (t test that alt allele more likely to lie in the last or first 11 bases than the ref).  
""")

parser.add_argument("--in","-i", dest="vcffile", help="vcf-file, tries stdin if not set use \"STDIN\" or nothing for piping; default \"False\"", default=False)

args = parser.parse_args()
vcf_file = vars(args)['vcffile']
# open vcf file
if not vcf_file or vcf_file == "STDIN":
    if not sys.stdin.isatty():
        inf = sys.stdin
    else:
        sys.exit("No vcf file or stdinput given")
else:
    if re.search("\.b?gz",vcf_file):
        inf = gzip.open(vcf_file,'rb')
    else:
        inf = open(vcf_file,"r")


for line in inf:
    line = line.rstrip()
    if (re.match("^\s*\#+",line)): # entry is comment/header
        print line
        continue
    entries = line.split("\t")
    if not ("N" in entries[4]):
        print line
        continue
    # there is an N in the alleles list but nothing else, or not a pure snp omit entry
    if not ("," in entries[4]) or (len(entries[3]) > 1):
        continue
    alleles = entries[4].split(",")
    maxal= 0
    for i in alleles:
        maxal = max(len(i),maxal)
    if maxal > 1:
        continue
    #print line
    try: 
        nallele = alleles.index("N")+ 1
    except: 
        sys.exit(line)
    alleles.remove("N")
    entries[4] = ",".join(alleles)
    # change all gentypes (all haplotype and < 10)
    for i in range(9,len(entries)) :
        try:
            gt = int(entries[i])
        except:
            continue
        if (gt < nallele):
            continue
        elif (gt == nallele):
            gt = "."
        else :
            gt = gt - 1
        entries[i] = str(gt) + entries[i][1:]
    #print line
    print "\t".join(entries)

    
          
 
    
