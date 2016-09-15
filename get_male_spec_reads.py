# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 12:07:18 2015

@author: lukasendler
"""

import pysam
from distutils.version import LooseVersion
#import numpy as np
#import copy
import gzip
#from scipy import stats
import sys, os
import re
from collections import defaultdict
import argparse

def read_snps(infile):
    """
    reads a cmh or sync like file and creates a dict with chrom->[ [ bps, allele ] ]  
    """
    #snp_dict=defaultdict(list)  #dictionary of chroms with positions and male specific allele
    snp_list=[]
    if re.search("\.b?gz",infile):
        inf = gzip.open(infile,'rb')
    else:
        inf = open(infile,"r")
    #load dictionary
    for line in inf:
        if re.match("\#",line):
            continue
        line=line.rstrip()
        fields=line.split()
        snp_list.append([ fields[0], int(fields[1]), fields[2] ])
    inf.close()
    return sorted(snp_list)

def get_qnames(snp,bamfh,read_set,map_qual=20,base_qual=20,v=False):
    """
    function for the new pysam synthax (pysam 0.8.0 and up), thanks to the maintainers for breaking :( 
    gets a SNP (chr, bps, allele), a bam file handle and a read_set and adds all qnames not in the set   
    """
    # create pileup and move iterator to the right position - problem due to pysam
    # if you use pysam 0.6, you will have to subtract a readlength and a bit from start (eg. start=bps-101) 
    # have to go through loop, I think
    (chrom,bps,allele)=snp
    for pile in bamfh.pileup(region=chrom,start=bps-1,end=bps):
        if pile.pos == bps-1:
            pile_col=pile.pileups
            break
    
    for pile_read in pile_col:
        # check whether read is in get_reads or col_reads
        qname=pile_read.alignment.query_name
        if qname in read_set or not pile_read.query_position:
            continue
        if ( pile_read.alignment.query_sequence[pile_read.query_position] == allele 
        and pile_read.alignment.mapping_quality >= map_qual 
        and pile_read.alignment.query_qualities[pile_read.query_position] >= base_qual ) :
            read_set.add(qname)
            #if v:
                #print "at "+ chrom + str(bps) + " allele found " + pile_read.alignment.query_sequence[pile_read.query_position] + "in read " + qname
    return read_set

################
###   MAIN   ###
################

### for testing try this bamfile
### bam_file="/Volumes/vetgrid10/Data/Vienna_2010/Realigned/pop1_earlylate_merged_RG_sanger_real.bam"
### chrom="2L"
### bps=40829


parser = argparse.ArgumentParser(description="""
gets a list of bam files and file containing specific snps and alleles and returns all reads having the specified allele
""") 
parser.add_argument("-b", dest="bam_files", help="bam files, list comma separated(eg. \"dark1bam,light1.bam\")", required=True)
parser.add_argument("-s", dest="snp_files", help="snp file", required=True)
parser.add_argument("--mq", dest="map_qual", help="mapping quality threshold", default=20)
parser.add_argument("--bq", dest="base_qual", help="base quality threshold", default=20)
parser.add_argument("-v","--verb", dest="verb", action="store_true", help="print verbose state on stderr", default=False)


args = parser.parse_args()
snp_files = vars(args)['snp_files'].split(",")
bam_files = vars(args)['bam_files'].split(",")
map_qual = int(vars(args)['map_qual'])
base_qual = int(vars(args)['base_qual'])
verb=vars(args)['verb']
snp_hash={}
for snp_file in snp_files:
    snp_list=read_snps(snp_file)
    snp_hash[snp_file]=snp_list

# outfhs dict of snpfile -> fh
outfh={}
for bam_file in bam_files:
    ## index BAM file if necessary
    if not os.path.exists(bam_file+".bai"):
        print "indexing "+ bam_file
        os.system("samtools index "+ bam_file)
    bamfh=pysam.AlignmentFile(bam_file,'rb')
    # hash of readnames -> [snpfile]
    read_hash=defaultdict(list)
    for snp_file in snp_hash.keys():
        if verb:
            print "processing snp file " + snp_file + " and " + bam_file 
        read_set = set()
        for snp in snp_hash[snp_file]:    
            read_set=get_qnames(snp,bamfh,read_set,map_qual=20,base_qual=20,v=verb)
        for qname in read_set:
            read_hash[qname].append(snp_file)
    bamfh.close()
    bamfh=pysam.AlignmentFile(bam_file,'rb')    
    if len(outfh) == 0:
        for snp_file in snp_hash.keys():
            outfh[snp_file]=pysam.Samfile(snp_file+".bam","wb",template=bamfh)
    index = 0
    for read in bamfh.fetch(until_eof=True):
        index += 1
        if verb and index%1000000 == 0:
            print "screened "+ str(index) + " reads of " + bam_file
        if read.query_name in read_hash:
            for snp_file in read_hash[read.query_name]:
                outfh[snp_file].write(read)
    bamfh.close()


