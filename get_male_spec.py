# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 10:30:00 2015

@author: lukasendler
"""

import pysam
from distutils.version import LooseVersion
import numpy as np
import copy
import gzip
#from scipy import stats
import sys
import re
from collections import defaultdict
import argparse

def read_snps(infile):
    """
    reads a cmh or sync like file and creates a dict with chrom->[ [ bps, allele ] ]  
    """
    snp_dict=defaultdict(list)  #dictionary of chroms with positions and male specific allele
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
        snp_dict[fields[0]].extend([ int(fields[1]), fields[2] ])
    inf.close()
    return snp_dict


def get_reads_and_qnames(snp,bamfh,get_reads,col_reads,read_list,map_qual=20,base_qual=20):
    """
    function for the new pysam synthax (pysam 0.8.0 and up), thanks to the maintainers for breaking :( 
    gets a SNP position (chrom & bps), SNP indx, a bam file handle and a read_dict ( defaultdict(list)), adds the reads with snp index and nuc to the read_dict and returns an array with the ACTG counts
    read_dict[readname]->[ [snp indx,nucleotide] ]    
    """
    # create pileup and move iterator to the right position - problem due to pysam
    # if you use pysam 0.6, you will have to subtract a readlength and a bit from start (eg. start=bps-101) 
    # have to go through loop, I think
    (chrom,bps,allele)=snp
    for pile in bamfh.pileup(region=chrom,start=bps-1,end=bps):
        if pile.pos == bps-1:
            pile_col=pile.pileups
            break
    #b=[x.alignment.seq[x.qpos] for x in pile_col ]
    #counts=[ b.count(x) for x in ["A","C","G","T"]]
    #read_dict=defaultdict(list)
    snp_nucs=[]
    for pile_read in pile_col:
        # check whether read is in get_reads or col_reads
        qname=pile_read.alignment.query_name
        if qname in get_reads and pile_read.alignment.is_read1 == get_reads[qname][2]:
            # add the read to read_list, remove from get_reads and add to col_reads
            add_read(pile_read.alignment,read_list,col_reads,get_reads)
            continue
        elif qname in col_reads and ( (pile_read.alignment.is_read1 and col_reads[qname][3] !=2 ) or (pile_read.alignment.is_read2 and col_reads[qname][3] != 1 )):
            # continue, has been collected already            
            continue
        # check whether matched, mapping qual and base qual alright
        elif  pile_read.alignment.query_sequence[pile_read.query_position] != allele or (pile_read.alignment.mapping_quality < map_qual) or pile_read.alignment.query_qualities[pile_read.query_position] < base_qual :
            continue
        # should be the right read now
        #add to col_reads, read_list, get mate postion for get_reads
        read_list[pile_read.alignment.query_name].append([snp_indx,pile_read.alignment.query_sequence[pile_read.query_position]])
        snp_nucs.append(pile_read.alignment.query_sequence[pile_read.query_position])
    return np.array([ snp_nucs.count(x) for x in ["A","C","G","T"]],dtype=int)

def add_read(align,read_list,col_reads,get_reads):
    # add read to readlist and col_reads
    read_list.append(align)
    readnum=1 if align.is_read1 else 2
    qname=align.query_name
    # col_reads value either 0 if new, else 1 or 2
    col_reads[qname] += readnum
    # check get_reads and remove entry if necessary
    if qname in get_reads:
        if get_read[qname][2] == align.is_read1:
            del get_read[qname]
    # if not both reads already, add get_reads entry
    if col_reads[qname] !=3:
        chrom=    
        # get mate reference id and position
        if not align.mate_is_unmapped:
        
    
def create_search_list(get_reads):
    # create a sorted list of positions from get_reads
    positions=[ get_reads[q][0:2] for q in get_reads.keys()]
    return sorted(positions)
    


#reads collected: read name -> 1:read1,2:read2,3:both reads 
col_reads=defaultdict(int)
#reads to gather:  qname->[chr,pos,isread1]
get_reads=defaultdict(list)





