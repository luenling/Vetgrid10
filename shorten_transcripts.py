#!/usr/bin/python

## functions: identify overlaps
def identify_overlaps(shortexons,genes,gene_coords,exons,stranded=False):
    """
    takes list of genomic regions (exons), checks for overlapping genes, and removes the intersecting bits or creates new exons out of the non-overlapping bits
    """
    nexs=[]
    for ex in shortexons:
        start=ex.iv.start
        stop=ex.iv.end
        chrom=ex.iv.chrom
        # find overlapping genes
        olps=[ genes[x] for x in gene_coords[chrom][ np.logical_and( gene_coords[chrom][:,1] >= start, gene_coords[chrom][:,0] <= start) ][:,2] ]
        olps=[ x[1] for x in olps if not ( x[0] == ex.attr['gene_id']
                                       or ( stranded and ex.iv.strand != exons[x[1][0]].iv.strand))]
        if len(olps) == 0:
            # no overlapping genes
            nexs.append(ex)
            continue
        gas=HTSeq.GenomicArrayOfSets([chrom],stranded)
        gas[ex.iv] += "ex"
        overlaps=0
        for olp in olps:
            # check whether same gene, or different strand (if stranded)
            for idx in olp:
                if not(exons[idx].iv.start <= stop and  exons[idx].iv.end >= start):
                    continue
                gas[exons[idx].iv] += "other"
                overlaps+=1
        if overlaps == 0:
            nexs.append(ex)
            continue
        new_iv=[ x[0] for x in gas.steps() if x[1] == {'ex'} ]
        for iv in new_iv:
            new_ex=copy.deepcopy(ex)
            new_ex.iv=iv
            new_ex.attr['overlaps']=str(overlaps)
            nexs.append(new_ex)
    return nexs

## functions: shorten transcripts
def shorten_transcripts(transcripts,exons,taglen):
    """
    takes a dict of transcript ids -> [exonidxs], list of exons (GenomicRegions), and a transcript length (integer)
    returns a list of shortened exons (without any order)
    """
    shortexons=[]
    for trns in transcripts.keys():
        # get exons and sort them from 3 prime
        exs=transcripts[trns]
        # # get strand for sorting and stuff
        # strand= 1 if (exons[exs[0]].iv.strand=="+") else -1
        # sort exons so 3 prime end comes first
        exs=sorted(exs,key=lambda x: -1*exons[x].iv.end_d)
        trlen=0
        for ex in exs:
            exon=exons[ex]
            exlen=exon.iv.end-exon.iv.start
            trlen +=exlen
            if (trlen >= taglen):
                exon=copy.deepcopy(exon)
                dist=trlen-taglen
                # set length to certain value
                exon.iv.length -= dist
                exon.attr['shortened']=str(dist)
                # if strand == 1:
                #     exon.iv.start += dist
                # else:
                #     exon.iv.end -= dist
                shortexons.append(exon)
                break
            shortexons.append(exon)
    return shortexons
    

import sys,re,os
import argparse
import copy
import numpy as np
from collections import defaultdict
#import select
#import gzip
import HTSeq

parser = argparse.ArgumentParser(description="""
this script should read a gtf file, shorten all transcripts to a maximal length from the 3 prime end and also remove all sequences overlapping with other genes
using htseq to read the gtf (as genomic regions)
""")

parser.add_argument("--gtf","-g", dest="gtffile", help="gtf file", required=True)
parser.add_argument("--stranded","-s", dest="stranded", action='store_true', help="use strand info for overlaps (default=False)", default=False)
parser.add_argument("--length","-l", dest="taglen", help="length for 3 prime cutoff in bases (int, default=1000)", default=1000)
args = parser.parse_args()
gtffile = vars(args)['gtffile']
stranded = vars(args)['stranded']
taglen=int(vars(args)['taglen'])

# list of exons
exons = []
# dictionary of genes to list of exons and strand
genedict = defaultdict(list)
# dictionary: transcript id -> exons
transcripts = defaultdict(list) 
# read gtf into data structures
index = 0
for gf in HTSeq.GFF_Reader(gtffile,end_included=True):
    if gf.type != "exon":
        continue
    exons.append(gf)
    genedict[gf.attr['gene_id']].append(index)
    transcripts[gf.attr['transcript_id']].append(index)
    index+=1
#list of genes
genes=genedict.keys()
#dict of chroms with list of gene start,stops and gene idx  
gene_coords=defaultdict(list)
for i,j in enumerate(genes):
    # get exons
    gexs=[ exons[x] for x in genedict[j]]
    # get chrom
    chrom=gexs[0].iv.chrom
    start=min([ x.iv.start for x in gexs])
    stop=max([ x.iv.end for x in gexs])
    gene_coords[chrom].append([start,stop,i])
genes=[ [x,genedict[x]]  for x in genes  ]
#get starts and stops
for chrom in gene_coords.keys():
    coords=np.array(gene_coords[chrom],dtype=int)
    gene_coords[chrom]=coords[coords[:,0].argsort()]
# shorten transcripts
shortexons=shorten_transcripts(transcripts,exons,taglen)
new_exons=identify_overlaps(shortexons,genes,gene_coords,exons,stranded=False)

# write gtf lines
for exon in new_exons:
    print HTSeq.GenomicFeature.get_gff_line(exon)

