import sys, os, re
import numpy as np
import argparse
from collections import defaultdict
#from scipy import stats
#Author: Lukas Endler
def read_snp_af_file(infile):
    """
    reads a cmh or sync like file and creates a dict with chrom->bps->[Alleles, rest of line]
    """
    snp_dict=defaultdict(defaultdict)  #dictionary of chroms with positions and pVals of snps
    inf = open(infile,"r")
    #load dictionary
    for line in inf:
        if re.match("\#",line):
            continue
        line=line.rstrip()
        fields=line.split()
        snp_dict[fields[0]][int(fields[1])]=[ fields[2] ,"\t".join(fields[3:]) ]
    inf.close()
    return snp_dict

def read_map_dict(snp_dict,mapfile):
    """
    reads a mapping file from mel to sim and creates a dictionary with entries for all positions of snp_dict a with mapped position:
    chrom->bps(sim)-> bps(mel), sim_reference
    if a position is not mapped, it is left out of the dictionary
    snp_dict: chr->bps->[ allele, rest, [ sim_bps, ref_allele, sim_allele, sim_rest ]  ]
    """
    map_dict=defaultdict(defaultdict)  #dictionary of chroms with positions and pVals of snps
    inf = open(mapfile,"r")
    for line in inf:
        if re.match("\#",line):
            continue
        line=line.rstrip()
        fields=line.split()
        chrom=fields[0]
        if (chrom in snp_dict) and (int(fields[1]) in snp_dict[chrom]):
            if re.match("\d+",fields[3]):
                map_dict[chrom][int(fields[3])]= int(fields[1])
                snp_dict[chrom][int(fields[1])].append( [ int(fields[3]), fields[4] ] )
    inf.close()
    return map_dict

def add_sim_polym(map_dict,snp_dict,asim):
    """
    adds simulans polymorphisms to snp_dict where entries exist in map_dict:
    chrom->bps-> [ [ bps(mel), sim_reference, [ alleles , rest of line ] ]
    """
    inf = open(asim,"r")
    for line in inf:
        if re.match("\#",line):
            continue
        line=line.rstrip()
        fields=line.split()
        chrom=fields[0].split("_")[-1]
        if (chrom in map_dict) and (int(fields[1]) in map_dict[chrom]):
            bps=map_dict[chrom][int(fields[1])]
            snp_dict[chrom][bps][2] = [ int(fields[1]), fields[2], fields[3][0:2], "\t".join(fields[4:])  ]
    inf.close()
    return snp_dict
    
    

parser = argparse.ArgumentParser(description="""
reads an allele frequency file in Sync Format (CHR BPS Maj/Minor Alleles), a file mapping mapping of Dmel positions to Simulans as provided by Ram\'s Mauve scripts, and an allele frequency file for the simulans strains.
The script gives on stdout a list with coordinates for mel and sim, an allele matching code, and the D mel. & D sim. alleles, followed by the major allele frequencies. The code encrypts whether the melanogaster and simulans alleles correspond:
no allele overlap: 0
and combinations of the major (1), or minor (2) allele of simulans matching the major and minor of mel (position):
eg. major mel = major sim, minor mel != minor sim : 10
 minor mel = minor sim, major mel != major sim : 2
 minor mel = major sim: 1
  major mel = minor sim: 20
""") 
parser.add_argument("-a","--amel", dest="affile", help="allele frequency file for melanogaster SNPs", required=True)
parser.add_argument("-s","--asim", dest="asim", help="allele frequency file for simulans", required=True)
parser.add_argument("-m","--map", dest="mapfile", help="comma seperated list of files with D. mel - D. sim. genome position mappings", required=True)
parser.add_argument("-v","--verb", dest="verb", action="store_true", help="print verbouse state on stderr", default=False)

args = parser.parse_args()
affile = vars(args)['affile']
asim = vars(args)['asim']
mapfile = vars(args)['mapfile']
verb = vars(args)['verb']
if verb:
    print >> sys.stderr, "reading mel af file"
snp_dict=read_snp_af_file(affile)
if verb:
    print >> sys.stderr, "reading mel sim mapping file"
map_dict=read_map_dict(snp_dict,mapfile)
if verb:
    print >> sys.stderr, "adding sim polymorphism data to snp dict"
snp_dict=add_sim_polym(map_dict,snp_dict,asim)
if verb:
    print >> sys.stderr, "writing output file"

count=0
for chrom in sorted(snp_dict.keys()):
    for bps in sorted(snp_dict[chrom].keys()):
        count += 1
        if verb and count%500000 == 0:
            print >> sys.stderr, "read "+str(count)+" SNPs"
        if len(snp_dict[chrom][bps]) < 3: # not mapped
            continue
        sim_rest=False
        ( sim_bps, sim_all) = snp_dict[chrom][bps][2][0:2]
        if len(snp_dict[chrom][bps][2]) > 2:
            ( sim_all, sim_rest) = snp_dict[chrom][bps][2][2:4]
        mel_all = snp_dict[chrom][bps][0]
        mel_rest = snp_dict[chrom][bps][1]
        if len(mel_all) > 2:
            # only take biallelic SNPs / get rid of comparison in joined af files
            mel_all = mel_all[0:2]
        if len(sim_all) > 2:
            # only take biallelic SNPs / get rid of comparison in joined af files
            if sim_all.find('|') > 0:
                sim_all = sim_a[0] + sim_all[sim_all.find('|')+1]
            else:
                sim_all = sim_all[0:2]
        all_cons = [] # allele conservation (1: major, 0: none, 2: minor is ancestral)
        for al in mel_all:
            all_cons.append(str(sim_all.find(al)+1))
        if sim_rest:
            mel_rest=mel_rest+"\t"+sim_rest
        print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrom, bps, sim_bps, mel_all, sim_all,len(sim_all),"".join(all_cons),mel_rest)
