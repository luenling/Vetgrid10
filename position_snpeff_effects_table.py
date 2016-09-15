from collections import defaultdict
import sys, os, re
import argparse
parser = argparse.ArgumentParser(description="reads vcf file created by snpeff and adds 0 or 1 for each effect at each position, zero if not found") 

parser.add_argument("--in","-i", dest="infile", help="vcf file", required=True)
args = parser.parse_args()
infile = vars(args)['infile']
field_names="CODON_CHANGE,DOWNSTREAM,EXON,INTERGENIC,INTRON,NON_SYNONYMOUS_CODING,NON_SYNONYMOUS_START,START_GAINED,START_LOST,STOP_GAINED,STOP_LOST,SYNONYMOUS_CODING,SYNONYMOUS_STOP,UPSTREAM,UTR_3_PRIME,UTR_5_PRIME".split(",")
effects=dict.fromkeys(field_names,0)
inf=open(infile,"r")
print "CHR\tBPS\t","\t".join(field_names)
for line in inf:
    if re.match('^\#', line):
        continue
    effects=dict.fromkeys(field_names,0)
    line = line.rstrip()
    entries= line.split()
    effect_list=re.split("=|\(.*?\),?",entries[-1])[1:]
    for i in effect_list:
        if i in effects.keys():
            effects[i] = 1
    effect_map="\t".join([ str(effects[x]) for x in field_names])
    print "\t".join(entries[0:2]),"\t",effect_map

    
    
