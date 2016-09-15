import sys, os, re
import numpy as np
import argparse
from collections import defaultdict
#from scipy import stats
#Author: Lukas Endler
def read_snp_af_file(infile):
    """
    reads a cmh or sync like file and creates a dict with chrom->bps->[Ref,Alleles, fileentri, avg.af, [ afs ] ]
    """
    snp_dict=defaultdict(defaultdict)  #dictionary of chroms with positions and pVals of snps
    inf = open(infile,"r")
    #load dictionary
    for line in inf:
        if re.match("\#",line):
            continue
        line=line.rstrip()
        fields=line.split()
        afs = np.array(fields[4:],dtype=float,ndmin=1) 
        snp_dict[fields[0]][int(fields[1])]=[ fields[2] , fields[3], 0 ,afs.mean(), afs ]
    inf.close()
    return snp_dict

def add_to_snp_dict(snp_dict,infile):
    """
    appends af and joins polymorphisms:
    chrom->bps->[0: Ref,1:Alleles, 2:file entry, 3:avg.af, 4:[ afs ] ]
    if a position is not mapped, it is left out of the dictionary
    """
    inf = open(infile,"r")
    for line in inf:
        if re.match("\#",line):
            continue
        line=line.rstrip()
        fields=line.split()
        bps=int(fields[1])
        afs = np.array(fields[4:],dtype=float,ndmin=1)
        if (fields[0] in snp_dict and bps in snp_dict[fields[0]]): # polymorphism in other file
            # check whether polym same, just add afs
            all_A=snp_dict[fields[0]][bps][1]
            mean_afs_A=snp_dict[fields[0]][bps][3]
            if fields[3][0:2] == all_A[0:2]:
                snp_dict[fields[0]][bps][4]=np.append(snp_dict[fields[0]][bps][4],afs)
                snp_dict[fields[0]][bps][2]=2
            if (fields[3] != all_A):
                if ( fields[3][1] == all_A[0] and fields[3][0] == all_A[1] ): # major of A is minor B (swap B)
                    afs = 1 - afs
                    snp_dict[fields[0]][bps][4]=np.append(snp_dict[fields[0]][bps][4],afs)
                    snp_dict[fields[0]][bps][2]=2
                    if snp_dict[fields[0]][bps][4].mean() < 0.5: # if mean afs < 0.5 reverse poly and change afs
                        snp_dict[fields[0]][bps][4] = 1 - snp_dict[fields[0]][bps][4]
                        snp_dict[fields[0]][bps][1] = all_A[0:2][::-1] 
                elif ( fields[3][0] == all_A[0]): # maj A equal maj B
                    if afs.sum() < snp_dict[fields[0]][bps][4].sum(): # minor B greater than A
                        snp_dict[fields[0]][bps][1] = fields[3][0:2]+all_A[1]
                    else:
                        snp_dict[fields[0]][bps][1] = all_A[0:2] + fields[3][1]
                    snp_dict[fields[0]][bps][4]=np.append(snp_dict[fields[0]][bps][4],afs)
                    snp_dict[fields[0]][bps][2]=2
                elif (fields[3][0] == all_A[1]): # major of B matches minor of A
                    if afs.sum() > snp_dict[fields[0]][bps][4].sum(): # major B greater than major A
                        snp_dict[fields[0]][bps][4] = 1 - snp_dict[fields[0]][bps][4]
                        snp_dict[fields[0]][bps][1] = all_A[0:2][::-1] + fields[3][1]
                    else:
                        snp_dict[fields[0]][bps][1] = all_A[0] + fields[3][1]
                    snp_dict[fields[0]][bps][4]=np.append(snp_dict[fields[0]][bps][4],afs)
                    snp_dict[fields[0]][bps][2]=2                    
                elif (fields[3][1] == all_A[0]): # minor of B matches major A
                    if afs.sum() > snp_dict[fields[0]][bps][4].sum(): # major B greater than major A
                        snp_dict[fields[0]][bps][4] = 1 - snp_dict[fields[0]][bps][4]
                        snp_dict[fields[0]][bps][1] = fields[3][0:2]+all_A[1]
                    else:
                        snp_dict[fields[0]][bps][1] = fields[3][0:2][::-1] + all_A[1]
                    snp_dict[fields[0]][bps][4]=np.append(snp_dict[fields[0]][bps][4],afs)
                    snp_dict[fields[0]][bps][2]=2                    
                    
                else:
                    snp_dict[fields[0]][bps][1] = all_A +"|"+fields[3]
                    snp_dict[fields[0]][bps][4]=np.append(snp_dict[fields[0]][bps][4],afs)
                    snp_dict[fields[0]][bps][2]=2
        else:
            snp_dict[fields[0]][int(bps)]=[ fields[2] , fields[3], 1 , afs.mean(), afs ] 
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
parser.add_argument("-a", dest="file_a", help="allele frequency file A", required=True)
parser.add_argument("-b", dest="file_b", help="allele frequency file B", required=True)
parser.add_argument("-p","--pops", dest="pops", help="comma seperated list of populations in each file (def: \"2,2\")", required=True)
parser.add_argument("-v","--verb", dest="verb", action="store_true", help="print verbouse state on stderr", default=False)

args = parser.parse_args()
file_a = vars(args)['file_a']
file_b= vars(args)['file_b']
pops = vars(args)['pops'].split(",")
pops=[ int(x) for x in pops ]
verb = vars(args)['verb']
if verb:
    print >> sys.stderr, "reading file A"
snp_dict=read_snp_af_file(file_a)
if verb:
    print >> sys.stderr, "reading file B"
snp_dict=add_to_snp_dict(snp_dict,file_b)

if verb:
    print >> sys.stderr, "writin joined file"
    # chrom->bps->[0: Ref,1:Alleles, 2:file entry, 3:avg.af, 4:[ afs ] ]
for chrom in sorted(snp_dict.keys()):
    for bps in sorted(snp_dict[chrom].keys()):
        if snp_dict[chrom][bps][2] < 2: # not both, needs padding with reference
            ref_pos = snp_dict[chrom][bps][1].find(snp_dict[chrom][bps][0])
            #ref_pos == 0: # all fine
            if ref_pos == 1: # invert major minor
                snp_dict[chrom][bps][4] = 1 - snp_dict[chrom][bps][4]
                snp_dict[chrom][bps][1] = snp_dict[chrom][bps][1][0:2][::-1]
            elif ref_pos < 0:
                snp_dict[chrom][bps][4] = 1 - snp_dict[chrom][bps][4]
                snp_dict[chrom][bps][1] = snp_dict[chrom][bps][0] + snp_dict[chrom][bps][1][1]
            if snp_dict[chrom][bps][2] == 1:
                snp_dict[chrom][bps][4] = np.append(np.repeat(1,pops[0]),snp_dict[chrom][bps][4])
            else:
                snp_dict[chrom][bps][4] = np.append(snp_dict[chrom][bps][4],np.repeat(1,pops[1]))
        print "{}\t{}\t{}\t{}".format(chrom, bps, "\t".join(snp_dict[chrom][bps][0:2]),"\t".join([ str(x) for x in snp_dict[chrom][bps][4]]))

