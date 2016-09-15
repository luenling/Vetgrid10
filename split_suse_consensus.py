import sys
import argparse
import numpy as np
from collections import defaultdict

def count_string(haplos,nucs):
    return np.array([ haplos.count(x) for x in nucs ])


nucs=["A","T","C","G","N"]


#Author: Lukas Endler
parser = argparse.ArgumentParser(description="""Reads Suse\'s Haplotype consensus file and writes a file with:
Chr Pos Ref Refstrain HaplotypeCounts(A:T:C:G:N) estimated_Afs_for_n_indv(A:T:C:G) afs_wo_Ns maj.allelefreq Number_of_Ns""") 
parser.add_argument("--infile","-i", dest="infile", help="consensus file to be read", required=True)
parser.add_argument("--out","-o", dest="outfile", help="output sync file", required=True)
parser.add_argument("--ind","-n", dest="indv", help="number of heterozygous individuals for allele frequency, default=length of haplotype string", default=False)

args = parser.parse_args()
infile = vars(args)['infile']
outfile = vars(args)['outfile']
indv=int(vars(args)['indv'])

inf = open(infile,"r")
out = open(outfile,"w")


for line in inf:
    entries=line.rstrip().split()
    # 0: chrom 1:pos 2:ref 3:refstrain 4: haplos 5: identificationtype 6:Nucs found
    if not indv:
        indv = len(entries[4])
    haps=count_string(entries[4],nucs)
    hap_count=":".join([str(x) for x in haps])    
    tot_af= (haps*0.5)/indv
    indv_cov = indv - entries[4].count("c")
    if indv_cov == 0:
        afs_wo_Ns = haps[:4]*0.0
    else:
        afs_wo_Ns = haps[:4]*0.5/indv_cov
    if (len(entries[3]) == 1 and entries[3] != "N"):
        tot_af[nucs.index(entries[3])] += 0.5
        # could also use the allele freq without the Ns due to c ?
        afs_wo_Ns[nucs.index(entries[3])] += 0.5
    elif (len(entries[3]) > 1): # more than one allele in reference strain
        ref_alleles = entries[3].split("/")
        entries[3]="".join(ref_alleles)
        for i in ref_alleles:
            tot_af[nucs.index(i)] += 0.5/len(ref_alleles)
            afs_wo_Ns[nucs.index(i)] += 0.5/len(ref_alleles)
    alleles = [ nucs[x] for x in np.flatnonzero(afs_wo_Ns) ] 
    alleles = "".join([ alleles[x] for x in np.argsort(afs_wo_Ns[ np.flatnonzero(afs_wo_Ns) ])[::-1] ])
    maf = "{:.3f}".format(np.max(afs_wo_Ns))
    afs = ":".join([ "{:.3f}".format(x) for x in tot_af])
    afs_wo_Ns =  ":".join([ "{:.3f}".format(x) for x in afs_wo_Ns])
    print >> out, "\t".join(entries[:4])+"\t"+hap_count+"\t"+afs+"\t"+afs_wo_Ns+"\t"+alleles+"\t"+maf+"\t"+str(haps[4])

inf.close()
out.close()

    
    
