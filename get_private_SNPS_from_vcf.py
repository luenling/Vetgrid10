def sample_allele(pop,fields,gq):
    # gets array with population indeces and fields (starting from idx 8 on)
    # returns common allele, -1 if none found or if not the same in population or heterozygous
    all=-1
    for idx in pop:
        all_gq=fields[idx].split(":")
        if len(all_gq) < 3 or int(all_gq[2]) < gq:
            continue
        alleles=re.split("[|/]",all_gq[0]) 
        if len(alleles) != 2:
            continue
        elif alleles[0] !=  alleles[1] or not (alleles[0] in "01") :
            all=-1
            break
        if all == -1:
            all = int(alleles[0])
        elif all != int(alleles[0]):
            all=-1
            break
    return all

import sys, re
import argparse
parser = argparse.ArgumentParser(description='read a vcf file and find the SNPs private to certain populations, only for bianry snps, not for snps with both alleles different from reference')

parser.add_argument("--in","-i", dest="infile", help="vcf file", required=True)
parser.add_argument("--pops","-p", dest="pops", help="comma separated list of samples of a certain population, populations separated by colons (eg: \"1,2,3:4,5,6:7,8,9\") ", required=True)
parser.add_argument("--clades","-c", dest="clades", help="comma separated list of names for clades  (eg: \"I-III,V,VI\") ", default="I-III,V,VI")
parser.add_argument("--gq", dest="gq", help="minimal genotype quality (default: 30) ", default=30)
parser.add_argument("-n", dest="num", help="maximal number of clades a gentoype can be private to (default: 1) ", default=1)


args = parser.parse_args()
infile = vars(args)['infile']
pops= vars(args)['pops'].split(":")
pops = [[int(y) for y in x.split(",")] for x in pops ]
num = int(vars(args)['num'])
clades=vars(args)['clades'].split(",")
gq=int(vars(args)['gq'])
print "#chrom\tBPS\tRef:Alt\t"+"\t".join(clades)
with open(infile,"r") as f:
    for line in f:
        if re.match("\#",line):
            continue
        line=line.rstrip()
        fields=line.split("\t")
        if ( len(fields[3]) != len(fields[4]) ) or len(fields[3]) != 1 or not ( fields[4] in "ACTG") :
            continue
        chrom=fields[0]
        bps=fields[1]
        alleles=fields[3:5]
        pop_all=[]
        for pop in pops:
            pop_all.append(sample_allele(pop,fields[8:],gq))
        if -1 in pop_all:
            continue
        if min(pop_all.count(0),pop_all.count(1)) <= num and ( min(pop_all.count(0),pop_all.count(1)) >= 1 ):
            print "{}\t{}\t{}:{}\t{}".format(chrom,bps,alleles[0],alleles[1],"\t".join([ alleles[x] for x in pop_all]))
