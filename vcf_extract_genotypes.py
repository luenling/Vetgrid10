def read_maf_file(infile,pops,maxfreq):
    """
    reads a maf sync like file and creates a dict with chrom->bps->[alleles,avgmaf]
    """
    snp_dict=defaultdict(defaultdict)  #dictionary of chroms with positions and pVals of snps
    if re.search("\.b?gz",infile):
        inf = gzip.open(infile,'rb')
    else:
        inf = open(infile,"r")
    #fill dictionary
    for line in inf:
        if re.match("\#",line):
            continue
        line=line.rstrip()
        fields=line.split()
        maf=sum([ float(fields[x]) for x in pops ])/len(pops)
        if maf > maxfreq:
            continue
        snp_dict[fields[0]][int(fields[1])]=[fields[2][:2],maf]
    inf.close()
    return snp_dict

import vcf
import sys, re
import os 
import argparse
import select
import gzip
from collections import defaultdict


parser = argparse.ArgumentParser(description="""Go through a VCF file and write out the genotype counts for all biallelic SNPs:
CHR BPS ALLELES TotalGT HomoRef Hetero HomoAlt Fis freqRef freqAlt
if a file with population allelefreqs is given the Fit is calculated using those frequencies additionally: HomExpected HetExpected freqRefPop freqAltPop Fit
For piping use STDIN as a vcf filename
""")

parser.add_argument("--in","-i", dest="vcffile", help="vcf-file", required=True)
parser.add_argument("--exclude","-x", dest="exclude", help="Samples to exclude, comma seperated list def.: \"NONE\"", default="NONE")
parser.add_argument("--allelefreq","-f", dest="maf", help="Sync file with allele frequencies default=False", default=False)
parser.add_argument("--pop","-p", dest="pops", help="population in sync file to use, comma separated list for averaging multiple entries (default=1)", default="1")
parser.add_argument("--minfreq","-m", dest="minfreq",type=float, help="minimal minor allele frequency to consider for Fit calculation (default=0.15)", default=0.15)
args = parser.parse_args()
vcf_file = vars(args)['vcffile']
exclude = vars(args)['exclude'].split(",")
pops=[ int(x)+2 for x in vars(args)['pops'].split(",")]
maf=vars(args)['maf']
minfreq=vars(args)['minfreq']
maxfreq=1-minfreq
if maf:
    maf_dict=read_maf_file(maf,pops,maxfreq)   

# open vcf file
if vcf_file == "STDIN":
    inf = sys.stdin
    vcf_reader = vcf.Reader(inf)
else:
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
# get set of samples to consider
samples=set(vcf_reader.samples)
for exc in exclude:
    if exc in samples:
        samples.remove(exc)
for entry in vcf_reader:
    if not entry.is_snp or len(entry.alleles) > 2:
        # not snp or more than two alleles
        continue
    #   
    alleles = "".join([ str(x) for x in entry.alleles])
    if maf:
        # check if SNP in dict and if same major/minor
        if not( entry.CHROM in maf_dict.keys() and entry.POS in maf_dict[entry.CHROM].keys()):
            continue
        if not ( alleles[0] in maf_dict[entry.CHROM][entry.POS][0] and alleles[1] in maf_dict[entry.CHROM][entry.POS][0]):
            continue
        # check if major is reference allele:
        if alleles != maf_dict[entry.CHROM][entry.POS][0]:
            maf_dict[entry.CHROM][entry.POS][1]=1-maf_dict[entry.CHROM][entry.POS][1]
    homR=0 # homozygous reference
    homA=0 # homozygous alternative
    nogt=0 # not called
    for sam in samples:
        gt=entry.genotype(sam).data.GT
        if not gt:
            nogt+=1
        elif gt == '0/0':
            homR+=1
        elif gt == '1/1':
            homA+=1
    tot_cal=len(samples)-nogt # total number of called samples
    if tot_cal == 0:
        continue
    het=tot_cal-(homR+homA)
    fR=(homR+het/2.0)/tot_cal
    fA=(homA+het/2.0)/tot_cal
    assert abs(1 -(fR + fA)) < 1e-3, "fR and fA do not sum to 1 at"+ entry.CHROM+"\t"+str(entry.POS)
    # reduction of heterozygosity from hardy weinberg, F coefficient of inbreeding
    if het == 0:
        Fi = 1.0
    else:
        Fi=1.0-(float(het)/tot_cal)/(2.0*fR*fA)
    tot_pop=""
    if maf:
        # calculate population total
        fRt=maf_dict[entry.CHROM][entry.POS][1]
        fAt=1-fRt
        HomExp=tot_cal*(fRt**2+fAt**2)
        HetExp=tot_cal*(2*fRt*fAt)
        try:
            Fit=1.0-(float(het)/tot_cal)/(2.0*fRt*fAt)
        except:
            Fit=1.0
        tot_pop="\t{:.2}\t{:.2}\t{:.3}\t{:.3}\t{:.3}".format(HomExp,HetExp,fRt,fAt,Fit)
    print entry.CHROM+"\t"+str(entry.POS)+"\t"+alleles+"\t"+"{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}".format(tot_cal,homR,het,homA,Fi,fR,fA)+tot_pop




    

