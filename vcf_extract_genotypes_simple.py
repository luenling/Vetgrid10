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
CHR BPS ALLELES TotalGT HomoRef Hetero
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
    #vcf_reader = vcf.Reader(inf)
else:
    if re.search("\.b?gz",vcf_file):
        inf = gzip.open(vcf_file,'rb')
    else:
        inf = open(vcf_file,"r")
    #vcf_reader = vcf.Reader(open(vcf_file, 'r'))
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NONE    ind_1   ind_10  ind_11  ind_12  ind_2   ind_3   ind_4   ind_5   ind_6   ind_7   ind_8   ind_9
# get set of samples to consider
samples=False
sample_idx=[]
for entry in inf:
    if re.match("\#\#",entry):
        continue
    entry=entry.rstrip()
    if not samples and re.match("\#CHROM",entry):
        samp_dict=dict(zip(entry.split()[9:],range(9,len(entry.split()))))
        for exc in exclude:
            if exc in samp_dict.keys():
                del(samp_dict[exc])
        #samples=sorted(samples.values()) # or sorted([ samples[x] for x in samples.keys()])
        samples=sorted(samp_dict.keys())
        sample_idx=[ samp_dict[x] for x in samples ]
        print "CHR\tBPS\tALLELES\ttot_cal\thomR\thet\thomA\tFi\tfR\tfA\tHomExp\tHetExp\tfRt\tfAt\tFit\t"+"\t".join(samples)
        continue
    if re.match("\#",entry):
        continue
    # now work on the vcf entries
    entries=entry.split()
    if ( len(entries[3]) != 1 ) or ( len(entries[4]) != 1 ):
        # not snp or more than two alleles
        continue
    chrom=entries[0]
    bps=int(entries[1])
    alleles = entries[3] + entries[4]
    if maf:
        # check if SNP in dict and if same major/minor
        if not( chrom in maf_dict.keys() and bps in maf_dict[chrom].keys()):
            continue
        if not ( alleles[0] in maf_dict[chrom][bps][0] and alleles[1] in maf_dict[chrom][bps][0]):
            continue
        # check if major is reference allele:
        if alleles != maf_dict[chrom][bps][0]:
            maf_dict[chrom][bps][1]=1-maf_dict[chrom][bps][1]
    homR=0 # homozygous reference
    homA=0 # homozygous alternative
    het=0 # heterozygous
    nogt=0 # not called
    gts=[]
    for idx in range(0,len(sample_idx)):
        gt=entries[sample_idx[idx]].split(":")[0]
        if gt == "1/1":
            homA += 1
            gts.append("1.0")
        elif gt == "0/1" or gt == "1/0":
            gts.append("0.5")
            het +=1
        elif gt == "0/0":
            gts.append("0.0")
            homR +=1
        else:
            gts.append("NA")
            nogt+=1
    tot_cal=len(samples)-nogt # total number of called samples
    if tot_cal == 0:
        continue
    assert het==tot_cal-(homR+homA), "some thing went wrong counting the genotypes at"+ chrom+"\t"+str(bps)
    fR=(homR+het/2.0)/tot_cal
    fA=(homA+het/2.0)/tot_cal
    assert abs(1 -(fR + fA)) < 1e-3, "fR and fA do not sum to 1 at"+ chrom+"\t"+str(bps)
    # reduction of heterozygosity from hardy weinberg, F coefficient of inbreeding
    if het == 0:
        Fi = 1.0
    else:
        Fi=1.0-(float(het)/tot_cal)/(2.0*fR*fA)
    tot_pop=""
    if maf:
        # calculate population total
        fRt=maf_dict[chrom][bps][1]
        fAt=1-fRt
        HomExp=tot_cal*(fRt**2+fAt**2)
        HetExp=tot_cal*(2*fRt*fAt)
        try:
            Fit=1.0-(float(het)/tot_cal)/(2.0*fRt*fAt)
        except:
            Fit=1.0
        tot_pop="\t{:.2}\t{:.2}\t{:.3}\t{:.3}\t{:.3}".format(HomExp,HetExp,fRt,fAt,Fit)
    genotypes="\t".join(gts)
    print chrom+"\t"+str(bps)+"\t"+alleles+"\t"+"{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}".format(tot_cal,homR,het,homA,Fi,fR,fA)+tot_pop+"\t"+genotypes

# for problem of forgotten tab:
#  gsed 's/\([\.e\-]\+[0-9]\+\|nan\)\([01]\.[05]\|NA\)/\1\t\2/' BGI_106_snps_unified_gt_recal99.het_b.sync > BGI_106_snps_unified_gt_recal99.het_b_tab.sync 

    

