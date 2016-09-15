import sys, os, re
import pybedtools
from collections import defaultdict 
import argparse
parser = argparse.ArgumentParser(description='reads a gff and a vcf file and creates two sync like files, one with the genes containing the count of SNPs and one with SNPS and the genes they are in.') 

parser.add_argument("--gff","-g", dest="gff", help="gff file with genes/CDS", required=True)
parser.add_argument("--vcf","-v", dest="vcf", help="vcf file with SNPs", required=True)
parser.add_argument("--type","-t", dest="feat_types", help="types to use for overlapping (comma seperated list)", default="gene")
parser.add_argument("--verb", action="store_true", dest="verb", help="set verbouse output", default=False)
args = parser.parse_args()
gff_fn = vars(args)['gff']
vcf_fn=vars(args)['vcf']
feat_types = set(vars(args)['feat_types'].split(","))
gff = pybedtools.BedTool(gff_fn)
vcf = pybedtools.BedTool(vcf_fn)
verb =vars(args)['verb']
if verb:
    print "intersecting "+vcf_fn+" "+gff_fn
gff_in_vcf=gff.intersect(vcf,c=True).filter(lambda x: int(x.fields[-1]) > 0 and x.fields[2] in feat_types ).saveas()

genes_sync=vcf_fn+".genes.sync"
of = open(genes_sync,'w')
if verb:
    print "writing file " + genes_sync
print >>of, "CHR\tSTART\tSTOP\tGENE\tSNPs"
for feature in gff_in_vcf:
    print >>of, "{}\t{}\t{}\t{}\t{}".format(feature.chrom,feature.start+1,feature.end,feature.attrs['ID'],feature[-1])
of.close()

if verb:
    print "intersecting "+vcf_fn+" "+gff_fn
gff_in_vcf=gff.intersect(vcf).filter(lambda x: x.fields[2] in feat_types )
snps_sync=vcf_fn+".genes_per_SNP.sync"
snp_dict=defaultdict(lambda: defaultdict(list))
for feature in gff_in_vcf:
    snp_dict[feature.chrom][feature.end].append(feature.attrs['ID'])
if verb:
    print "writing "+snps_sync
    
of = open(snps_sync,"w")
for chrom in sorted(snp_dict.keys()):
    for pos in sorted(snp_dict[chrom].keys()):
        print >>of, "{}\t{}\t{}".format(chrom,pos,"\t".join(sorted(set(snp_dict[chrom][pos]))))
of.close()
