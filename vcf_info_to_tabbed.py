import vcf
import sys, re
import os 
import argparse
import select

parser = argparse.ArgumentParser(description="""Go through a VCF file and write a tabbed version of SNPs consisting of:
CHR BPS ALLELES INFO1 INFO2 INFO3
Attention: only takes first entry of lists, does not print more than two alleles
For piping use STDIN as a vcf filename
""")

parser.add_argument("--in","-i", dest="vcffile", help="vcf-file", required=True)
parser.add_argument("--tags","-t", dest="tags", help="comma separated list of INFO fields to be printed, to get specific fields addcolon sperated entries indexes (zero based, eg. PV4:0:1:2:3, by default only first entry is output) (def: \"AF,BaseQRankSum,DP,Dels,FS,MQ,MQRankSum,QD,ReadPosRankSum\")", default="AF,BaseQRankSum,DP,Dels,FS,MQ,MQRankSum,QD,ReadPosRankSum")

args = parser.parse_args()
vcf_file = vars(args)['vcffile']
tags= vars(args)['tags'].split(",")

# open vcf file
if vcf_file == "STDIN":
    inf = sys.stdin
    vcf_reader = vcf.Reader(inf)
else:
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

# create dict of fields to output
tag_dict={ x.split(":")[0]:[int(z) for z in x.split(":")[1:]]  for x in tags if (len(x.split(":")) > 1)}
# create header string
tags_disp= [ x if len(x.split(":")) == 1 else "\t".join([ x.split(":")[0]+"_"+z for z in x.split(":")[1:]] )  for x in tags]
# clean list of tags to remove :
tags = [ x.split(":")[0] for x in tags ]
print "#CHR\tBPS\tALL\t"+"\t".join(tags_disp)

for entry in vcf_reader:
    #if not entry.is_snp:
    #    continue
    alleles = "".join([ str(x) for x in entry.alleles[:2]])
    info_fields=[]
    for tag in tags:
        try:
            ifield=entry.INFO[tag]
        except:
            if tag in tag_dict.keys():
                ifield="\t".join(["NA" for x in  tag_dict[tag]])
            else:
                ifield="NA"
        if type(ifield)==list:
            if tag in tag_dict.keys():
                ifield="\t".join([ str(ifield[x]) for x in tag_dict[tag]])
            #ifield=",".join([ str(x) for x in ifield])
            else:
                ifield=str(ifield[0])
        info_fields.append(str(ifield))
    print entry.CHROM+"\t"+str(entry.POS)+"\t"+alleles+"\t"+"\t".join(info_fields)

    
    
    
