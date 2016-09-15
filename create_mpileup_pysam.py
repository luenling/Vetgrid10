import pysam
import numpy as np
from scipy import stats
import re,sys
#from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description="""
gets a list of bam files and a sync file with snps (CHR BPS AB) and prints a file with the following fields for each SNP:
CHR BPS AB F_A:R_A:F_B:R_B:mPos_A:mPos_B:mIns_A:mIns_B:SC_A:SC_B
""") 
parser.add_argument("-b", dest="bam_files", help="bam files, comma separated list", required=True)
parser.add_argument("-s", dest="snp_file", help="snp file", required=True)
parser.add_argument("-v","--verb", dest="verb", action="store_true", help="print verbouse state on stderr", default=False)
parser.add_argument("-q","--bq", dest="bqual", help="minimal base quality", default=20)
parser.add_argument("-m","--mq", dest="mqual", help="minimal mapping quality", default=20)
#parser.add_argument("--bcf", dest="bcf", action="store_true", help="create bcf file", default=False)

args = parser.parse_args()
bam_files = vars(args)['bam_files'].split(",")
snp_file= vars(args)['snp_file']
verb = vars(args)['verb']
bqual = int(vars(args)['bqual'])
mqual = int(vars(args)['mqual'])
#bcf = vars(args)['bcf']

#bam_files=" ".join(bam_files)
snp_fh = open(snp_file,"r")
# if bcf:
#     first=True

for line in snp_fh:
    if re.match("\#",line):
        continue
    line=line.rstrip()
    fields=line.split()
    chrom=fields[0]
    bps=int(fields[1])
    maj=fields[2][0].capitalize()
    min=fields[2][1].capitalize()
    snp_region=chrom+":"+fields[1]+"-"+fields[1]
    # if bcf:
    #     if first:
    #         for i in pysam.mpileup('-u','-s','-q',str(mqual),'-Q',str(bqual),'-r',snp_region,*bam_files):
    #             sys.stdout.write(i)
    #         first=False
    #     else:
    #         sys.stdout.write("\n"+pysam.mpileup('-u','-s','-q',str(mqual),'-Q',str(bqual),'-r',snp_region,*bam_files)[-1])
    mpileup_line=pysam.mpileup('-q',str(mqual),'-Q',str(bqual),'-r',snp_region,*bam_files)[0].split("\t")
    seqs=[ mpileup_line[4+3*i] for i in range(0,len(bam_files))]
    fw_rev = np.array([ [ [x.count(maj),  x.count(min) ] , [ x.count(maj.lower()),  x.count(min.lower()) ] ] for x in seqs],dtype=int)
    fw_rev_tot = fw_rev.sum(axis=0)
    pop_pv = "\t".join([ "{:.3g}".format(stats.fisher_exact(x)[1]) for x in fw_rev ])
    tot_pv = "{:.3g}".format(stats.fisher_exact(fw_rev_tot)[1])
    min_supp_read = str(fw_rev_tot.min())
    print chrom+"\t"+fields[1]+"\t"+fields[2]+"\t"+pop_pv+"\t"+tot_pv+"\t"+min_supp_read
#python /Volumes/Temp/Lukas/Tools/Scripts/create_mpileup_pysam.py -s females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_mct8_mcv25_fdr.sync.strandbias.af  -b "/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/9167.picnodupl.filtered.mq20_chrfilt_RG_sanger_real.bam,/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/9168.picnodupl.filtered.mq20_chrfilt_RG_sanger_real.bam,/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/9390.picnodupl.filtered.mq20_chrfilt_RG_sanger_real.bam,/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/pop1_earlylate_merged_RG_sanger_real.bam,/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/pop2_earlylate_merged_RG_sanger_real.bam,/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/pop3_earlylate_merged_RG_sanger_real.bam" 
