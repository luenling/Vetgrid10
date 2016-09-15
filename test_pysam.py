import pysam
import numpy as np
from scipy import stats
import re
from collections import defaultdict
import argparse

def get_pos_allele_stats(read, pos, q_enc=33):
    cigars=read.cigar;
    seq_pos=0
    soft_clip=0 
    pos_diff = pos - (read.pos + 1) 
    strand=1 - 2*int(read.is_reverse) # 1 for forward, -1 for reverse
    insert=abs(read.tlen) 
    #print str(pos_diff)
    #cig_types=np.zeros() # [S,H,]
    for (cig_type,cig_len) in cigars:
        if pos_diff > 0:
            if (cig_type in (0,4,5,7,8) ): # M, S, H, =, X 
                if pos_diff > cig_len:
                    pos_diff -= cig_len
                    seq_pos += cig_len
                else: 
                    seq_pos += pos_diff
                    pos_diff = 0
            elif (cig_type == 1 ): #insertions
                    seq_pos += cig_len
                    #rseq=rseq[0:seq_pos]+rseq[seq_pos+cig_len:]
            elif (cig_type in (2,3)): #deletion, skip (D,N)
                pos_diff -= cig_len
                if pos_diff == 0:
                    pos_diff = -1
                #rseq=rseq[0:seq_pos]+cig_len*"-"+rseq[seq_pos:]
        if cig_type == 4:
            soft_clip = 1
        #print str(pos_diff)
    read_frac=0.0
    read_pos=0
    if pos_diff == 0:
        allele=read.seq[seq_pos]
        qual=ord(read.qual[seq_pos])
        read_pos=seq_pos
        read_frac=(seq_pos+1.0)/read.rlen
        if (strand == -1):
            read_pos = read.rlen - read_pos
            read_frac = 1 - read_frac 
    elif pos_diff < 0: # position at deletion
        allele="-"
        qual=0
    else: # pos_diff > 0 - not reached with read
        allele="n"
        qual=0
    return (allele,qual-q_enc,strand,soft_clip,read_frac,abs(read.tlen),read_pos)


def get_snp_stats(read_stats,ref,alt):
    snp_vals={}
    snp_stats={}
    snp_vals[ref]=[]
    snp_vals[alt]=[]
    for i in read_stats:
        if i[0] in (ref,alt):
            snp_vals[i[0]].append(i[1:])
    for x in snp_vals.keys():
        if snp_vals[x]:
            snp_vals[x]=np.array(snp_vals[x])
            # fw, rev, ( mean qual, sd) , num SC, (mean read frac, sd), (mean insert size, sd),(mean read pos, sd) 
            snp_stats[x]=[ (snp_vals[x].shape[0] + snp_vals[x][:,1].sum())/2,  (snp_vals[x].shape[0] - snp_vals[x][:,1].sum())/2 , ( snp_vals[x][:,0].mean() , snp_vals[x][:,0].std() ),  snp_vals[x][:,2].sum(), ( snp_vals[x][:,3].mean() , snp_vals[x][:,3].std() ) ,  ( snp_vals[x][:,4].mean(),  snp_vals[x][:,4].std() ) , ( snp_vals[x][:,5].mean() , snp_vals[x][:,5].std() )]
        else:
           snp_stats[x] = [] 
    return snp_stats


def strand_bias(maj_fw,maj_rev,min_fw,min_rev):
    try:
        sb=abs(float(min_fw)/(maj_fw+min_fw) - float(min_rev)/(maj_rev+min_rev) )/float(min_fw+min_rev)*(maj_fw+maj_rev+min_fw+min_rev)
    except ZeroDivisionError :
        sb="Inf" 
    return sb

def strand_bias_gatk(maj_fw,maj_rev,min_fw,min_rev):   
    try:
        sb= max( float(min_fw) /(maj_fw+min_fw)*float(maj_rev)/(maj_rev+min_rev), float(min_rev)/(maj_rev+min_rev)*float(maj_fw)/(maj_fw+min_fw) ) / (maj_fw+maj_rev) * (maj_fw+min_fw+maj_rev+min_rev)
    except ZeroDivisionError :
        sb="Inf" 
    return sb

def strand_bias_fisher(maj_fw,maj_rev,min_fw,min_rev):
    return stats.fisher_exact([[maj_fw,min_fw],[maj_rev,min_rev]])[1]
    

parser = argparse.ArgumentParser(description="""
gets a list of bam files and a sync file with snps (CHR BPS AB) and prints a file with the following fields for each SNP:
CHR BPS AB F_A:R_A:F_B:R_B:mPos_A:mPos_B:mIns_A:mIns_B:SC_A:SC_B
""") 
parser.add_argument("-b", dest="bam_files", help="bam files, comma separated list", required=True)
parser.add_argument("-s", dest="snp_file", help="snp file", required=True)
parser.add_argument("-v","--verb", dest="verb", action="store_true", help="print verbouse state on stderr", default=False)

args = parser.parse_args()
bam_files = vars(args)['bam_files'].split(",")
snp_file= vars(args)['snp_file']
verb = vars(args)['verb']

bam_fh = [ pysam.Samfile(x,"rb") for x in bam_files ]

snp_fh = open(snp_file,"r")

for line in snp_fh:
    if re.match("\#",line):
        continue
    line=line.rstrip()
    fields=line.split()
    chrom=fields[0]
    bps=int(fields[1])
    maj=fields[2][0]
    min=fields[2][1]
    outline=""
    snp_stat_tot={}
    read_stats_tot=[]
    #print chrom+"\t"+str(bps)+"\t"+fields[2]
    for pop in bam_fh:
        snp_reads=[ x for x in pop.fetch(reference=chrom, start=bps-1, end=bps)]
        readstats=[ get_pos_allele_stats(x,bps) for x in snp_reads]
        #print str(len(readstats))
        #snp_stats=get_snp_stats(readstats,maj,min)
        read_stats_tot.extend(readstats)
    snp_stats_tot=get_snp_stats(read_stats_tot,maj,min)
    if not ( snp_stats_tot[maj] and  snp_stats_tot[min] ) :
        continue
    ( maj_fw, maj_rev)=snp_stats_tot[maj][0:2]
    ( min_fw, min_rev)=snp_stats_tot[min][0:2]
    sb1=strand_bias(maj_fw,maj_rev,min_fw,min_rev)
    sb2=strand_bias_gatk(maj_fw,maj_rev,min_fw,min_rev)
    sbf=strand_bias_fisher(maj_fw,maj_rev,min_fw,min_rev)
    (pos_maj,pos_maj_std)=snp_stats_tot[maj][4]
    (pos_min,pos_min_std)=snp_stats_tot[min][4]
    sc=snp_stats_tot[maj][3]+snp_stats_tot[maj][3]
    maj_mean_insert=snp_stats_tot[maj][5][0]
    min_mean_insert=snp_stats_tot[min][5][0]
    print "{}\t{}\t{}{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrom,bps,maj,min,sb1,sb2,sbf,pos_maj,pos_maj_std,pos_min,pos_min_std,sc,maj_mean_insert,min_mean_insert)
    
        # print the different fields into a string
    #
        
    

# bam_file="/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/9167.picnodupl.filtered.mq20_chrfilt_RG_sanger_real.bam"
# samFile = pysam.Samfile(bam_file,"rb")
# ref_genome="/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa"
# snp_region="2L:414041-414041"
# snp_chr="2L"
# snp_pos=414041
# # get list of reads for each snp
# snp_reads=[ x for x in samFile.fetch(reference=snp_chr, start=snp_pos-1, end=snp_pos)]
# readstats=[ get_pos_allele_stats(x,snp_pos) for x in snp_reads]
# maj="C"
# min="G"
# get_snp_stats(readstats,maj,min)


# for rd in snp_reads:
#     if rd.alen != rd.rlen:
#         rd_lendiff.append(rd)
#         print rd.cigarstring

