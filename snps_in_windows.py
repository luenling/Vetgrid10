def get_B_in_win(winpos,A,B,winsize):
    winhalf=float(winsize)/2.0
    chrom=""
    bp=0
    for x in sorted(A.keys()):
        if winpos < A[x]['upper']:
            chrom=x
            bp=A[x]['bps'][ winpos - A[x]['lower'] ]
            break
    winstop = bp + winhalf
    winstart = bp - winhalf
    num_B = len(B[chrom]['sig'][ (B[chrom]['sig'] >= winstart )  & (B[chrom]['sig'] <= winstop)])
    return num_B

def get_significant(snp_dict,thr):
    """
    adds an array of significant SNPs (p Value < thr) to the chromosome dict. returns number of all significants
    """
    sig_num = 0
    for x in snp_dict.keys():
        signifs = snp_dict[x]['bps'][ snp_dict[x]['pVals'] <= thr ]
        sig_num += len(signifs)
        snp_dict[x]['sig'] = signifs
    return sig_num

def read_cmh_file(infile):
    """
    reads a cmh file and returns a dictionary of chromosomes with a key pointing at a numpy array of bp positions and one at the pValues for each chromsome key:
    eg cmhdict['2L']['bps'] = array[BPs],  cmhdict['2L']['pVals'] = array[pV]
    """
    chroms=[]
    # intermediate lists of values
    bps=[]
    pVals=[] 
    cmh_dict={}  #dictionary of chroms with positions and pVals of snps
    inf = open(infile,"r")
    #load dictionary
    for line in inf:
        line.rstrip()
        fields=line.split()
        if not (fields[0] in chroms):
            chroms.append(fields[0])
            if (len(chroms) > 1):
                cmh_dict[chroms[-2]]={}
                cmh_dict[chroms[-2]]['bps']=np.array(bps)
                cmh_dict[chroms[-2]]['pVals']=np.array(pVals)
                cmh_dict[chroms[-2]]['len']=len(cmh_dict[chroms[-2]]['pVals'])
                bps=[]
                pVals=[]
        pVals.append(float(fields[-1]))
        bps.append(int(fields[1]))
    # set last chromosome
    cmh_dict[chroms[-1]]={}
    cmh_dict[chroms[-1]]['bps']=np.array(bps)
    cmh_dict[chroms[-1]]['pVals']=np.array(pVals)
    cmh_dict[chroms[-1]]['len']=len(cmh_dict[chroms[-1]]['pVals'])
    inf.close()
    return cmh_dict


import sys, os, re
import numpy as np
from scipy import stats
from collections import defaultdict
#np.seterr(all=None, divide='ignore', over=None, under=None, invalid='ignore')
import argparse
parser = argparse.ArgumentParser(description='reads files with SNPs from two tests/populations A and B and looks for enrichment of  significant SNPs of B in windows around significant SNPs of A. For the background propability it uses windows around random SNPs of A') 

parser.add_argument("--inA", dest="inA", help="cmh file of pop A", required=True)
parser.add_argument("--inB", dest="inB", help="cmh file of pop B (assumed to be only significant ones)", required=True)
parser.add_argument("--tA", dest="thrA",  help="float of significance level pValue threshold for A, if 0 uses Bonferroni (default=0)",default="0")
parser.add_argument("--tB", dest="thrB",  help="float of significance level pValue threshold for B, if 0 uses Bonferroni (default=0)",default="0")
parser.add_argument("--win", dest="windowsize",  help="window size (default 1000)",default="1000")
parser.add_argument("--rnd_wins", dest="rnd_wins",  help="number of times to do random windows for getting background frequencies (default 1000)", default="1000")

args = parser.parse_args()
inA = vars(args)['inA']
inB = vars(args)['inB']
thrA = float(vars(args)['thrA'])
thrB = float(vars(args)['thrB'])
winsize = int(vars(args)['windowsize'])
win_half = float(winsize)/2.0
rnd_wins = int(vars(args)['rnd_wins'])
all_snps_A = read_cmh_file(inA)
all_snps_B = read_cmh_file(inB)
numSNPs_A= 0
numSNPs_B= 0
for x in all_snps_A.keys():
    numSNPs_A += len(all_snps_A[x]['bps'])
for x in all_snps_B.keys():
    numSNPs_B += len(all_snps_B[x]['bps'])
if (thrA == 0):
    thrA = 0.05/numSNPs_A
if (thrB == 0):
    thrB = 0.05/numSNPs_B
sig_A = get_significant(all_snps_A,thrA)
sig_B = get_significant(all_snps_B,thrB)
#if (rnd_wins == 0):
#    rnd_wins = 1000 * sig_A
# do window for each significant in A and search for occurences of B
sigB_wins = []
for chrom in all_snps_A.keys():
    for bpos in all_snps_A[chrom]['sig']:
        winstart=bpos-win_half
        winstop=bpos+win_half
        sigB_win = len(all_snps_B[chrom]['sig'][ (all_snps_B[chrom]['sig'] >= winstart )  & (all_snps_B[chrom]['sig'] <= winstop)])
        # set an entry in the array of 
        sigB_wins.append(sigB_win)
## add up the numbers of SNPs found for averaging
#num_wins = 0
#num_B_in_wins = 0
#for i in sorted(sigB_wins.keys(),reverse=True):
#    num_wins += sigB_wins[i]
#    num_B_in_wins += sigB_wins[i]*i
sigB_wins = np.array(sigB_wins)
avg_B_per_win = np.mean(sigB_wins)

# do random windows to get distribution:

# get function to map random number to SNPs in popA
# numbers of SNPs per chromosome
tot_SNPs = 0
for chrom in sorted(all_snps_A.keys()):
    all_snps_A[chrom]['lower']=tot_SNPs
    tot_SNPs += all_snps_A[chrom]['len']
    all_snps_A[chrom]['upper']=tot_SNPs
tot_SNPs_B = 0
for chrom in sorted(all_snps_B.keys()):
    tot_SNPs_B += all_snps_B[chrom]['len']
aver_B_rand = np.zeros(rnd_wins)
for i in range(0,len(aver_B_rand)):
    # create array of random numbers from 0..tot_SNPs length sig_A
    randsnps=np.random.random_integers(0,tot_SNPs-1,sig_A)
    for j in randsnps:
        aver_B_rand[i] += get_B_in_win(j,all_snps_A,all_snps_B,winsize)
aver_B_rand /= sig_A
p_B_in_win = float(len(aver_B_rand[aver_B_rand >= avg_B_per_win]))/len(aver_B_rand)
print "InfileA:{}\tSNPsA:{}\tsign. SNPSsA:{}\tsign. thresholdA:{}".format(inA,tot_SNPs,sig_A,thrA)
print "Infile2:{}\tSNPsB:{}\tsign. SNPSsB:{}\tsign. thresholdB:{}".format(inB,tot_SNPs_B,sig_B,thrB)
print "Windowsize:{}\trandom loops:{}".format(winsize,rnd_wins)
print "average number of significant SNPs from B around significant A:{}\tP:{f}\taverage sig around random A:{}".format(avg_B_per_win,p_B_in_win,np.mean(aver_B_rand))