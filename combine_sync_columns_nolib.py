#need sync parser
class SyncLineParser(object):
    def __init__(self,syncline, replicates=None, includedel=False, strands=False, min_count = 0, min_fract = 0):
        '''
        SyncLineParser
        pass it a single line from a sync file
        optional options:
            which (default=None) will indicate which columns have the phenotype of interest
            includedel  (default=False) set to True indicates that the last two columns of sync should be included
        other attributes set on init:
                self.whichcols #
                self.min_count : all alleles having a count < mincount over all populations are set to zero (default:0) 
                self.min_fract : all alleles having a frequency below min_fract in a population are set to zero in that population (default: 0)
                self.chr 
                self.pos 
                self.ref 
                self.seqs #sync columns with bp information
                self.cmhp
        functions:
            self.get_two_major_alleles()
                sets self.majminor, gets the two major alleles, ignoring anything higher than a di-allellic site
            self.get_pop_allele_freqs(two=True)
                gets overall allele freqs from seq information across all replicates
                if two is set to True, major and minor alleles only
            self.get_overall_allele_freqs(two=True)
                gets overall allele freqs from seq information across all replicates
                if two is set to True, major and minor alleles only
            self.get_reduced_allele_freq_dict(two=True)
                gets overall allele freqs from seq information across all replicates
                if two is set to True, major and minor alleles only
            
        '''
        self.includedel = includedel
        self.strands = strands
        if self.includedel:
            self.ncol = 6
        else:
            self.ncol = 4
        self.whichcols=replicates #labels for replicates
        #include last two columns if True
        #parse syncline
        sline = syncline.split()
        self.min_count = min_count
        self.min_fract = min_fract
        self.chr =sline.pop(0)
        self.pos =sline.pop(0)
        self.ref =sline.pop(0)
        self.seqs = [ np.array([int (x) for x in y.split(':')[:self.ncol] ]) for y in sline if ':' in y]
        if self.strands:
            self.fw_strand = [ np.array([int (x) for x in y.split(':')[6:6+4] ]) for y in sline if ':' in y]
        # get rid of alleles with count less than min_count
        if min_count:
            self.set_seqs_mc()
        if min_fract:
            self.set_seqs_mf()
        if ':' not in sline[-1]: #CMH output
            self.cmhp =sline.pop()
        else:
            self.cmhp =None
        #make dictionary with information for phenotype or replicate
        self.majminor = None
        self.coverages = None
        self.pop_allele_freqs =None
        self.overall_allele_freqs=None
        self.seq_dict =None
        self.reduced_seq = None
        self.reduced_seq_dict = None
        self.reduced_af_dict = None
        

    def set_seqs_mc(self):
        seq_counts = [ np.array([ x[y] for x in self.seqs ]).sum() for y in range(0, self.ncol) ]
        for seq in self.seqs:
            for i,count in enumerate(seq_counts):
                if count > 0 and count < self.min_count:
                    seq[i] = 0
        return self.seqs

    def set_seqs_mf(self):
        fracts = [ x[:self.ncol]/float(x[:self.ncol].sum()) for x in self.seqs ]
        for i,seq in enumerate(self.seqs):
            for j in range(0, self.ncol) :
                if seq[j] > 0 and ( fracts[i][j] < self.min_fract):
                    seq[j] = 0
        return self.seqs

        
    def get_seq_dict(self):
        if self.seq_dict == None:
            self.seq_dict = collections.defaultdict(list)
            for r in range(0, len(self.seqs)):
                self.seq_dict[self.whichcols[r]].append(self.seqs[r])
        return self.seq_dict

    def get_two_major_alleles(self):
        if not self.majminor:
            whichcols = range(0, self.ncol)
            allele_totals = np.array([sum([y[x] for y in self.seqs]) for x in whichcols])
            # get the highest ranked columns in reverse order
            if sorted(allele_totals)[-2] > 0:
                self.majminor =list(allele_totals.argsort()[-1:-3:-1])
            else:
                self.majminor = [ allele_totals.argsort()[-1] ]
        return self.majminor
    def get_pop_coverages(self, two = True):
        if not self.coverages:
            if two and not self.majminor:
                self.get_two_major_alleles()
            if two:
                whichcols = self.majminor
            else:
                whichcols = range(0, self.ncol)
            reduced_seq = np.array([[float(x[y]) for y in whichcols] for x in self.seqs]) #reduce
            self.coverages =  [ x.sum() for x in reduced_seq]
        return self.coverages

    def get_pop_allele_freqs(self, two = True):            
        if not self.pop_allele_freqs:
            if two and not self.majminor:
                self.get_two_major_alleles()
            if two:
                whichcols = self.majminor
            else:
                whichcols = range(0, self.ncol)
                #print whichcols
            reduced_seq = np.array([[float(x[y]) for y in whichcols] for x in self.seqs]) #reduce
            #pop_totals =  [ x.sum() for x in reduced_seq]
            pop_totals = self.get_pop_coverages(two)
            self.pop_allele_freqs =[]
            for i in range(len(self.seqs)):
                if (pop_totals[i] > 0):
                    freq=reduced_seq[i]/pop_totals[i]
                else:
                    freq=np.zeros(len(reduced_seq[i]))+np.nan
                self.pop_allele_freqs.append(freq)
        return self.pop_allele_freqs
    
    def get_overall_allele_freqs(self, two=True):   
        if not self.overall_allele_freqs:
            self.overall_allele_freqs=[]
            if not self.pop_allele_freqs:
                self.get_pop_allele_freqs(two)
            num_pop=len(self.pop_allele_freqs)
            self.overall_allele_freqs = [sum([y[x] for y in self.pop_allele_freqs])/num_pop for x in range(0, len(self.pop_allele_freqs[0]))]
        return self.overall_allele_freqs
    
    def get_reduced_seq (self, two = True):
        if self.reduced_seq == None:
            if not self.majminor:
                self.get_two_major_alleles()
            self.reduced_seq = np.array([x[self.majminor] for x in self.seqs]) #reduce 
        return self.reduced_seq

    def get_reduced_seq_dict(self):
        if self.reduced_seq_dict == None:
            if self.reduced_seq == None:
                self.get_reduced_seq()
            self.reduced_seq_dict = {}
            for r,seq in enumerate(self.reduced_seq):
                rep = self.whichcols[r]
                if rep in self.reduced_seq_dict:
                    self.reduced_seq_dict[rep]=np.vstack([self.reduced_seq_dict[rep],seq])
                else:
                    self.reduced_seq_dict[rep]=seq
        return self.reduced_seq_dict
            
    def get_reduced_allele_freq_dict(self,two=True):
        if self.reduced_af_dict == None:
            if not self.pop_allele_freqs:
                self.get_pop_allele_freqs(two)
            self.reduced_af_dict = collections.defaultdict(list)
            for r,af in enumerate(self.pop_allele_freqs):
                rep = self.whichcols[r]
                self.reduced_af_dict[rep].append(af)
        return self.reduced_af_dict
    

        # whichcols = syncline.majminor
        # num_pop=len(syncline.seqs)
        # cont_table=np.zeros((2,2),dtype=int)
        # for indx in range(0,num_pop):
        #     cont_table += np.array([ syncline.fw_strand[indx][whichcols], syncline.seqs[indx][whichcols] - syncline.fw_strand[indx][whichcols] ])      
# for line in inf:
#     syncline = Sync_parser.SyncLineParser(line,strands=True)
#     syncline.get_strand_bias()
#     print str(syncline.strand_bias)
#     print str(syncline.seqs)
#     print str(syncline.fw_strand)

def combine_pops(pops,cov_p_ind,seqs):
    cov_p_ind=np.array([ cov_p_ind[x] for x in pops])
    seqs=[ seqs[x] for x in pops]
    ratio_cpi=cov_p_ind.min()/cov_p_ind
    scaled_counts=[ seqs[x]*ratio_cpi[x] for x in range(0,len(ratio_cpi)) ]
    return(np.array(sum(scaled_counts),dtype=int))

import sys, os, re
import numpy as np
#import Sync_parser
from numpy import array, zeros
import math
import gzip
import argparse
#import rpy2.robjects as robjects
#from rpy2.robjects.packages import importr

parser = argparse.ArgumentParser(description="""
read sync file and combine two columns each, depending on the coverage per individual.
For each population combination it calculates the coverage per individual and the ratio of that to the minimal coverage per individual and then adds the nucleotide counts times that scaling factor.
This should give you the individual weighted mean of the allele frequencies at 2 times the coverage per indiviudal of the population with the lower coverage per individual.  
It prints a sync file with combined populations.
""")

parser.add_argument("--pops","-p", dest="pops", help="the population pairs to combine \"|\" (eg. \"1+2,2+3,4+5\")",default=False)
parser.add_argument("--ind","-n", dest="inds", help="individuals per population (eg. \"65,75,100,140,110,202\")",default=False)

parser.add_argument("--in","-i", dest="infile", help="sync file to be read", required=True)
args = parser.parse_args()
infile = vars(args)['infile']
inds= np.array([ int(x) for x in vars(args)['inds'].split(",") ])
pops = [ [ int(y) -1 for y in x.split("+") ]  for x in vars(args)['pops'].split(",") ]

if re.search("\.b?gz",infile):
    inf = gzip.open(infile,'rb')
else:
    inf = open(infile,"r")

for line in inf:
    if re.match("\s*\#",line):
        # comment lines
        continue
    line=line.rstrip()
    syncline = SyncLineParser(line)
    cov_p_ind=syncline.get_pop_coverages()/inds
    pop_counts = [ combine_pops(x,cov_p_ind,syncline.seqs)  for x in pops ]
    pop_countstr = "\t".join([  ":".join([ str(y) for y in x]) for x in pop_counts ])
    print syncline.chr+"\t"+str(syncline.pos)+"\t"+syncline.ref+"\t"+pop_countstr
inf.close()
        
                                                
