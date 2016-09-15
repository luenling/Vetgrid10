#import Sync_parser
import sys, os, re
import gzip
import numpy as np
from scipy.stats.stats import nanmean
import collections
from scipy.stats import fisher_exact

import argparse

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
            self.strand_bias=None
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
    
    def get_strand_bias(self, two=True):
        """
        returns an array with three entries: first and array with all individual strandbias fisher pV, second the total  strandbias fisher pV, and last an array with the ratios of the supporting reads on the less common strand to the more common strand for the major and minor alleles. If there is anything going wrong with the first fisher test - maybe not a true snp or something messy, it returns False.
        """
        if not self.strands or not two:
            # no stranded information read
            return False
        if not self.majminor:
            self.get_two_major_alleles()
        whichcols = self.majminor
        num_pop=len(self.seqs)
        fw_rev = np.array([ [ self.fw_strand[indx][whichcols], self.seqs[indx][whichcols] - self.fw_strand[indx][whichcols] ] for indx in range(0,num_pop)],dtype=int)
        fw_rev_tot = fw_rev.sum(axis=0)
        try:
            pop_pv = [ fisher_exact(x)[1] for x in fw_rev ]
        except:
            self.strand_bias=False
            return self.strand_bias
        tot_pv = fisher_exact(fw_rev_tot)[1]
        min_supp_read = np.divide(fw_rev_tot.min(axis=0),fw_rev_tot.max(axis=0),dtype=float)
        self.strand_bias = [ pop_pv, tot_pv, min_supp_read ]
        return self.strand_bias

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



#from scipy import stats
#Author: Lukas Endler
parser = argparse.ArgumentParser(description='read a sync or cmhout file and create a file with the allele frequencies. the major/minor alleles are calculated  over all populations. the resulting file will contain the major allele freq and the coverage for each population.') 
parser.add_argument("-i","--infile", dest="infile", help="sync/cmhout file", required=True)
parser.add_argument("--all", dest="two",  action="store_false", help="output all allelefreqs and not just the major allele freq (calculated from just the major and minor counts) (default: True)", default=True)
parser.add_argument("-c","--cov", dest="cov",  action="store_true", help="output coverages too (default: False)", default=False)
parser.add_argument("-o","--out", dest="outfile", help="output sync file, if \"stdout\" redircet to STDOUT (default: infile.af)", default=None)
parser.add_argument("-p","--pV", dest="pV",  action="store_true", help="output P values (default: False)", default=False)
parser.add_argument("-a","--avg", dest="avg", help="create averages of afs in populations in each files eg. 1-3,4-6 (default: None)", default=None)
parser.add_argument("--pol", dest="pol", help="polarise afs. according rising allele according to population pairs, only effect if not averaging or outputting all allelefrwqs,  eg \"1:4,2:5,3:6\" (default: None)", default=None)
parser.add_argument("--ref", dest="ref", action="store_true", help="also print reference allele (default: False)", default=False)
parser.add_argument("--mc","--min_count", dest="min_count", help="integer indicating the minimal count over all populations above which to consider allels (default: 0)", default=0)
parser.add_argument("--mf", "--min_frac", dest="min_fract", help="float indicating the minimal allele frequency above which an allele is considered in each individual population (default: 0)", default=0)

args = parser.parse_args()
infile = vars(args)['infile']
outfile = vars(args)['outfile']
two = vars(args)['two']
pV = vars(args)['pV']
avg= vars(args)['avg']
ref= vars(args)['ref']
pol= vars(args)['pol']
cov= vars(args)['cov']
min_count= int(vars(args)['min_count'])
min_fract= float(vars(args)['min_fract'])


if pol:
    pol=[ x.split(":") for x in pol.split(",")]
    pol=[ [int(x[0])-1,int(x[1])-1]  for x in pol]

if avg:
    avg=[ x.split("-") for x in avg.split(",")]
    avg=[ [int(x[0])-1,int(x[1])-1]  for x in avg]

if outfile == None:
    outfile = infile + ".af"
if re.search("\.b?gz",infile):
    inf = gzip.open(infile,'rb')
else:
    inf = open(infile,"r")
if outfile == "stdout" or outfile == "STDOUT":
    out=sys.stdout
else:
    out = open(outfile,"w")
nucs = [x for x in "ATCG"]
for line in inf:
    entry = SyncLineParser(line,min_count=min_count,min_fract=min_fract)
    afs = entry.get_pop_allele_freqs(two)
    maj_min = "".join([ nucs[x] for x in  entry.get_two_major_alleles() ])
    all_pos=0 # 0: major, 1: minor
    if cov:
        coverages = entry.get_pop_coverages(two) 
    if pol:
        # see which direction of change is most frequent over replicates in the major allele (1: increasing,0,-1: decreasing)
        polarisation = np.median(np.sign(np.array([ afs[x[0]][0]-afs[x[1]][0]  for x in pol ])))
        #freqs = "\t".join([ ":".join([ str(y) for y in x ]) for x in afs ])
        if polarisation < 0:
            maj_min=maj_min[::-1]
            all_pos=1                
    if avg:
        # get the averages by first creating a matrix, averaging, and then recreating vectors
        if cov:
            coverages2 =  [ np.matrix(coverages[pos[0]:(pos[1]+1)]) for pos in avg]
            coverages =  np.array([ x.mean() for x in coverages2])
        #print afs
        try:
            afs2 = [ np.matrix(afs[pos[0]:(pos[1]+1)]) for pos in avg]
        except:
            print >> sys.stderr, "Problem with generating freq matrix for\n"+line+"\n"+str(afs)+"\nSkipping line"
            continue
        #print afs2
        afs = [ np.array(nanmean(x,0)).flatten() for x in afs2]
        #print afs
        #print avg
    if two:
        if cov: 
            freqs = "\t".join([ str(x[all_pos]) + "\t" + str(y) for x,y in zip(afs,coverages) ])
        else:
            freqs = "\t".join([ str(x[all_pos]) for x in afs ])
    else:
        freqs = "\t".join([ ":".join([ str(y) for y in x ]) for x in afs ])
    if pV and entry.cmhp:
        freqs += "\t"+str(entry.cmhp)
    if ref:
        print >>out, "{0}\t{1}\t{2}\t{3}\t{4}".format(entry.chr,entry.pos,entry.ref,maj_min,freqs)
    else:
        print >>out, "{0}\t{1}\t{2}\t{3}".format(entry.chr,entry.pos,maj_min,freqs)
inf.close()
out.close()

