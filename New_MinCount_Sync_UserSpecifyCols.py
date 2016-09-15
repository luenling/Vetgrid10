import sys
import collections
import operator
import numpy as np
import math
import itertools
#from optparse import OptionParser, OptionGroup
import argparse

#Author: Ray Tobler

##python /Volumes/Temp/Ray/Scripts_Ray/Scripts_other/New_MinCount_Sync_UserSpecifyCols.py --sync /Volumes/Temp/Ray/Seasonal_data/seasonal_filtered.sync --header "N" --cutoff 10 --popcols 1,2,3,4,5,6 --output /Volumes/Temp/Ray/Seasonal_data/seasonal_filtered_mincount10.sync


#########################################################   HELP   #########################################################################
#print
parser = argparse.ArgumentParser(description=
"""
    			D E S C R I P T I O N
Purges all SNPs with minimum count below user specified value. Min count is taken across all populations.
""") 

##########################################################  FUNCTIONS ######################################################################
def divide_allelecounts(allelestring,cov): # function which splits strings by : symbol and returns dictionary with counts associated with a key
	values = map(int,allelestring.split(":")[:4]) # apply int to member after allelestring.split (map is python equivalent of apply in R)
	values=np.array(values)
	min_cov=True
	if (np.sum(values) < cov):
		min_cov = False
	#return dict(zip(["A", "T", "C", "G"], values)) # the output of the function is a dictionary
	return ( values , min_cov ) # return values as array and true or false if < cov


def get_major_and_minor_alleles(group): # get the two most frequent alleles from all pops in a group and returns a tuple of allele and count
	totcov = collections.defaultdict(lambda:0) # we use lambda function in place of actual keys - these will be provided within the loop
	for pop in group:
		for base,values in pop.items():
			if values != "0":
				totcov[base] += pop[base]
	minor=sorted(totcov.iteritems(),key=operator.itemgetter(1))[2]
	major=sorted(totcov.iteritems(),key=operator.itemgetter(1))[3]		
	return [minor,major]

#########################################################   CODE   #########################################################################

parser.add_argument('--sync','-s', dest="sync", help='Sync file with p-values from one or more CMH test.', required=True)
parser.add_argument('--cutoff',"-c", dest="cutoff", help="Minimum count threshold - all SNPs with counts lower than this value are removed. Taken over all populaitons.", required=True)
parser.add_argument('--cov', dest="cov", help="Minimum coverage threshold on a per population base. All populations need to have a coverage above this. (default=0)", default=0)
parser.add_argument('--popcols', '-p',dest="popcols", help="Columns with allele count data to be outputted - note: column 1 = first column with count data, default=all columns")
parser.add_argument('--verbose','-v' ,dest="verb", action="store_true", help="print progress default: no output", default=False)
parser.add_argument('--header', dest="header",  action="store_true", help="Does the file contain a header? (def.: False)", default=False)
parser.add_argument("--output","-o", dest="output", help="An output file", required=True)
args = parser.parse_args()

popcols=vars(args)['popcols']
cov = int(vars(args)['cov'])
#columns with counts
if (not popcols): # if not set empty allpopcols array
	allpopcols = []
else:     # translate population columns to indeces
	allpopcols = popcols.split(",")
	allpopcols = [int(i)+2 for i in allpopcols ]

verb=vars(args)['verb']

# create output file
output_file =  vars(args)['output']
ofh = open(output_file,"w")

# parse input file and generate output
f = open(vars(args)['sync'])
cutoff=int(vars(args)['cutoff'])
cc=0
kk=0
print "Parsing file ...."
for r in f: # goes through input file 1 row at a time; r = row
	cc+=1
	r=r.rstrip()
	if ( verb and cc%1000000==0):
		print "Parsing file: processed "+str(cc)+" rows ...."
	if cc==1 and vars(args)['header']: # for first row (header in input) create new header row for output file
		print "Skipping header..."
		print >> ofh, r
		continue
	else:
		#Parse the current row and store the allelecounts in a list of dictionaries
		vals = r.split() # splits each row into columns separated by tab
		if (not allpopcols):
			allpopcols=range(3,len(vals))
		allele_counts=np.zeros(4,dtype=int)
		min_cov = True
		for i in allpopcols:
			if vals[i] != "-": 
				( ret_val, min_cov )= divide_allelecounts(vals[i],cov)
				if min_cov:
					allele_counts += ret_val
				else: # coverage too low
					min_cov = False
					break
		if not min_cov:
			next
		min_count=np.sort(allele_counts)[-2]
		#apc = [vals[i] for i in allpopcols]
		#all_count_cols = [col for col in apc if col!="-"]
		#all_pop_allelecounts = [divide_allelecounts(pop) for pop in all_count_cols]
		#alleles=get_major_and_minor_alleles(all_pop_allelecounts)		
		#Get coverages and frequencies
		# allele1=0
		# allele2=0
		# for group in all_pop_allelecounts:
		# 	allele1+=int(group[alleles[0]])
		# 	allele2+=int(group[alleles[1]])
		#min_count = alleles[0][1]
		#Print the output formatted
		if min_count >= cutoff:
			kk+=1
			print >> ofh, r
print "Job Complete: Processed "+str(cc)+" SNPs, "+str(kk)+"had minimum count >= "+ str(cutoff) + "and coverage >= " + str(cov)
