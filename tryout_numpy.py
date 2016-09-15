#!/usr/bin/env python

def create_hist(pVals,bins=1000,log=False,max_val=False):
	"""
	gets a numpy array of pValues and returns an array of bins with counts. if log=True, p values are first transformed to their negative decadic logarithm and then binned. If max_val provided binning performed from 0 to max_val.
	"""
	# gets an array of pValues and returns an array of bins with counts
	hist_count = np.zeros(bins)
	if (log):
		nlp = -1.0*np.log10(pVals)
		if (max_val):
			nlp_max=max_val
		else:
			nlp_max = nlp.max()
		pVals = pVals/nlp_max
	for i in pVals:
		j=int(bins*i)
		if j == bins:
			j -= 1
		hist_count[j] += 1
	return hist_count

def chi2distance_not_vect (obs_pval_hist,null_pval_hist):
	"""
	takes two histograms with identical bin numbers and calculates the chisquare distance between them.
	returns a float. not vectorized
	"""
	chisum = 0
	for i in range(len(obs_pval_hist)):
			if (obs_pval_hist[i] + null_pval_hist[i] > 0):
				chisum += (obs_pval_hist[i]-null_pval_hist[i])**2/(obs_pval_hist[i]+null_pval_hist[i])/2
	return chisum
def chi2distance (obs_pval_hist,null_pval_hist):
	"""
	takes two histograms with identical bin numbers and calculates the chisquare distance between them.
	returns a float. 
	"""
        chi2 = (obs_pval_hist-null_pval_hist)**2/(obs_pval_hist+null_pval_hist) * 1/2
	chisum = np.sum(chi2)
        return chisum

def kull_leib_distance_not_vect (obs_pval_hist,null_pval_hist):
	"""
	takes two histograms with identical bin numbers and calculates the average of the Kullback_Leibler divergences in both directions KL(hist1-hist2) and KL(hist2,hist1) between them. returns a float. not using the numpy funtions for vectorization
	"""
	kl_on = 0
	kl_no = 0
	tot_o=0
	tot_n=0
	#calculate total number of points for normalisation
	for i in range(len(obs_pval_hist)):
		tot_o+=obs_pval_hist[i]
		tot_n+=null_pval_hist[i]
	# go two ways as not symmetric and take the mean
	for i in range(len(obs_pval_hist)):
			if (obs_pval_hist[i] !=0 and null_pval_hist[i] != 0):
				obs_p = obs_pval_hist[i]/tot_o
				null_p = null_pval_hist[i]/tot_n
				kl_on += obs_p * np.math.log( obs_p / null_p )
				kl_no += null_p * np.math.log( null_p / obs_p )
	return (kl_on+kl_no)/2

def kull_leib_distance (obs_pval_hist,null_pval_hist):
	"""
	takes two histograms (numpy arrays) with identical bin numbers and calculates the average of the Kullback_Leibler divergences in both directions KL(hist1-hist2) and KL(hist2,hist1) between them. returns a float
	"""
	kl_on = 0
	kl_no = 0
	#calculate total number of points for normalisation
	tot_o=np.sum(obs_pval_hist)
	tot_n=np.sum(null_pval_hist)
	# normalise
        obs_p_norm = obs_pval_hist / tot_o
        null_p_norm = null_pval_hist / tot_n
        # use masked arrays
        # masked arrays allow to leave out undefined values such as nan and inf
        kl_on = np.sum( obs_p_norm * np.ma.log(obs_p_norm / null_p_norm))
        kl_no = np.sum( null_p_norm * np.ma.log(null_p_norm / obs_p_norm))
	# go two ways as not symmetric and take the mean
	return (kl_on+kl_no)/2

if __name__ == '__main__':
	import numpypy
	import sys, os
        import math
        import numpy as np
	#       from  matplotlib import pylab
        import timeit
        #create some np arrays
        a = [1.0,2.0,3.0]
        a = np.array(a)
        # fast vectorized logarithm
        b = np.log(a)
        # slow normal array transformation
        c = [ np.math.log(x) for x in a ]
        # random already gives np arrays
        obs_pval = np.random.uniform(0,1,size=10000)
        null_pval = np.random.uniform(0,1,size=10000)
        obs_hist = create_hist(obs_pval,bins=250)        
        null_hist = create_hist(null_pval,bins=250)
        chi2dist = chi2distance(obs_hist,null_hist)
        chi2dist_not_vect = chi2distance_not_vect(obs_hist,null_hist)
        kl_dist =  kull_leib_distance(obs_hist,null_hist)
        kl_dist_not_vect = kull_leib_distance_not_vect(obs_hist,null_hist)
        print "timeing chi2:"
        print 1000*(timeit.timeit(stmt="chi2distance(obs_hist,null_hist)", setup="from __main__ import chi2distance ; from __main__ import obs_hist; from __main__ import null_hist ;",number=3)) 
        print "timeing chi2_nonvect:"
        print 1000*(timeit.timeit(stmt="chi2distance_not_vect(obs_hist,null_hist)", setup="from __main__ import chi2distance_not_vect ; from __main__ import obs_hist; from __main__ import null_hist ;",number=3)) 
        # plotting a line:
        a = np.arange(0,9,0.1)
        b = np.sin(a)
        # pylab.figure()
        # pylab.plot(a,b,label="sin x")
        # pylab.legend()
        # pylab.show()
