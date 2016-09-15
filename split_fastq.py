# -*- coding: utf-8 -*-
"""
Created on Mon May 18 10:59:51 2015

@author: lukasendler
"""

def create_fh(base,barcode,infile,gzout):
    # get rid of (b)gz
    fn=re.sub("\.b?gz$","",os.path.basename(infile))
    fn=base+"_"+barcode+"_"+fn
    if gzout:
        fn+=".gz"
        return gzip.open(fn,mode="wb")
    return open(fn,mode="wb")


import re,os
import argparse
import gzip
#from collections import defaultdict
from Bio.SeqIO.QualityIO import FastqGeneralIterator
 
parser = argparse.ArgumentParser(description="""reads a fastq file and splits by barcodes.  
resulting files are named prefix_barcode_original_name
reads (b)gzipped files andby default produces gzipped files
""")

parser.add_argument("--files","-f", dest="infile", help="fastq filenames separated by commas", required=True)
parser.add_argument("--out","-o", dest="base", help="outfile prefix (default: split )", default="split")
parser.add_argument("--gzout","-g", dest="gzout", action='store_false', help="outfile basename (default: True )", default=True)
args = parser.parse_args()
infiles = vars(args)['infile'].split(",")
gzout = vars(args)['gzout']
base = vars(args)['base']

# regular expression for catching barcodes
repreg = re.compile("\#+([ACGTN]+)(?:[_/120]*)$")
for infile in infiles:
    if re.search("\.b?gz",infile):
        inf = gzip.open(infile,'rb')
    else:
        inf = open(infile,"r")

    #barcodes=set()
    fh_dict={}
    print "Working through file"+infile
    # read fastq entries
    for name,seq,qual in FastqGeneralIterator(inf):
        # find barcode
        query=re.search(repreg,name)
        if not(query):
            continue
        barcode=query.groups()[0]
        if not(barcode in fh_dict.keys()):
            # create file and filehandle
            fh_dict[barcode]=create_fh(base,barcode,infile,gzout)
        print >> fh_dict[barcode], "@"+name+"\n"+seq+"\n+\n"+qual 
    for barcode in fh_dict:
        fh_dict[barcode].close()
 
