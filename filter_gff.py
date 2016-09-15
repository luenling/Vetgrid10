import sys,re
import collections
import argparse
parser = argparse.ArgumentParser(description='removes some entry types from a gff file outputting a slightly cleaned version') 
parser.add_argument("--in", dest="infile", help="gff file to read", required=True)
parser.add_argument("--out", dest="outfile", help="gff output",required=True)
parser.add_argument("--fields", dest="ftypes", type=str, help="fields to keep (default: "genes,mIRNA")",default=5000)
parser.add_argument("--keep-same", dest="keep", action='store_true', help="keep more than one fields starting at same position apart from genes (default: "False")",default=False)
args = parser.parse_args()
infile = vars(args)['infile']
outfile = vars(args)['outfile']
winsize = vars(args)['windowsize']

