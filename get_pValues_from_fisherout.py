#!/usr/bin/env python

import sys, os, glob, re, math

import argparse
parser = argparse.ArgumentParser(description='calculate histogram distances between a file with observed PValues and a group of files each containing the same number of P values') 
parser.add_argument("--in", dest="infile", help="fisher file with P values", required=True)
parser.add_argument("--out", dest="outfile", help="output file",default="")
parser.add_argument("--comp", dest="comps", help="comparisons to retrieve (eg: \"1:4,2:5,3:6\")",required=True)
args = parser.parse_args()
comps = vars(args)['comps'].split(",")
infile = vars(args)['infile']
outfile = vars(args)['outfile']
if (outfile == ""):
    outfile = re.sub(r'\.\w+$','',infile)+".fishcomp"
inf = open(infile,"r")
f_pos={}
out = open(outfile,"w")
for line in inf:
    line.rstrip()
    fields=line.split()
    if (len(f_pos) == 0):
        for i in range(0,len(fields)):
            for j in comps:
                if (re.match(j+"=",fields[i])):
                    f_pos[j]=i
        diff_xor = set(f_pos) ^ set(comps)
        if (len(diff_xor) > 0):
            sys.exit("Field"+str(diff_xor)+" not found")
        print >>out,"\#"+"\t".join(comps)
        pos=[ f_pos[x] for x in comps]
    out_array=[ fields[x].split("=")[1] for x in pos]
    for i in range(0,len(out_array)):
        if (out_array[i] != "na"):
            out_array[i] = str(10**(-1*float( out_array[i])))
    print >>out,"\t".join(out_array)
inf.close()
out.close()
