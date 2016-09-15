import sys, re
import os
import argparse
import select
import gzip

parser = argparse.ArgumentParser(description="""Reads replacement list and replaces all occurences of term in first or other given columns in none (framed by whitespace or word boundaries) in a file with the other. Usually replaces column 1 with 2, can be changed.  
""")

parser.add_argument("--in","-i", dest="infile", help="text file", required=True)
parser.add_argument("--rep","-r", dest="replacements", help="replacement file", required=True)
parser.add_argument("--dir","-d", dest="direction", help="which column (one based) hold terms to be replaced and replacement (default: 1:2)", default="1:2")
args = parser.parse_args()
infile = vars(args)['infile']
replacements = vars(args)['replacements']
direction = [ int(x) - 1 for x in vars(args)['direction'].split(":")]
# load replacement file
reps=[]
with open(replacements,'r') as inf:
    for line in inf:
        if re.match(".*\#",line):
            continue
        line=line.rstrip()
        fields=line.split()
        reps.add([fields[direction[0]],fields[direction[0]] ])

# replace stuff:
# create replacement regexp:
repreg = re.compile("|".join('((?<![.:\w])%s(?![.:\w]))' % re.escape(x[0]) for x in reps) )
if re.search("\.b?gz",infile):
    inf = gzip.open(infile,'rb')
else:
    inf = open(infile,"r")

for line in inf:
    if re.match(".*\#",line):
        print line
        continue
    line=line.rstrip()
    line=re.sub(repreg,lambda m: reps[m.lastindex - 1][1],line)
    print line

