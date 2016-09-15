import sys, os, re 
import Sync_parser
import argparse
parser = argparse.ArgumentParser(description='takes a sync file with strand information and writes a sync file with different strand bias estimates:\n either fisher pV for all populations, the total pV and the ratios of the less common to the more common strand for major and minor alleles; \n with --total it adds only the fisher pV of the total strand sums to the line') 
parser.add_argument("--in", "-i",dest="infile", help="sync file to be read (required)", required=True)
parser.add_argument("--out","-o", dest="outfile", help="sync file with bias to be written (optional)", default=False)
parser.add_argument("--total","-t", dest="total", action="store_true", help="print only the total fisher pV at the end of the full line", default=False)
parser.add_argument("--header", dest="header", action="store_true", help="print header (not for total)", default=False)
args = parser.parse_args()
infile = vars(args)['infile']
outfile= vars(args)['outfile']
header= vars(args)['header']
total= vars(args)['total']
if total:
    header=False
if not outfile:
    outfile=infile+".strandbias"
inf=open(infile,"r")
outf=open(outfile,"w")

for line in inf:
    if line[0] == '#':
        continue
    syncline = Sync_parser.SyncLineParser(line,strands=True)
    biases=syncline.get_strand_bias()
    if biases == False:
        print >> sys.stderr, "Problem calculating strand bias at:\n"+line
        continue
    if header:
        head_str="#CHR\tBPS\tREF\t"
        head_str+="\t".join([ "fisherP_pop_"+str(x) for x in range(1,len(syncline.seqs)+1) ])
        head_str+="\tfisherP_total\tfrac_maj\tfrac_min"
        print >>outf,head_str
        header=False
    if total:
        pop_pv = "\t"+"{:.3g}".format(biases[1])
        line=line.rstrip()
    else:
        line="\t".join(line.split("\t")[:3])
        pop_pv = "\t".join([ "{:.3g}".format(x) for x in biases[0] ])
        pop_pv += "\t"+"{:.3g}".format(biases[1])
        pop_pv += "\t"+ "\t".join([ "{:.3g}".format(x) for x in biases[2] ])
    print >>outf,line + "\t"+ pop_pv
    
