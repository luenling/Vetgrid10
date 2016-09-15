from collections import defaultdict
import re


trdict=defaultdict(dict)
for idx,line in enumerate(a):
    if re.match("\#",line):
        continue
    fields=line.split("\t")
    len=abs(int(fields[3])-int(fields[4]))+1
    regres=re.search("transcript_id\s+\"(\w+)\".*exon_number\s+\"(\d+)\"",fields[-1])
    if not regres:
        continue
    trid=regres.group(1)
    exnmbr=int(regres.group(2))
    trdict[trid][exnmbr]=[len,idx]
