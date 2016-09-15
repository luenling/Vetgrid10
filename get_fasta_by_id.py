import Bio
from Bio import Entrez
import sys
for i in sys.argv[1:]:
    Entrez.email = "lukasendler@gmail.com"
    handle = Entrez.efetch(db="nucleotide", id=i, rettype="fasta", retmode="text")
    fh=open(i+".fa","w")
    print >> fh, handle.read()
