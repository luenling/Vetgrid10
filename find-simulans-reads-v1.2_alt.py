import argparse
import sys
import os
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import re
import operator

# Author: Ram Vinay Pandey
# Date: 31-03-2012
# NOTE: This script expects samtools and gsnap should be installed in main path (/usr/bin/). It also requires python 2.7.2 and biopython

parser = argparse.ArgumentParser(description="Map melanogater reads with melanogaster and simulans genome to seperate simulans reads.)")

parser = argparse.ArgumentParser(description="""
Description
-----------
    simulans reads filtering pipeline: Map melanogater reads with melanogaster and simulans genome to seperate simulans reads.""",
    formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Prerequisites
-------------
     1) python 2.7.2
     2) biopython 1.59  [biopython.org/DIST/biopython-1.59.tar.gz]


Algorithm
---------
1) Mapping all reads against combined genome (D.mel + D.sim)
2) Mapping with GSNAP done with option --split-output=gsnap_25000reads. With split-output option GSNAP will create following outputs:
    gsnap_25000reads.concordant_mult
    gsnap_25000reads.concordant_uniq
    gsnap_25000reads.halfmapping_mult
    gsnap_25000reads.halfmapping_transloc
    gsnap_25000reads.halfmapping_uniq
    gsnap_25000reads.nomapping
    gsnap_25000reads.paired_mult
    gsnap_25000reads.paired_uniq_inv
    gsnap_25000reads.paired_uniq_long
    gsnap_25000reads.paired_uniq_scr
    gsnap_25000reads.unpaired_mult
    gsnap_25000reads.unpaired_transloc
    gsnap_25000reads.unpaired_uniq

3) in "gsnap_25000reads.concordant_uniq output" (read aligned as proper pair and unique hit) reads which got hit on D.mel contig are kept in D.mel fastq file and reads which got hit on D.sim contig are kept in D.sim fastq file.
4) in "gsnap_25000reads.concordant_mult output" (read aligned as proper pair but pair got multiple hits) if mapping quality is higher when read is mapped on D.mel contig than read is mapped on D.sim contig then this read is kept in D.mel fastq file.
   If mapping quality is higher with D.sim contig alignment then read is kept in D.sim fastq file.
5) in case of ambigius hits if a read got same mapping quality with D.mel as well as with D.sim then this read always goes to D.mel fastq file.
6) Reads in rest of the output files are discarded. Neither kept in mel fastq nor in sim fastq



Disclaimer!!
------------
We strongly discourage to use any SAM output file which was created by this script for further downstream analysis. To avoid re-mapping Martin's add-ons script can be used.


Authors
-------
Ram Vinay Pandey
Christian Schloetterer
Viola Nolte



Example of a usage
------------------
To see help manual
------------------
python find-simulans-reads-v1.2.py -h

    OR
    
python2.7 find-simulans-reads-v1.2.py -h

To run on real dataset
----------------------
python find-simulans-reads-v1.2.py --mel-reference-fasta-file dmel_Phix_ref_genome.fasta --sim-reference-fasta-file Kib32_CLC_trimmed_all14+C167+wMel+wPel+wRi.fasta --read1 L3_25000_seqs_read_1.fq --read2 L3_25000_seqs_read_2.fq --fastq-type illumina --output-dir outdir --thread 20  --mean-fragment-size 400 --maximum-fragment-size 600

    OR
    
python2.7 find-simulans-reads-v1.2.py --mel-reference-fasta-file dmel_Phix_ref_genome.fasta --sim-reference-fasta-file Kib32_CLC_trimmed_all14+C167+wMel+wPel+wRi.fasta --read1 L3_25000_seqs_read_1.fq --read2 L3_25000_seqs_read_2.fq --fastq-type illumina --output-dir outdir --thread 20  --mean-fragment-size 400 --maximum-fragment-size 600


""")

parser.add_argument("--mel-reference-fasta-file", type=str, required=True, dest="melref", default=None, help="Provide the melanogaster reference sequence multi-fasta file")
parser.add_argument("--sim-reference-fasta-file", type=str, required=True, dest="simref", default=None, help="Provide the simulans reference sequence multi-fasta file")
parser.add_argument("--read1", type=str, required=True, dest="read1", default=None, help="Provide the read1 fastq file")
parser.add_argument("--read2", type=str, required=True, dest="read2", default=None, help="Provide the read2 fastq file")
parser.add_argument("--fastq-type", required=True, dest="fastqtype", type=str, help="FASTQ type illumina or sanger")
parser.add_argument("--output-dir", type=str, required=True, dest="outputdir", default=None, help="Provide the output directory path")
parser.add_argument("--thread", type=int, required=True, dest="thread", default=None, help="Provide the number of processors to be used")
#parser.add_argument("--insert-size", type=int, required=True, dest="insert_size", default=None, help="Provide the insert size")
parser.add_argument("--mean-fragment-size", type=int, required=True, dest="mean_fragment_size", default=None, help="Provide the mean fragment size read1 + gap + read2")
parser.add_argument("--maximum-fragment-size", type=int, required=True, dest="maximum_fragment_size", default=None, help="Provide the maximum fragment size read1 + gap + read2")

args = parser.parse_args()




def write_param_file(args):
    
    melref = args.melref
    simref = args.simref
    read1 = args.read1
    read2 = args.read2
    fastqtype = args.fastqtype
    outputdir = args.outputdir
    thread = args.thread
    #insert_size = args.insert_size
    mean_fragment_size = args.mean_fragment_size
    maximum_fragment_size = args.maximum_fragment_size
    
    #print melref,simref,read1,read2,fastqtype,executables_dir,outputdir,thread,insert_size
    #sys.exit()
    if not os.path.exists(outputdir):
		os.mkdir(outputdir)
    
    
	    
    paramfile = outputdir+"/find-imulans-reads.param"
    pfh=open(paramfile,"w")
    print >> pfh, "Using reference genome melanogaster reference fasta file\t"+str(melref)
    print >> pfh, "Using reference genome simulans reference fasta file\t"+str(simref)
    print >> pfh, "Using read1 file\t"+str(read1)
    print >> pfh, "Using read2 file\t"+str(read2)
    print >> pfh, "Using mean fragment size read1 + gap + read2 in GSNAP mapping\t"+str(mean_fragment_size)
    print >> pfh, "Using maximum fragment size read1 + gap + read2 in GSNAP mapping\t"+str(maximum_fragment_size)
    #print >> pfh, "Using insert sie \t"+str(insert_size)
    print >> pfh, "Using number of threads (CPU)\t"+str(thread)
    print >> pfh, "Using fastq-type\t"+str(fastqtype)
    print >> pfh, "Using output directory\t"+str(outputdir)
    pfh.close()

    arg_dict={}
    arg_dict["minmqual"] = 0
    arg_dict["outputdir"] = outputdir
    arg_dict["refseqfile"] = str(outputdir)+"/combined-genomes.fa"
    arg_dict["melref"] = melref
    arg_dict["simref"] = simref
    arg_dict["read1"] = read1
    arg_dict["read2"] = read2
    arg_dict["thread"] = int(thread)
    arg_dict["fastqtype"] = fastqtype
    arg_dict["mean_fragment_size"] = int(mean_fragment_size)
    arg_dict["maximum_fragment_size"] = int(maximum_fragment_size)
    #arg_dict["insert_size"] = int(insert_size)
    arg_dict = create_genome_index_dir(arg_dict)
    return arg_dict



def merge_two_genome(arg_dict):
    os.system("cat "+str(arg_dict["melref"])+" "+str(arg_dict["simref"])+" > "+str(arg_dict["refseqfile"]))
    
def create_genome_index_dir(arg_dict):
    gsnap_dbname = "Ref_noSNP"
    
    ### create a sub directory to write simulated reads fastq files
    reference_dir = arg_dict["outputdir"]+"/gsnap_reference"
    if not os.path.exists(reference_dir):
	os.mkdir(reference_dir)
	    
    reference_dir1 = reference_dir+"/gsnap_ref"
    if not os.path.exists(reference_dir1):
	os.mkdir(reference_dir1)
    arg_dict["gsnap_dbdir"] = reference_dir1
    arg_dict["gsnap_dbname"] = gsnap_dbname
	
    return arg_dict

def create_initial_genome_GSNAP_index(arg_dict):
    refseqfile = arg_dict["refseqfile"]
    gsnap_dbname = arg_dict["gsnap_dbname"]
    gsnap_dbdir = arg_dict["gsnap_dbdir"]
	
    #print "gmap_build -d "+gsnap_dbname+" -D "+gsnap_dbdir+" "+refseqfile
    os.system("gmap_build -d "+gsnap_dbname+" -D "+gsnap_dbdir+" "+refseqfile)


def mapping_reads(arg_dict,simchr):
    
    gsnap_dbname = arg_dict["gsnap_dbname"]
    gsnap_dbdir = arg_dict["gsnap_dbdir"]
    fastqtype = arg_dict["fastqtype"]
    thread = arg_dict["thread"]
    mapping_quality = int(arg_dict["minmqual"])
    #insert_size = int(arg_dict["insert_size"])
    mean_fragment_size = arg_dict["mean_fragment_size"]
    maximum_fragment_size = arg_dict["maximum_fragment_size"]
    output_dir = arg_dict["outputdir"]
    read1=arg_dict["read1"]
    read2=arg_dict["read2"]
    output1=""
    output2=""
    output3=""
    gsnap_command=""
    mapping_quality_filtering_command=""
    read1file = os.path.basename(read1)
    read2file = os.path.basename(read2)
    
    output1=output_dir+"/gsnap_pe.sam"
    output2=output_dir+"/"+read1file+"_"+read2file+"_gsnap_pe_mq"+str(mapping_quality)+"_pp.sam"
    output2=output_dir+"/unambi-mapped.sam"
    output3=output_dir+"/ambi-mapped.sam"   
    
    output_name=output_dir+"/gsnap_mapping_output"
    #gsnap_command = "gsnap -d "+ gsnap_dbname +" -D "+ gsnap_dbdir +" "+read1+" "+read2+" --no-sam-headers --pairmax-dna="+str(int(insert_size))+" --quality-protocol="+fastqtype+" -t "+str(int(thread))+" -A sam > "+ output1
    gsnap_command = "gsnap -d "+ gsnap_dbname +" -D "+ gsnap_dbdir +" "+read1+" "+read2+" --no-sam-headers --split-output="+str(output_name)+" --pairmax-dna="+str(int(maximum_fragment_size))+" --pairexpect="+str(int(mean_fragment_size))+" --quality-protocol="+fastqtype+" -t "+str(int(thread))+" -A sam"
   
    os.system(gsnap_command)
    
    sys.stderr.write("\n\tMESSAGE --> Filtering reads\n\n")
    all_reads = {}
    discarded_reads = {}
    for f in os.listdir(output_dir):
	file = str(output_dir)+"/"+str(f)
	if re.search(output_name,file):
	    if re.search(r'concordant_uniq$',file):
		#print f
		sim_reads_unambigious = read_sam_unambigious(file ,simchr)
		#print "unambi-proppair: "+str(len(sim_reads_unambigious))
		all_reads.update(sim_reads_unambigious)
	    elif re.search(r'concordant_mult$',file):
		#print "multi: "+str(f)
		sim_reads_ambigious = read_sam_ambigious(file,simchr)
		#print "ambi-proppair: "+str(len(sim_reads_ambigious))
		all_reads.update(sim_reads_ambigious)
	    else:
		#print "all_reads: "+str(len(all_reads))
		other_reads = read_sam_ambigious(file,simchr)
		discarded_reads.update(other_reads)
		#all_reads.update(other_reads)
		
    #discarded_reads_extra = dict([(k,all_reads[k]) for k in set(all_reads)-set(discarded_reads)])
    #in_both = dict([(k,all_reads[k]) for k in set(all_reads) & set(discarded_reads)])
    #print discarded_reads_extra
    #print len(discarded_reads_extra)
    #print in_both
    
    print "Number of read pair of simulans found: "+str(len(all_reads))
    print "Number of discarded read pairs: "+str(len(discarded_reads))
    read_fastq(arg_dict["read1"],arg_dict["read2"],all_reads,discarded_reads,output_dir)
    
    ### Deleting the intermediate SAM files
    
    for f in os.listdir(output_dir):
	file = str(output_dir)+"/"+str(f)
	if re.search(output_name,file):
	    if re.search(r'concordant_uniq$',file):
		#os.system("rm -r "+str(file))
		pass
	    elif re.search(r'concordant_mult$',file):
		#os.system("rm -r "+str(file))
		pass
	    else:
		os.system("rm -r "+str(file))
		
    
    
    
def read_fasta(input):
    genome_dict = SeqIO.to_dict(SeqIO.parse(input, "fasta"))
    genome = {}

    for chro, r in genome_dict.items():
        genome[chro] = chro
            
    return genome

def read_sam_unambigious(input,simchr):
    fh = open(input,"r")
    sim_reads={}
    for line in fh:
        if line:
            line = line.rstrip().lstrip()
            headerm = re.search(r'^@',line)
            if headerm==None:
                #print line
                a = line.split("\t")
                read_id = a[0].rstrip().lstrip()
                contig_id = a[2].rstrip().lstrip()
                
                if contig_id in simchr:
                    #print read_id,contig_id
                    sim_reads[read_id]=contig_id
                    
    
    
    fh.close()
    return sim_reads

def read_sam_ambigious(input,simchr):
    fh = open(input,"r")
    sim_reads=[]
    sim_reads_final={}
    for line in fh:
        if line:
            line = line.rstrip().lstrip()
            headerm = re.search(r'^@',line)
            if headerm==None:
                #print line
                a = line.split("\t")
                read_id = a[0].rstrip().lstrip()
                contig_id = a[2].rstrip().lstrip()
                mapping_quality = int(a[4].rstrip().lstrip())
                if len(sim_reads) == 0:
                    #initialize sim_read array
                    sim_reads = [read_id,mapping_quality,set([contig_id ])]
                if sim_reads[0] == read_id:
                    # same read id -> check mapping quality 
                    if (sim_reads[1] < mapping_quality):
                        sim_reads[1] = mapping_quality
                        sim_reads[2] = set([ contig_id ])
                    elif (sim_reads[1] == mapping_quality):
                        sim_reads[2].add( contig_id )
                else:
                    # a new read id coming up
                    # check whether previous read id mapped to melanogaster
                    melcount = 0
                    for contig in sim_reads[2]:
                        if not (contig in simchr):
                            # not a simulans chromosome
                            melcount += 1
                    if (melcount == 0):
                        sim_reads_final[sim_reads[0]]=sim_reads_final.setdefault(sim_reads[0],sim_reads[1])
                        if (sim_reads_final[sim_reads[0]] < sim_reads[1]):
                            sim_reads_final[sim_reads[0]] == sim_reads[1]
                    else:
                        #if there exists a previous entry in sim_reads for read id and its quality is equal or less than the current one, remove it
                        if (sim_reads[0] in sim_reads_final) and (sim_reads_final[sim_reads[0]] <= sim_reads[1]):
                            del sim_reads_final[sim_reads[0]]
                    sim_reads=[read_id,mapping_quality,set([contig_id ])]
    fh.close()
    return sim_reads_final


def read_fastq(read1,read2,sim_reads,discarded_reads,outputdir):
    fh1=open(read1,"r")
    fh2=open(read2,"r")
    
    read1file = os.path.basename(read1)
    read2file = os.path.basename(read2)
    
    mel_read1 = outputdir+"/"+str(read1file)+"_mel_read_1.fq"
    mel_read2 = outputdir+"/"+str(read2file)+"_mel_read_2.fq"
    
    sim_read1 = outputdir+"/"+str(read1file)+"_sim_read_1.fq"
    sim_read2 = outputdir+"/"+str(read2file)+"_sim_read_2.fq"
    
    #discarded_read1 = outputdir+"/"+str(read1file)+"_discarded_read_1.fq"
    #discarded_read2 = outputdir+"/"+str(read2file)+"_discarded_read_2.fq"
    
    melofh1 = open(mel_read1,"w")
    melofh2 = open(mel_read2,"w")
    
    simofh1 = open(sim_read1,"w")
    simofh2 = open(sim_read2,"w")
    
    #disofh1 = open(discarded_read1,"w")
    #disofh2 = open(discarded_read1,"w")
    
    
    while(True):
        
        h1 = fh1.readline().rstrip()
        s1 = fh1.readline().rstrip()
        he1 = fh1.readline().rstrip()
        qual1 = fh1.readline().rstrip()
            
        h2 = fh2.readline().rstrip()
        s2 = fh2.readline().rstrip()
        he2 = fh2.readline().rstrip()
        qual2 = fh2.readline().rstrip()
        if (h1=="" or h2==""):
            break
        if h1[1:-2] in discarded_reads:
	    #print >> disofh1, h1
	    #print >> disofh1, s1
	    #print >> disofh1, he1
	    #print >> disofh1, qual1
		
	    #print >> disofh2, h2
	    #print >> disofh2, s2
	    #print >> disofh2, he2
	    #print >> disofh2, qual2
	    pass
	
	else:
	    if h1[1:-2] in sim_reads:
    
		print >> simofh1, h1
		print >> simofh1, s1
		print >> simofh1, he1
		print >> simofh1, qual1
		
		print >> simofh2, h2
		print >> simofh2, s2
		print >> simofh2, he2
		print >> simofh2, qual2
		
		
	    else:
		print >> melofh1, h1
		print >> melofh1, s1
		print >> melofh1, he1
		print >> melofh1, qual1
		
		print >> melofh2, h2
		print >> melofh2, s2
		print >> melofh2, he2
		print >> melofh2, qual2
            
            
    fh1.close()
    fh2.close()
    melofh1.close()
    melofh2.close()
    simofh1.close()
    simofh2.close()
    #disofh1.close()
    #disofh2.close()
    sim_reads.clear()
    discarded_reads.clear()
    
arg_dict = write_param_file(args)

### Merging two genomes fasta file
merge_two_genome(arg_dict)

### Creating GSNAP index for mappin
create_initial_genome_GSNAP_index(arg_dict)


output_dir = arg_dict["outputdir"]

### Reading Simulans Fasta files to get a dictonary of Simulans contigs
simchr=read_fasta(arg_dict["simref"])

sys.stderr.write("\n\tMESSAGE --> Mapping reads to combined genomes\n\n")
mapping_reads(arg_dict,simchr)

sys.exit()

