#----------------------------------------------------------
##Part 1: Map and extract new F15 replicate data
#----------------------------------------------------------

##Fasta files from BGI
#/Volumes/vetgrid12/Downloads/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_1.fq
#/Volumes/vetgrid12/Downloads/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_2.fq

##Key for renamed files
# index_4	TGACCAAT = Pablo´s D. mel cage  HOT F15  r. 2 females
# index_5	ACAGTGAT = Pablo´s D. mel cage  HOT F15  r. 3 females
# index_6	GCCAATAT = Pablo´s D. mel cage  COLD F15  r. 6 females
# index_7	CAGATCAT = Pablo´s D. mel cage  COLD F15  r. 7 females

#Everything mapped on vetgrid11 unless otherwise noted
#output files being temporarily stored in:
#/Volumes/disk5s1/F15HotCold_new_reps 
#samtools 0.1.9
#bwa 0.5.8c
#popoolation v. 208
#popoolation2 v. 175

##Trim combined reads
perl /Volumes/Temp/popoolation/basic-pipeline/trim-fastq.pl --input1 /Volumes/vetgrid12/Downloads/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_1.fq --input2 /Volumes/vetgrid12/Downloads/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_2.fq --output /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed --quality-threshold 18 --fastq Illumina --min-length 50 > /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed.stats

##Map combined reads
#Vetgrid 11
bwa aln -t 20 -l 200 -n 0.01 -e 12 -d 12 -o 2 /Volumes/Temp/Ray/new-ref/dmelwolb_UminusMito.fasta /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed_1 > /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed_1.sai

#Vetgrid 12
bwa aln -t 20 -l 200 -n 0.01 -e 12 -d 12 -o 2 /Volumes/Temp/Ray/new-ref/dmelwolb_UminusMito.fasta /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed_2 > /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed_2.sai

#Run sampe
bwa sampe /Volumes/Temp/Ray/new-ref/dmelwolb_UminusMito.fasta /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed_1.sai /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed_2.sai /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed_1 /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed_2 > /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed.sam

#Make bam
samtools view -q 20 -S -h -b -F 0x0008 -T /Volumes/Temp/Ray/new-ref/dmelwolb_UminusMito.fasta /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed.sam | samtools view -hb -f 0x0002 - | samtools view -hb -F 0x0004 - > /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed_PropPairs.bam

#Get header
samtools view -H /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed_PropPairs.bam > /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_inh.sam

#Extract population reads
samtools view /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed_PropPairs.bam | grep '#TGACCAAT' | cat /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_inh.sam - | samtools view -Shb - > /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_HotF15r2_PropPairs.bam &
samtools view /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed_PropPairs.bam | grep '#ACAGTGAT' | cat /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_inh.sam - | samtools view -Shb - > /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_HotF15r3_PropPairs.bam &
samtools view /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed_PropPairs.bam | grep '#GCCAATAT' | cat /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_inh.sam - | samtools view -Shb - > /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_ColdF15r6_PropPairs.bam &
samtools view /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_trimmed_PropPairs.bam | grep '#CAGATCAT' | cat /Volumes/disk5s1/F15HotCold_new_reps/BGI_32a_120619_I321_FCC0YELACXX_L2_CHKPEI12060003_inh.sam - | samtools view -Shb - > /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_ColdF15r7_PropPairs.bam
 
wait

samtools sort /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_HotF15r2_PropPairs.bam /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_HotF15r2_PropPairs_sorted2 &
samtools sort /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_HotF15r3_PropPairs.bam /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_HotF15r3_PropPairs_sorted2 &
samtools sort /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_ColdF15r6_PropPairs.bam /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_ColdF15r6_PropPairs_sorted2 &
samtools sort /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_ColdF15r7_PropPairs.bam /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_ColdF15r7_PropPairs_sorted2

#--------------------------------------------------------------
##Part 2: Remove extra contigs and reheader F33 and F59 data
#--------------------------------------------------------------

samtools view -h /Volumes/Temp-4/martin/Projects/exp-ev/novel_timepoints/data/F33R8_sort.bam | awk '!($3=="Acetobacter" || $3=="Lactobacillus" || $3=="gi|9626372|ref|NC_001422.1|")' | samtools view -Sbh - > /Volumes/disk5s1/BAMs/F33R8_sort_RemovedExtraGenomes.bam &
samtools view -h /Volumes/Temp-4/martin/Projects/exp-ev/novel_timepoints/data/F33R9_sort.bam | awk '!($3=="Acetobacter" || $3=="Lactobacillus" || $3=="gi|9626372|ref|NC_001422.1|")' | samtools view -Sbh - > /Volumes/disk5s1/BAMs/F33R9_sort_RemovedExtraGenomes.bam &
samtools view -h /Volumes/Temp-4/martin/Projects/exp-ev/novel_timepoints/data/F33R10_sort.bam | awk '!($3=="Acetobacter" || $3=="Lactobacillus" || $3=="gi|9626372|ref|NC_001422.1|")' | samtools view -Sbh - > /Volumes/disk5s1/BAMs/F33R10_sort_RemovedExtraGenomes.bam &
samtools view -h /Volumes/Temp-4/martin/Projects/exp-ev/novel_timepoints/data/F59R1_sort.bam | awk '!($3=="Acetobacter" || $3=="Lactobacillus" || $3=="gi|9626372|ref|NC_001422.1|")' | samtools view -Sbh - > /Volumes/disk5s1/BAMs/F59R1_sort_RemovedExtraGenomes.bam &
samtools view -h /Volumes/Temp-4/martin/Projects/exp-ev/novel_timepoints/data/F59R4_sort.bam | awk '!($3=="Acetobacter" || $3=="Lactobacillus" || $3=="gi|9626372|ref|NC_001422.1|")' | samtools view -Sbh - > /Volumes/disk5s1/BAMs/F59R4_sort_RemovedExtraGenomes.bam &
samtools view -h /Volumes/Temp-4/martin/Projects/exp-ev/novel_timepoints/data/F59R5_sort.bam | awk '!($3=="Acetobacter" || $3=="Lactobacillus" || $3=="gi|9626372|ref|NC_001422.1|")' | samtools view -Sbh - > /Volumes/disk5s1/BAMs/F59R5_sort_RemovedExtraGenomes.bam

wait

samtools sort /Volumes/disk5s1/BAMs/F33R8_sort_RemovedExtraGenomes.bam /Volumes/disk5s1/BAMs/F33R8_sort2_RemovedExtraGenomes &
samtools sort /Volumes/disk5s1/BAMs/F33R9_sort_RemovedExtraGenomes.bam /Volumes/disk5s1/BAMs/F33R9_sort2_RemovedExtraGenomes &
samtools sort /Volumes/disk5s1/BAMs/F33R10_sort_RemovedExtraGenomes.bam /Volumes/disk5s1/BAMs/F33R10_sort2_RemovedExtraGenomes &
samtools sort /Volumes/disk5s1/BAMs/F59R1_sort_RemovedExtraGenomes.bam /Volumes/disk5s1/BAMs/F59R1_sort2_RemovedExtraGenomes &
samtools sort /Volumes/disk5s1/BAMs/F59R4_sort_RemovedExtraGenomes.bam /Volumes/disk5s1/BAMs/F59R4_sort2_RemovedExtraGenomes &
samtools sort /Volumes/disk5s1/BAMs/F59R5_sort_RemovedExtraGenomes.bam /Volumes/disk5s1/BAMs/F59R5_sort2_RemovedExtraGenomes

wait

#replace header in F33 and F59 bams
samtools reheader /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/F33R8_sort2_RemovedExtraGenomes.bam > /Volumes/disk5s1/BAMs/F33R8_sort2_RemovedExtraGenomes_NewHeader.bam &
samtools reheader /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/F33R9_sort2_RemovedExtraGenomes.bam > /Volumes/disk5s1/BAMs/F33R9_sort2_RemovedExtraGenomes_NewHeader.bam &
samtools reheader /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/F33R10_sort2_RemovedExtraGenomes.bam > /Volumes/disk5s1/BAMs/F33R10_sort2_RemovedExtraGenomes_NewHeader.bam &
samtools reheader /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/F59R1_sort2_RemovedExtraGenomes.bam > /Volumes/disk5s1/BAMs/F59R1_sort2_RemovedExtraGenomes_NewHeader.bam &
samtools reheader /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/F59R4_sort2_RemovedExtraGenomes.bam > /Volumes/disk5s1/BAMs/F59R4_sort2_RemovedExtraGenomes_NewHeader.bam &
samtools reheader /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/F59R5_sort2_RemovedExtraGenomes.bam > /Volumes/disk5s1/BAMs/F59R5_sort2_RemovedExtraGenomes_NewHeader.bam 

#--------------------------------------------------------------
##Part 3: Make BAMs from Pablo mapped populations
#--------------------------------------------------------------

##Remake filtered BAMs
#Date 03.07.12

####All original bam files stored in blackpupa
####Files temporarily saved on portable disk in vetgrid11

###merge pops with different sequence runs
samtools merge -h /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/BR3_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/BR3_run1_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/BR3_run2_bwa058c_sorted.bam

samtools merge -h /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/BR9_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/BR9_s1_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/BR9_s5_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/BR9_s6_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/BR9_s7_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/BR9_s8_bwa058c_sorted.bam

samtools merge -h /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/F15r4_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F15r4_s1_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F15r4_s8_run1_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F15r4_s8_run2_bwa058c_sorted.bam

samtools merge -h /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/F15r5_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F15r5_chicken_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F15r5_low_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F15r5_s6_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F15r5_s7_bwa058c_sorted.bam

samtools merge -h /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/F15r8_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F15r8_s1_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F15r8_s2_bwa058c_sorted.bam

samtools merge -h /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/F15r9_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F15r9_s3_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F15r9_s4_bwa058c_sorted.bam

samtools merge -h /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/F15r10_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F15r10_s1_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F15r10_s2_bwa058c_sorted.bam

samtools merge -h /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/F23r1_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F23r1_s1_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F23r1_s6_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F23r1_s7_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F23r1_s8_bwa058c_sorted.bam

samtools merge -h /Volumes/Temp/Ray/BR9_inh.sam /Volumes/disk5s1/BAMs/F27r5_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F27r5_s5_bwa058c_sorted.bam /Volumes/Temp-1/Ray/Dmel_pipeline_all_unfiltered_BAMs/F27r5_s6_bwa058c_sorted.bam

###remove ambiguous reads
samtools view -q 20 -hb -F 0x0008 /Volumes/disk5s1/BAMs/BR3_bwa058c_sorted.bam | samtools view -hb -F 0x0004 - | samtools view -hb -f 0x0002 - | samtools view -hb - > /Volumes/disk5s1/BAMs/BR3_bwa058c_sorted_PropPairs.bam &

samtools view -q 20 -hb -F 0x0008 /Volumes/disk5s1/Dmel_pipeline_all_unfiltered_BAMs/BR4_bwa058c_NEW_sorted.bam | samtools view -hb -F 0x0004 - | samtools view -hb -f 0x0002 - | samtools view -hb - > /Volumes/disk5s1/BAMs/BR4_bwa058c_sorted_PropPairs.bam &

samtools view -q 20 -hb -F 0x0008 /Volumes/disk5s1/BAMs/BR9_bwa058c_sorted.bam | samtools view -hb -F 0x0004 - | samtools view -hb -f 0x0002 - | samtools view -hb - > /Volumes/disk5s1/BAMs/BR9_bwa058c_sorted_PropPairs.bam &

samtools view -q 20 -hb -F 0x0008 /Volumes/disk5s1/BAMs/F15r4_bwa058c_sorted.bam | samtools view -hb -F 0x0004 - | samtools view -hb -f 0x0002 - | samtools view -hb - > /Volumes/disk5s1/BAMs/F15r4_bwa058c_sorted_PropPairs.bam &

samtools view -q 20 -hb -F 0x0008 /Volumes/disk5s1/BAMs/F15r5_bwa058c_sorted.bam | samtools view -hb -F 0x0004 - | samtools view -hb -f 0x0002 - | samtools view -hb - > /Volumes/disk5s1/BAMs/F15r5_bwa058c_sorted_PropPairs.bam &

samtools view -q 20 -hb -F 0x0008 /Volumes/disk5s1/BAMs/F15r8_bwa058c_sorted.bam | samtools view -hb -F 0x0004 - | samtools view -hb -f 0x0002 - | samtools view -hb - > /Volumes/disk5s1/BAMs/F15r8_bwa058c_sorted_PropPairs.bam &

samtools view -q 20 -hb -F 0x0008 /Volumes/disk5s1/BAMs/F15r9_bwa058c_sorted.bam | samtools view -hb -F 0x0004 - | samtools view -hb -f 0x0002 - | samtools view -hb - > /Volumes/disk5s1/BAMs/F15r9_bwa058c_sorted_PropPairs.bam &

samtools view -q 20 -hb -F 0x0008 /Volumes/disk5s1/BAMs/F15r10_bwa058c_sorted.bam | samtools view -hb -F 0x0004 - | samtools view -hb -f 0x0002 - | samtools view -hb - > /Volumes/disk5s1/BAMs/F15r10_bwa058c_sorted_PropPairs.bam &

samtools view -q 20 -hb -F 0x0008 /Volumes/disk5s1/BAMs/F23r1_bwa058c_sorted.bam | samtools view -hb -F 0x0004 - | samtools view -hb -f 0x0002 - | samtools view -hb - > /Volumes/disk5s1/BAMs/F23r1_bwa058c_sorted_PropPairs.bam &

samtools view -q 20 -hb -F 0x0008 /Volumes/disk5s1/BAMs/F27r5_bwa058c_sorted.bam | samtools view -hb -F 0x0004 - | samtools view -hb -f 0x0002 - | samtools view -hb - > /Volumes/disk5s1/BAMs/F27r5_bwa058c_sorted_PropPairs.bam &

samtools view -q 20 -hb -F 0x0008 /Volumes/disk5s1/Dmel_pipeline_all_unfiltered_BAMs/F37r1_bwa058c_sorted.bam | samtools view -hb -F 0x0004 - | samtools view -hb -f 0x0002 - | samtools view -hb - > /Volumes/disk5s1/BAMs/F37r1_bwa058c_sorted_PropPairs.bam &

samtools view -q 20 -hb -F 0x0008 /Volumes/disk5s1/Dmel_pipeline_all_unfiltered_BAMs/F37r4_bwa058c_sorted.bam | samtools view -hb -F 0x0004 - | samtools view -hb -f 0x0002 - | samtools view -hb - > /Volumes/disk5s1/BAMs/F37r4_bwa058c_sorted_PropPairs.bam &

samtools view -q 20 -hb -F 0x0008 /Volumes/disk5s1/Dmel_pipeline_all_unfiltered_BAMs/F37r5_bwa058c_sorted.bam | samtools view -hb -F 0x0004 - | samtools view -hb -f 0x0002 - | samtools view -hb - > /Volumes/disk5s1/BAMs/F37r5_bwa058c_sorted_PropPairs.bam

wait

###sort bam
samtools sort /Volumes/disk5s1/BAMs/BR3_bwa058c_sorted_PropPairs.bam /Volumes/disk5s1/BAMs/BR3_bwa058c_sorted2_PropPairs &

samtools sort /Volumes/disk5s1/BAMs/BR4_bwa058c_sorted_PropPairs.bam /Volumes/disk5s1/BAMs/BR4_bwa058c_sorted2_PropPairs &

samtools sort /Volumes/disk5s1/BAMs/BR9_bwa058c_sorted_PropPairs.bam /Volumes/disk5s1/BAMs/BR9_bwa058c_sorted2_PropPairs &

samtools sort /Volumes/disk5s1/BAMs/F15r4_bwa058c_sorted_PropPairs.bam /Volumes/disk5s1/BAMs/F15r4_bwa058c_sorted2_PropPairs &

samtools sort /Volumes/disk5s1/BAMs/F15r5_bwa058c_sorted_PropPairs.bam /Volumes/disk5s1/BAMs/F15r5_bwa058c_sorted2_PropPairs &

samtools sort /Volumes/disk5s1/BAMs/F15r8_bwa058c_sorted_PropPairs.bam /Volumes/disk5s1/BAMs/F15r8_bwa058c_sorted2_PropPairs &

samtools sort /Volumes/disk5s1/BAMs/F15r9_bwa058c_sorted_PropPairs.bam /Volumes/disk5s1/BAMs/F15r9_bwa058c_sorted2_PropPairs &

samtools sort /Volumes/disk5s1/BAMs/F15r10_bwa058c_sorted_PropPairs.bam /Volumes/disk5s1/BAMs/F15r10_bwa058c_sorted2_PropPairs &

samtools sort /Volumes/disk5s1/BAMs/F23r1_bwa058c_sorted_PropPairs.bam /Volumes/disk5s1/BAMs/F23r1_bwa058c_sorted2_PropPairs &

samtools sort /Volumes/disk5s1/BAMs/F27r5_bwa058c_sorted_PropPairs.bam /Volumes/disk5s1/BAMs/F27r5_bwa058c_sorted2_PropPairs &

samtools sort /Volumes/disk5s1/BAMs/F37r1_bwa058c_sorted_PropPairs.bam /Volumes/disk5s1/BAMs/F37r1_bwa058c_sorted2_PropPairs &

samtools sort /Volumes/disk5s1/BAMs/F37r4_bwa058c_sorted_PropPairs.bam /Volumes/disk5s1/BAMs/F37r4_bwa058c_sorted2_PropPairs &

samtools sort /Volumes/disk5s1/BAMs/F37r5_bwa058c_sorted_PropPairs.bam /Volumes/disk5s1/BAMs/F37r5_bwa058c_sorted2_PropPairs

#--------------------------------------------------------------
##Part 4: Make pileup, mask pileup, make sync and cmh files
#--------------------------------------------------------------


#all BAMs except new replicates now in /BAMs folder
cd /Volumes/disk5s1/BAMs/

#Make mpileup 
samtools mpileup -Bf /Volumes/Temp/Ray/new-ref/dmelwolb_UminusMito.fasta BR3_bwa058c_sorted2_PropPairs.bam BR4_bwa058c_sorted2_PropPairs.bam BR9_bwa058c_sorted2_PropPairs.bam /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_HotF15r2_PropPairs_sorted.bam /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_HotF15r3_PropPairs_sorted.bam F15r4_bwa058c_sorted2_PropPairs.bam F15r5_bwa058c_sorted2_PropPairs.bam /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_ColdF15r6_PropPairs_sorted.bam /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_ColdF15r7_PropPairs_sorted.bam F15r8_bwa058c_sorted2_PropPairs.bam F15r9_bwa058c_sorted2_PropPairs.bam F15r10_bwa058c_sorted2_PropPairs.bam F23r1_bwa058c_sorted2_PropPairs.bam F27r5_bwa058c_sorted2_PropPairs.bam F33R8_sort2_RemovedExtraGenomes_NewHeader.bam F33R9_sort2_RemovedExtraGenomes_NewHeader.bam F33R10_sort2_RemovedExtraGenomes_NewHeader.bam F37r1_bwa058c_sorted2_PropPairs.bam F37r4_bwa058c_sorted2_PropPairs.bam F37r5_bwa058c_sorted2_PropPairs.bam F59R1_sort2_RemovedExtraGenomes_NewHeader.bam F59R4_sort2_RemovedExtraGenomes_NewHeader.bam F59R5_sort2_RemovedExtraGenomes_NewHeader.bam > /Volumes/disk5s1/AllDmelCagePops_F0toF59.mpileup &


#Make pileup (first merge all bams)
samtools merge -h /Volumes/Temp/Ray/BR9_inh.sam AllDmelCagePops_F0toF59.bam BR3_bwa058c_sorted2_PropPairs.bam BR4_bwa058c_sorted2_PropPairs.bam BR9_bwa058c_sorted2_PropPairs.bam /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_HotF15r2_PropPairs_sorted.bam /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_HotF15r3_PropPairs_sorted.bam F15r4_bwa058c_sorted2_PropPairs.bam F15r5_bwa058c_sorted2_PropPairs.bam /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_ColdF15r6_PropPairs_sorted.bam /Volumes/disk5s1/F15HotCold_new_reps/BGI_new_reps_ColdF15r7_PropPairs_sorted.bam F15r8_bwa058c_sorted2_PropPairs.bam F15r9_bwa058c_sorted2_PropPairs.bam F15r10_bwa058c_sorted2_PropPairs.bam F23r1_bwa058c_sorted2_PropPairs.bam F27r5_bwa058c_sorted2_PropPairs.bam F33R8_sort2_RemovedExtraGenomes_NewHeader.bam F33R9_sort2_RemovedExtraGenomes_NewHeader.bam F33R10_sort2_RemovedExtraGenomes_NewHeader.bam F37r1_bwa058c_sorted2_PropPairs.bam F37r4_bwa058c_sorted2_PropPairs.bam F37r5_bwa058c_sorted2_PropPairs.bam F59R1_sort2_RemovedExtraGenomes_NewHeader.bam F59R4_sort2_RemovedExtraGenomes_NewHeader.bam F59R5_sort2_RemovedExtraGenomes_NewHeader.bam

wait

samtools pileup -Bf /Volumes/Temp/Ray/new-ref/dmelwolb_UminusMito.fasta AllDmelCagePops_F0toF60.bam > /Volumes/disk5s1/AllDmelCagePops_F0toF60.pileup

#Identify indels
perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/indel_filtering/identify-indel-regions.pl --min-count 2 --indel-window 5 --input /Volumes/disk5s1/AllDmelCagePops_F0toF59.pileup --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_pileup_INDELMASK_mc2.gtf &

perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/indel_filtering/identify-indel-regions.pl --min-count 2 --indel-window 5 --input /Volumes/disk5s1/AllDmelCagePops_F0toF59.mpileup --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_mpileup_INDELMASK_mc2.gtf &

perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/indel_filtering/identify-indel-regions.pl --min-count 5 --indel-window 5 --input /Volumes/disk5s1/AllDmelCagePops_F0toF59.mpileup --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_mpileup_INDELMASK_mc5.gtf &

perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/indel_filtering/identify-indel-regions.pl --min-count 10 --indel-window 5 --input /Volumes/disk5s1/AllDmelCagePops_F0toF59.mpileup --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_mpileup_INDELMASK_mc10.gtf &

perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/indel_filtering/identify-indel-regions.pl --min-count 15 --indel-window 5 --input /Volumes/disk5s1/AllDmelCagePops_F0toF59.mpileup --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_mpileup_INDELMASK_mc15.gtf

#make combined indel and repeat gtf for masking
sed 's/\;/ /g' /Volumes/Temp/Ray/Masking_files/Repeat_gff/dmel5.38-clean_repeat-masked_SSR-NOT-masked.fa.out.gff | sed 's/\ID=/ /g' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,"gene_id \""$9"\";"}' > /Volumes/Temp/Ray/Masking_files/Repeat_gff/dmel5.38-clean_repeat-masked_SSR-NOT-masked.fa.out.gtf

cat /Volumes/Temp/Ray/Masking_files/Repeat_gff/dmel5.38-clean_repeat-masked_SSR-NOT-masked.fa.out.gtf /Volumes/Temp/Ray/Masking_files/Indel_gtf/AllDmelCagePops_F0toF59_mpileup_INDELMASK_mc10.gtf > /Volumes/Temp/Ray/Masking_files/Merged_repeat_indel.gtf

#Mask indels and repeats
perl /Volumes/Temp/Ray/software_ray/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59.mpileup --gtf /Volumes/Temp/Ray/Masking_files/Merged_repeat_indel.gtf --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.mpileup
 
#Make synchronized file
java -Xmx20g -jar /Volumes/Temp/Ray/software_ray/popoolation2-svn/mpileup2sync.jar --threads 20 --min-qual 20 --fastq-type illumina --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.mpileup --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync

#Run CMH tests
#run on vetgrid11
#1. BaseHot
perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 1-13,2-6,3-7 --remove-temp &

perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-5_3-4.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 1-5,3-4 --remove-temp &

#2. BaseCold
perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 1-10,2-12,3-11 --remove-temp &

perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-8_3-9.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 1-6,3-9 --remove-temp &

#3. HotCold
perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15HotF15Cold_13-10_6-12_7-11.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 13-10,6-12,7-11 --remove-temp &


#4. Hot-Hot
perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15HotFDR_4-5_13-6.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 4-5,13-6 --remove-temp &

perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15HotFDR_6-7_5-13.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 6-7,5-13 --remove-temp &

perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15HotFDR_5-6_4-7.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 5-6,4-7 --remove-temp

wait

perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15HotFDR_7-13_6-4.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 7-13,6-4 --remove-temp &

perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15HotFDR_13-4_7-5.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 13-4,7-5 --remove-temp &

#5. Cold-Cold
#run on blackpupa
perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15ColdFDR_9-10_8-11.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 9-10,8-11 --remove-temp &

perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15ColdFDR_11-12_10-8.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 11-12,10-8 --remove-temp &

perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15ColdFDR_10-11_9-12.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 10-11,9-12 --remove-temp &

perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15ColdFDR_12-8_11-9.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 12-8,11-9 --remove-temp &

perl /Volumes/Temp/Ray/software_ray/popoolation2-svn/cmh-test.pl --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --output /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15ColdFDR_8-9_12-10.cmhout --min-count 15 --min-coverage 15 --max-coverage 2% --population 8-9,12-10 --remove-temp

#--------------------------------------------------------------
##Part 5: Run post-cmh filtering and create concatenated cmh
#--------------------------------------------------------------

#concatenate cmh files
paste <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11.cmhout) <(cut -f27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-8_3-9.cmhout) <(cut -f27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7.cmhout) <(cut -f27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-5_3-4.cmhout) <(cut -f27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15HotF15Cold_13-10_6-12_7-11.cmhout) <(cut -f27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15ColdFDR_9-10_8-11.cmhout) <(cut -f27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15ColdFDR_11-12_10-8.cmhout) <(cut -f27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15ColdFDR_10-11_9-12.cmhout) <(cut -f27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15ColdFDR_12-8_11-9.cmhout) <(cut -f27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15ColdFDR_8-9_12-10.cmhout) <(cut -f27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15HotFDR_4-5_13-6.cmhout) <(cut -f27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15HotFDR_6-7_5-13.cmhout) <(cut -f27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15HotFDR_5-6_4-7.cmhout) <(cut -f27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15HotFDR_7-13_6-4.cmhout) <(cut -f27 /Volumes/disk5s1/AllDmelCagePops_F0toF59_F15HotFDR_13-4_7-5.cmhout) > /Volumes/disk5s1/AllDmelCagePops_MergedCMH.cmhout

#Remove unused contigs
grep -v "^U_minus" /Volumes/disk5s1/AllDmelCagePops_MergedCMH.cmhout | grep -v "^Uextra" | grep -v "^dmel" | grep -v "^gi" > /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH.cmhout

#1,797,300 SNPs prior to contig removal, 1,793,193 after removal

#Put header on merged file
cat /Volumes/Temp/Ray/mergedCMH_header.txt /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH.cmhout > /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH_header.cmhout

#Make FREQ and COV file
python /Volumes/Temp/Ray/Make_sync_with_allele_freqs_and_coverage_FAST2.py --syncext /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH_header.cmhout --common_pops 4,5,6:16,9,10:7,8:16,9,10,7,8:13,14,15:11,12:13,14,15,11,12 --common_pop_names Base,F15Hot_test,F15Hot_check,F15Hot,F15Cold_test,F15Cold_check,F15Cold --CMHpval 27,28,29,30,31,32,33,34,35,36,37,38,39,40,41 --output /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH_COVFREQ.cmhout | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,19,20,21,22,23,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,49,50,51,52,53,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80 - > /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH_COVFREQ.cmhout

###---DEPRECATED: based on old version of script which took major and minor alleles from user specified populations rather than all populations
#1st get all averages for test and check replicates separately
#python /Volumes/Temp/Ray/Make_sync_with_allele_freqs_and_coverage_FAST2.py --syncext /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH_header.cmhout --common_pops 4,5,6:16,9,10:7,8:13,14,15:11,12 --common_pop_names Base,F15Hot_test,F15Hot_check,F15Cold_test,F15Cold_check --CMHpval 27,28,29,30,31,32,33,34,35,36,37,38,39,40,41 --output /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH_COVFREQ_1.cmhout

#Next get averages for all F15 replicates combined
#python /Volumes/Temp/Ray/Make_sync_with_allele_freqs_and_coverage_FAST2.py --syncext /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH_header.cmhout --common_pops 4,5,6:16,9,10,7,8:13,14,15,11,12 --common_pop_names Base,F15Hot_all,F15Cold_all --CMHpval 27,28,29,30,31,32,33,34,35,36,37,38,39,40,41 --output /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH_COVFREQ_2.cmhout

#Last concatenate output and remove extra cols - need to take frequencies from file where minor allele is calculated from all 13 base and F15 pops counted only once (i.e. file2)
#Also remove header and replace with correct labels
#paste <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23 /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH_COVFREQ_1.cmhout) <(cut -f20,21,22,23,24,25,26,27,28,29,30,31,32,33,34 /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH_COVFREQ_2.cmhout) <(cut -f37,38,39,40,41 /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH_COVFREQ_1.cmhout) <(cut -f36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52 /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH_COVFREQ_2.cmhout) | sed '1d' | cat /Volumes/Temp/Ray/file_headers/mergedCMH_COVFREQ_header.txt - > /Volumes/disk5s1/AllDmelCagePops_MergedFilteredCMH_COVFREQ.cmhout
###---END DEPRECATED

#--------------------------------------------------------------
##Part 6: Generate igv and genes files
#--------------------------------------------------------------
#Make COVFREQ genes file
python /Volumes/Temp/Ray/Scripts_Martin/sync2vcf.py --sync /Volumes/Temp/Ray/sync/AllDmelCagePops_F0toF59_RepIndMasked.sync --snps /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ.cmhout > /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15.cmhout.vcf

cd /Volumes/Temp/Ray/Scripts_Martin/snpEff_2_0_3

#genes plus 200bp up/downstream
#v5.18
java -Xmx6g -jar ./snpEff.jar dm5.18 /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15.cmhout.vcf -inOffset 1 -outOffset 1 -ud 200 -s /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_snpeff_summary.html > /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15_dm518_UpDown200.snpeff
#v5.40
java -Xmx6g -jar ./snpEff.jar dm5.40 /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15.cmhout.vcf -inOffset 1 -outOffset 1 -ud 200 -s /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_snpeff_summary.html > /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15_dm540_UpDown200.snpeff

#Add to COVFREQ FILE
python /Volumes/Temp/Ray/Scripts_Martin/link_snpeff_203_NoAlleleCols.py --snpeff /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15_dm518_UpDown200.snpeff --input /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15.cmhout > /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15_dm518_UpDown200.cmhout.genes

python /Volumes/Temp/Ray/Scripts_Martin/link_snpeff_203_NoAlleleCols.py --snpeff /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15_dm540_UpDown200.snpeff --input /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15.cmhout > /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15_dm540_UpDown200.cmhout.genes

###DEPRECATED
##For individual comparisons
#Create igv files for CMH output
#grep -v "^U_minus" /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7.cmhout | grep -v "^Uextra" | grep -v "^dmel" | grep -v "^gi" | awk '$27!="Nan"' | sort -k 27,27g > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.cmhout
#1797300 before NaN removal and 1773627 after
#awk '{print $1"\t"$2-1"\t"$2"\tfeature\t"$27}' /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.cmhout | awk 'BEGIN{print "Chromosome\tstart\tend\tfeature\tcmh"}1' > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.cmhout.igv
#
#grep -v "^U_minus" /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11.cmhout | grep -v "^Uextra" | grep -v "^dmel" | grep -v "^gi" | awk '$27!="Nan"' | sort -k 27,27g > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.cmhout
#1797300 before NaN removal and 1779544 after
#awk '{print $1"\t"$2-1"\t"$2"\tfeature\t"$27}' /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.cmhout | awk 'BEGIN{print "Chromosome\tstart\tend\tfeature\tcmh"}1' > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.cmhout.igv

#Use SNPeff to assign gene info to all SNPs (creates gene file for all SNPs in each comparison)
#First, create vcf files for each comparison
#python /Volumes/Temp/Ray/Scripts_Martin/sync2vcf.py --sync /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --snps /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.cmhout > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.vcf &

#python /Volumes/Temp/Ray/Scripts_Martin/sync2vcf.py --sync /Volumes/disk5s1/AllDmelCagePops_F0toF59_RepIndMasked.sync --snps /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.cmhout > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.vcf
#
wait

#
#Second, generate SNPEFF output for each comparison
#cd /Volumes/Temp/Ray/Scripts_Martin/snpEff_2_0_3
#
#java -Xmx14g -jar ./snpEff.jar dm5.40 /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.vcf -inOffset 1 -outOffset 1 -ud 200 -s /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED_snpeff_summary.html > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.snpeff &
#
#java -Xmx14g -jar ./snpEff.jar dm5.40 /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.vcf -inOffset 1 -outOffset 1 -ud 200 -s /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED_snpeff_summary.html > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.snpeff
#
#wait
#
#Finally, link SNPEFF output to each file
#python /Volumes/Temp/Ray/Scripts_Martin/link_snpeff_203.py --snpeff /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.snpeff --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.cmhout > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.cmhout.genes
#
#python /Volumes/Temp/Ray/Scripts_Martin/link_snpeff_203.py --snpeff /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.snpeff --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.cmhout > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.cmhout.genes
#
##Second for top 2000 candidates
#1:Create candidate files with funny region
#head -2000 /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.cmhout > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CANDs2000.cmhout &
#
#head -2000 /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.cmhout > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CANDs2000.cmhout
#
#wait
#
#Create igv file for candidates
#awk '{print $1"\t"$2-1"\t"$2"\tfeature\t"$17}' /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CANDs2000.cmhout | awk 'BEGIN{print "Chromosome\tstart\tend\tfeature\tCMH"}1' > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CANDs2000.cmhout.igv &
#
#awk '{print $1"\t"$2-1"\t"$2"\tfeature\t"$17}' /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CANDs2000.cmhout | awk 'BEGIN{print "Chromosome\tstart\tend\tfeature\tCMH"}1' > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CANDs2000.cmhout.igv
#
#wait
#
#Create genes file for candidates
#python /Volumes/Temp/Ray/Scripts_Martin/link_snpeff_203.py --snpeff /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.snpeff --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CANDs2000.cmhout  > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CANDs2000.cmhout.genes &
#
#python /Volumes/Temp/Ray/Scripts_Martin/link_snpeff_203.py --snpeff /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.snpeff --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CANDs2000.cmhout > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CANDs2000.cmhout.genes
#
#wait
#
#2:Create candidate files with funny region
#awk '!($1=="3R" && $2>=9945225 && $2<=11521655)' /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.cmhout | head -2000 > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CANDs2000_FunnyOut.cmhout &
#
#awk '!($1=="3R" && $2>=9945225 && $2<=11521655)' /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.cmhout | head -2000 > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CANDs2000_FunnyOut.cmhout
#
#wait
#
#Create igv file for candidates
#awk '{print $1"\t"$2-1"\t"$2"\tfeature\t"$17}' /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CANDs2000_FunnyOut.cmhout | awk 'BEGIN{print "Chromosome\tstart\tend\tfeature\tCMH"}1' > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CANDs2000_FunnyOut.cmhout.igv &
#
#awk '{print $1"\t"$2-1"\t"$2"\tfeature\t"$17}' /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CANDs2000_FunnyOut.cmhout | awk 'BEGIN{print "Chromosome\tstart\tend\tfeature\tCMH"}1' > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CANDs2000_FunnyOut.cmhout.igv
#
#wait
#
#Create genes file for candidates
#python /Volumes/Temp/Ray/Scripts_Martin/link_snpeff_203.py --snpeff /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.snpeff --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CANDs2000_FunnyOut.cmhout  > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CANDs2000_FunnyOut.cmhout.genes &
#
#python /Volumes/Temp/Ray/Scripts_Martin/link_snpeff_203.py --snpeff /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.snpeff --input /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CANDs2000_FunnyOut.cmhout > /Volumes/disk5s1/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CANDs2000_FunnyOut.cmhout.genes

#--------------------------------------------------------------
##Part 6 extra: Create new COVFREQ file based on min count of 30
## get SNPeff gene info for gene only
## add to new covfreq file
#--------------------------------------------------------------
#done on vetgrid13
###Do min count = 30
python /Volumes/Temp/Ray/Scripts_Ray/Scripts_other/New_MinCount_Sync.py --sync /Volumes/Temp/Ray/CMH_files/BaseHot/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.cmhout --header "N" --cutoff 30 --output /Volumes/Temp/Ray/CMH_files/BaseHot/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED_mc30.cmhout

python /Volumes/Temp/Ray/Scripts_Ray/Scripts_other/New_MinCount_Sync.py --sync /Volumes/Temp/Ray/CMH_files/BaseCold/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.cmhout --header "N" --cutoff 30 --output /Volumes/Temp/Ray/CMH_files/BaseCold/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED_mc30.cmhout

#for making FREQCOV
python /Volumes/Temp/Ray/Scripts_Ray/Scripts_other/New_MinCount_Sync.py --sync /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_header.cmhout --header "Y" --cutoff 30 --output /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_header_mc30.cmhout

python /Volumes/Temp/Ray/Scripts_Ray/Scripts_other/Make_sync_with_allele_freqs_and_coverage_FAST2.py --syncext /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_header_mc30.cmhout --common_pops 4,5,6:16,9,10:7,8:16,9,10,7,8:13,14,15:11,12:13,14,15,11,12 --common_pop_names Base,F15Hot_test,F15Hot_check,F15Hot,F15Cold_test,F15Cold_check,F15Cold --CMHpval 27,28,29,30,31,32,33,34,35,36,37,38,39,40,41 --output /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc30.cmhout_temp 
cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,19,20,21,22,23,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,49,50,51,52,53,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80 /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc30.cmhout_temp > /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc30.cmhout

###Do min count = 45
python /Volumes/Temp/Ray/Scripts_Ray/Scripts_other/New_MinCount_Sync.py --sync /Volumes/Temp/Ray/CMH_files/BaseHot/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED.cmhout --header "N" --cutoff 45 --output /Volumes/Temp/Ray/CMH_files/BaseHot/AllDmelCagePops_F0toF59_BaseF15Hot_1-13_2-6_3-7_CMH-SORTED_mc45.cmhout

python /Volumes/Temp/Ray/Scripts_Ray/Scripts_other/New_MinCount_Sync.py --sync /Volumes/Temp/Ray/CMH_files/BaseCold/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED.cmhout --header "N" --cutoff 45 --output /Volumes/Temp/Ray/CMH_files/BaseCold/AllDmelCagePops_F0toF59_BaseF15Cold_1-10_2-12_3-11_CMH-SORTED_mc45.cmhout

#for making FREQCOV
python /Volumes/Temp/Ray/Scripts_Ray/Scripts_other/New_MinCount_Sync.py --sync /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_header.cmhout --header "Y" --cutoff 45 --output /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_header_mc45.cmhout

python /Volumes/Temp/Ray/Scripts_Ray/Scripts_other/Make_sync_with_allele_freqs_and_coverage_FAST2.py --syncext /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_header_mc45.cmhout --common_pops 4,5,6:16,9,10:7,8:16,9,10,7,8:13,14,15:11,12:13,14,15,11,12 --common_pop_names Base,F15Hot_test,F15Hot_check,F15Hot,F15Cold_test,F15Cold_check,F15Cold --CMHpval 27,28,29,30,31,32,33,34,35,36,37,38,39,40,41 --output /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc45.cmhout_temp
cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,19,20,21,22,23,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,49,50,51,52,53,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80 /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc45.cmhout_temp > /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc45.cmhout

###Do SNPeff with genes only
cd /Volumes/Temp/Ray/Scripts_Martin/snpEff_2_0_3

java -Xmx6g -jar ./snpEff.jar dm5.18 /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15.cmhout.vcf -inOffset 1 -outOffset 1 -ud 0 -s /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_snpeff_summary_dm518.html > /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15_dm518_GeneOnly.snpeff

java -Xmx6g -jar ./snpEff.jar dm5.40 /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15.cmhout.vcf -inOffset 1 -outOffset 1 -ud 0 -s /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_snpeff_summary_dm518.html > /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15_dm540_GeneOnly.snpeff

#add genes from
#genes plus/minus 200bp
python /Volumes/Temp/Ray/Scripts_Martin/link_snpeff_203_NoAlleleCols.py --snpeff /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15_dm518_UpDown200.snpeff --input /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc30.cmhout > /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc30_dm518_UpDown200.cmhout.genes
#v5.40
python /Volumes/Temp/Ray/Scripts_Martin/link_snpeff_203_NoAlleleCols.py --snpeff /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15_dm540_UpDown200.snpeff --input /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc30.cmhout > /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc30_dm540_UpDown200.cmhout.genes

#gene only - v5.18
python /Volumes/Temp/Ray/Scripts_Martin/link_snpeff_203_NoAlleleCols.py --snpeff /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15_dm518_GeneOnly.snpeff --input /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc30.cmhout > /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc30_dm518_GeneOnly.cmhout.genes
#v5.40
python /Volumes/Temp/Ray/Scripts_Martin/link_snpeff_203_NoAlleleCols.py --snpeff /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc15_dm540_GeneOnly.snpeff --input /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc30.cmhout > /Volumes/Temp/Ray/CMH_files/All_merged/AllDmelCagePops_MergedFilteredCMH_COVFREQ_mc30_dm540_GeneOnly.cmhout.genes
