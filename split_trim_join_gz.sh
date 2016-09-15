BASE_NAME=$1
FQ_DIR=$2
POPO=/Volumes/Temp/Lukas/Tools/popoolation/
if [ ! $FQ_DIR ]; then FQ_DIR=..; fi 
gunzip -c $FQ_DIR/${BASE_NAME}1.fastq.gz | split -a 1 -l 100000000 - ${BASE_NAME}_1_ &
gunzip -c $FQ_DIR/${BASE_NAME}2.fastq.gz | split -a 1 -l 100000000 - ${BASE_NAME}_2_ 
wait
PIDS=()
for i in ${BASE_NAME}_1_??; do { perl $POPO/basic-pipeline/trim-fastq.pl --input1 $i --input2 ${i/_1_/_2_} --output ${i/_1_/_trim_} --quality-threshold 20 --fastq-type sanger --discard-internal-N --min-length 50 --no-5p-trim 2>&1 >> ${i/_1_/_trim_}.log  & }; PIDS=( ${PIDS[@]} $! ) ; done

for i in ${PIDS[@]}; do wait $i; done

rm -f ${BASE_NAME}_[12]_??

cat ${BASE_NAME}_trim_??_1.gz | gunzip -c  | gzip -c > ${BASE_NAME}_trimmed_1.gz &
cat ${BASE_NAME}_trim_??_2.gz | gunzip -c  | gzip -c > ${BASE_NAME}_trimmed_2.gz 
wait

rm -f ${BASE_NAME}_trim_??_[12].gz

