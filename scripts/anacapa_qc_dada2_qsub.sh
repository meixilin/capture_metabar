#!/bin/bash
#$ -l highp,h_rt=04:00:00,h_data=24G
#$ -N capture_QC_dada2
#$ -cwd
#$ -m bea
#$ -M meixilin
#$ -o /u/project/rwayne/meixilin/capture_array/metabar/anacapa/qsub_log/capture_QC_dada2_10292018.out.txt
#$ -e /u/project/rwayne/meixilin/capture_array/metabar/anacapa/qsub_log/capture_QC_dada2_10292018.err.txt

# echo job info on joblog:
echo "Job $JOB_ID started on:    " `hostname -s`
echo "Job $JOB_ID started on:    " `date `
echo " "

source /u/local/Modules/default/init/bash
# define variables
WORK=/u/project/rwayne/meixilin/capture_array/metabar/anacapa
ANACAPA_QC=/u/home/m/meixilin/Anacapa/Anacapa_db/anacapa_QC_dada2.sh
ANACAPA_DB=/u/home/m/meixilin/Anacapa/Anacapa_db
INSEQ=${WORK/anacapa}raw_seq
OUTDIR=${WORK}/anacapa_out
PRIMER_F=${WORK}/config/forward_primers.txt
PRIMER_R=${WORK}/config/reverse_primers.txt
PRIMER_LEN=${WORK}/config/metabarcode_loci_min_merge_length.txt

cd ${WORK}

bash ${ANACAPA_QC} \
-i ${INSEQ} \
-o ${OUTDIR} \
-d ${ANACAPA_DB} \
-u meixilin \
-f ${PRIMER_F} \
-r ${PRIMER_R} \
-e ${PRIMER_LEN} \
-a nextera \
-t MiSeq \
-m 30 \
-x 0 \
-y 0 \
-g


# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
