#!/bin/bash
#$ -l highmem,highp,h_rt=48:00:00,h_data=48G
#$ -N build_crux_db_CO1_12042018
#$ -cwd
#$ -m bea
#$ -M meixilin
#$ -o /u/project/rwayne/meixilin/capture_array/metabar/anacapa/qsub_log/build_crux_db_CO1_12042018.out.txt
#$ -e /u/project/rwayne/meixilin/capture_array/metabar/anacapa/qsub_log/build_crux_db_CO1_12042018.err.txt

# echo job info on joblog:
echo "Job $JOB_ID started on:    " `hostname -s`
echo "Job $JOB_ID started on:    " `date `
echo " "

source /u/local/Modules/default/init/bash

# define variables
# please note that this primer was full length.
PRIMER=CO1
PRIMER_F=TNTCNACNAAYCAYAARGAYATTGG
PRIMER_R=TANACYTCNGGRTGNCCRAARAAYCA
MIN_LEN=50
MAX_LEN=1000
MAX_ERR=3

WORK=/u/project/rwayne/meixilin/capture_array/metabar/anacapa
CRUX_DB=/u/project/rwayne/software/Crux/crux_db
OUTDIR=/u/project/rwayne/meixilin/crux_db/${PRIMER}

# make directory if not made before
mkdir -p ${OUTDIR}

cd ${WORK}

bash ${CRUX_DB}/crux.sh \
-u meixilin \
-n ${PRIMER} \
-f ${PRIMER_F} \
-r ${PRIMER_R} \
-s ${MIN_LEN} \
-m ${MAX_LEN} \
-o ${OUTDIR} \
-d ${CRUX_DB} \
-e ${MAX_ERR}

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
