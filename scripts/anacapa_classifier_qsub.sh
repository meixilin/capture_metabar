#!/bin/bash
#$ -l highp,h_rt=04:00:00,h_data=24G
#$ -N capture_anacapa_classifier
#$ -cwd
#$ -m bea
#$ -M meixilin
#$ -o /u/project/rwayne/meixilin/capture_array/metabar/anacapa/qsub_log/capture_anacapa_classifier_11272018.out.txt
#$ -e /u/project/rwayne/meixilin/capture_array/metabar/anacapa/qsub_log/capture_anacapa_classifier_11272018.err.txt

# echo job info on joblog:
echo "Job $JOB_ID started on:    " `hostname -s`
echo "Job $JOB_ID started on:    " `date `
echo " "

source /u/local/Modules/default/init/bash
# define variables
WORK=/u/project/rwayne/meixilin/capture_array/metabar/anacapa
ANACAPA_CL=/u/home/m/meixilin/Anacapa/Anacapa_db/anacapa_classifier.sh
ANACAPA_DB=/u/home/m/meixilin/Anacapa/Anacapa_db
OUTDIR=${WORK}/anacapa_out

cd ${WORK}

bash ${ANACAPA_CL} \
-u meixilin \
-o ${OUTDIR} \
-d ${ANACAPA_DB}


# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
