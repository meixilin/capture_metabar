#!/bin/bash
#$ -l highp,h_rt=06:00:00,h_data=24G
#$ -N X12SV5_ecopcr_10302018
#$ -cwd #current working directory
#$ -m bea # send a email when it starts and done
#$ -o /u/home/m/meixilin/project-rwayne/capture_array/metabar/X12SV5_ecopcr_10302018.out.txt
#$ -e /u/home/m/meixilin/project-rwayne/capture_array/metabar/X12SV5_ecopcr_10302018.err.txt
#$ -M meixilin # username on hoffman

ECOPCR=/u/home/m/meixilin/ecoPCR/src/ecoPCR
ECOPCRDB=/u/home/m/meixilin/project-rwayne/embl_r133/embl_r133
FORWARD_SEQ=TTAGATACCCCACTATGC
REVERSE_SEQ=TAGAACAGGCTCCTCTAG
PRIMER_NAME=X12SV5
OUTDIR=/u/home/m/meixilin/project-rwayne/capture_array/metabar

# perform ecopcr on the embl_r133 database
${ECOPCR} -d ${ECOPCRDB} -l 0 -L 250 -e 3 ${FORWARD_SEQ} ${REVERSE_SEQ} > ${OUTDIR}/${PRIMER_NAME}_0_250_e3.ecopcr
