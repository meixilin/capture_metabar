#!/bin/bash
#$ -l highp,h_rt=05:00:00,h_data=16G
#$ -N obitools_on_qcoutput_10302018
#$ -cwd #current working directory
#$ -m bea # send a email when it starts and done
#$ -o /u/project/rwayne/meixilin/capture_array/metabar/obi_analysis/qsub_log/obitools_on_qcoutput_10302018.out.txt
#$ -e /u/project/rwayne/meixilin/capture_array/metabar/obi_analysis/qsub_log/obitools_on_qcoutput_10302018.err.txt
#$ -M meixilin # username on hoffman

# echo job info on joblog:
echo "Job $JOB_ID started on:    " `hostname -s`
echo "Job $JOB_ID started on:    " `date `
echo " "

# load the job environment:
source /u/local/Modules/default/init/bash
module load python/2.7.13

WORK=/u/home/m/meixilin/project-rwayne/capture_array/metabar/obi_analysis
cd $WORK

OUTDIR=${WORK}/obitools_out
ECOPCR=/u/home/m/meixilin/ecoPCR/src/ecoPCR
ECOPCRDB=/u/home/m/meixilin/project-rwayne/embl_r133/embl_r133

# list of primers to analyse
primers=(12SV5 MAMM2)

for ii in "${primers[@]}"; do
        PRIMER=${ii}
        INSEQ=/u/project/rwayne/meixilin/capture_array/metabar/anacapa/anacapa_out/${PRIMER}/${PRIMER}dada2_out/nochim_merged${PRIMER}.fasta
        ECOTAGDB=/u/home/m/meixilin/project-rwayne/capture_array/metabar/primer_info/${PRIMER}_db.fasta

## assign taxonomy using ecotags
ecotag -d $ECOPCRDB -R $ECOTAGDB ${INSEQ} > $OUTDIR/nochim_merged${PRIMER}.tag.fasta

# generate table
obitab -o $OUTDIR/nochim_merged${PRIMER}.tag.fasta > $OUTDIR/nochim_merged${PRIMER}.tag.tsv
done

echo "4. assign each seq to a taxon ended on" `date `

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
