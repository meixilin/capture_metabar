#!/usr/bin/env bash

# define the variables
OUT=/u/home/m/meixilin/project-rwayne/capture_array/metabar/anacapa/anacapa_out
MB=12SV5
FIL=l30
PERCENT="40 50 60 70 80 90 95 100"
WORK=${OUT}/${MB}/${MB}_taxonomy_tables/Summary_by_percent_confidence_${FIL}
mkdir -p ${WORK}

mkdir -p ${OUT}/Run_info/taxon_summary_logs_${FIL}
cd ${OUT}/${MB}/${MB}_taxonomy_tables

for per in ${PERCENT}
do
  mkdir -p ${WORK}/${per}
  echo "${per}"
  python ${DB}/scripts/reformat_summary_for_r.py ${MB}_ASV_taxonomy_detailed_l30.txt  ${WORK}/${per}/${MB}_ASV_raw_taxonomy_${per}.txt ${per}
  Rscript --vanilla ${DB}/scripts/sum_blca_for_R_by_taxon.R  ${WORK}/${per}/${MB}_ASV_raw_taxonomy_${per}.txt ${MB} ${WORK}/${per}/${MB}_ASV_sum_by_taxonomy_${per}.txt 2>> ${OUT}/Run_info/taxon_summary_logs_${FIL}/${MB}_ASV_sum_by_taxonomy.log.txt
done
date
