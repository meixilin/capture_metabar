#!/usr/bin/env bash

# This is a script to link databases
DB_18S=/u/project/rwayne/software/Anacapa/Anacapa_db/18S
DB_CO1=/u/project/rwayne/software/Anacapa/Anacapa_db/CO1
DB_PITS=/u/project/rwayne/software/Anacapa/Anacapa_db/PITS
DB_12SV5=/u/project/rwayne/meixilin/crux_db/12SV5/12SV5_db_filtered
DB_MAMM2=/u/project/rwayne/meixilin/crux_db/MAMM2/MAMM2_db_filtered
ANACAPA_DB=/u/home/m/meixilin/Anacapa/Anacapa_db

ln -s ${DB_18S} ${ANACAPA_DB}
ln -s ${DB_CO1} ${ANACAPA_DB}
ln -s ${DB_PITS} ${ANACAPA_DB}
ln -s ${DB_12SV5} ${ANACAPA_DB}
ln -s ${DB_MAMM2} ${ANACAPA_DB}

# rename the files
mv ${ANACAPA_DB}/12SV5_db_filtered ${ANACAPA_DB}/12SV5
mv ${ANACAPA_DB}/MAMM2_db_filtered ${ANACAPA_DB}/MAMM2
