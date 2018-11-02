module load python/2.7.13

ECOPCR=/u/home/m/meixilin/ecoPCR/src/ecoPCR
ECOPCRDB=/u/home/m/meixilin/project-rwayne/embl_r133/embl_r133
WORK=/u/home/m/meixilin/project-rwayne/capture_array/metabar/primer_info
primers=(X12SV5 Mamm2)

cd ${WORK}

for ii in "${primers[@]}"; do
        obigrep -d ${ECOPCRDB} --require-rank=species --require-rank=genus --require-rank=family ${ii}_0_250_e3.ecopcr > ${ii}_0_250_e3_clean.fasta
        obiuniq -d ${ECOPCRDB} ${ii}_0_250_e3_clean.fasta > ${ii}_0_250_e3_clean_uniq.fasta
        obigrep -d ${ECOPCRDB} --require-rank=family ${ii}_0_250_e3_clean_uniq.fasta > ${ii}_cuc.fasta
        obiannotate --uniq-id ${ii}_cuc.fasta > ${ii}_db.fasta
done

# then manually renamed files
