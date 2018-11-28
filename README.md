> Author: Meixi Lin 
>
> Date: November 2018
> 

This is a repo using a combination of `anacapa` and `obitools` to generate taxonomical classification for some NGS sequences. 

# scripts 

1. move files: [copy_rename_seq_10292018.sh](scripts/copy_rename_seq_10292018.sh)
2. qc and demultiplex using anacapa: [anacapa_qc_dada2_qsub.sh](scripts/anacapa_qc_dada2_qsub.sh)
3. ecotag classification: 
   1.  build & clean ref database: 
       1.  [test_12SV5_ecopcr.sh](scripts/test_12SV5_ecopcr.sh) 
       2.  [clean_ecopcr.sh](scripts/clean_ecopcr.sh)
   2.  run classification: [obitools_on_qcoutput_10302018.sh](scripts/obitools_on_qcoutput_10302018.sh)
4. anacapa classification: 
   1. build `crux` database: [build_crux_db_MAMM2.sh](scripts/build_crux_db_MAMM2.sh) & a similar one 
   2. create symlink for database used: [create_symlink_crux_db.sh](scripts/create_symlink_crux_db.sh)
   3. run classification: [anacapa_classifier_qsub.sh](scripts/anacapa_classifier_qsub.sh) 

# The primer sequence info

please refer to the config file

# ref

1. for obitools: https://git.metabarcoding.org/explore/projects
2. for anacapa: https://github.com/limey-bean/Anacapa

# Side notes

> Date: 11/27/2018

Remember to set `BLCAperMINlen="0.8"`in the `anacapa_vars.sh`

