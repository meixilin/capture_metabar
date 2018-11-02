SEQSOURCE=/u/home/m/meixilin/project-rwayne/rawrawseq/RW_ML-98671768/FASTQ_Generation_2018-10-15_17_26_24Z-130565915
WORK=/u/home/m/meixilin/project-rwayne/capture_array/metabar

mkdir -p ${WORK}/raw_seq

INSEQ=${WORK}/raw_seq

# copy and rename files

# filename list
cd ${SEQSOURCE}
seqnames=(extract_blank FOD_DNA_1 MVP_DNA_2 MVP_DNA_3 pcr_blank positive_mammal PPB_DNA_1 SMM_DNA_1 SMV_DNA_1 SSA_DNA_1)
for ii in "${seqnames[@]}"; do
        cp ${ii}*/*R1_001.fastq.gz ${INSEQ}/${ii}_R1_001.fastq.gz
        cp ${ii}*/*R2_001.fastq.gz ${INSEQ}/${ii}_R2_001.fastq.gz
done

# untar all fastq files
gunzip ${INSEQ}/*.fastq.gz
