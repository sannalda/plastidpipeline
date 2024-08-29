#!/bin/bash

#SBATCH --mail-user=sjannalda@zedat.fu-berlin.de
#SBATCH --job-name=Downsampling
#SBATCH --mail-type=end
#SBATCH --mem=50000
#SBATCH --time=0-12:00:00
#SBATCH --qos=standard
#SBATCH --output=Downsampling.out
#SBATCH --error=Downsampling.err
#SBATCH --cpus-per-task=16

module load BBMap

DATADIR=/scratch/reichel/pipelinetest/Spain/FASTQ
WORKDIR=/scratch/sjannalda/bgbm/data/DownloadedData/Spain
cat ${DATADIR}/HFYM3DSX7_1_379UDI-idt-UMI_1.fastq.gz ${DATADIR}HFYM3DSX7_2_379UDI-idt-UMI_1.fastq.gz > ${WORKDIR}/merged1_1.fastq.gz
cat ${DATADIR}/HFYM3DSX7_3_379UDI-idt-UMI_1.fastq.gz ${DATADIR}HFYM3DSX7_4_379UDI-idt-UMI_1.fastq.gz > ${WORKDIR}/merged2_1.fastq.gz
cat ${WORKDIR}/merged1_1.fastq.gz ${WORKDIR}/merged2_1.fastq.gz > ${WORKDIR}/merged_1.fastq.gz 
rm ${WORKDIR}/merged1_1.fastq.gz
rm ${WORKDIR}/merged2_1.fastq.gz

cat ${DATADIR}/HFYM3DSX7_1_379UDI-idt-UMI_2.fastq.gz ${DATADIR}HFYM3DSX7_2_379UDI-idt-UMI_2.fastq.gz > ${WORKDIR}/merged1_2.fastq.gz
cat ${DATADIR}/HFYM3DSX7_3_379UDI-idt-UMI_2.fastq.gz ${DATADIR}HFYM3DSX7_4_379UDI-idt-UMI_2.fastq.gz > ${WORKDIR}/merged2_2.fastq.gz
cat ${WORKDIR}/merged1_2.fastq.gz ${WORKDIR}/merged2_2.fastq.gz > ${WORKDIR}/merged_2.fastq.gz 
rm ${WORKDIR}/merged1_2.fastq.gz
rm ${WORKDIR}/merged2_2.fastq.gz

INF1=${WORKDIR}/merged_1.fastq.gz
INF2=${WORKDIR}/merged_2.fastq.gz

INF1_DS=${WORKDIR}/merged_1.downsampled.fastq.gz
INF2_DS=${WORKDIR}/merged_2.downsampled.fastq.gz
bbnorm.sh in=$INF1 in2=$INF2 out=${INF1_DS} out2=${INF2_DS} target=200 min=5 usejni=t threads=16

#seqtk sample -s100 ${WORKDIR}/merged_1.fastq.gz 0.10 > /scratch/sjannalda/projects/arnica/test/sample_10/sub2.fq

#seqtk sample -s100 /scratch/sjannalda/projects/arnica/test/data/Am0057_1.paired.qualtrimmed.downsampled.small.sorted.fastq.gz 0.10 > /scratch/sjannalda/projects/arnica/test/sample_10/sub1.fq