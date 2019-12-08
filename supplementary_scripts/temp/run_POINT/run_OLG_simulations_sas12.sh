#! /usr/bin/bash

#PBS -V
#PBS -N OLG_R2_100k_sas12
#PBS -l ncpus=25,walltime=999:00:00
#PBS -M cwnelson88@gmail.com
#PBS -m abe
#PBS -j oe

### Remember to match ncpus to NCPUS value below
NBOOTSTRAPS=1000
NCPUS=25
FRAME_TYPE=sas12
DISTANCE=0.05854733
R=2
WORK=${PBS_O_WORKDIR}

# Enter working directory if not there already (may be redundant)
cd ${WORK}

# Run program
date # to get a sense of how long it took
/nas3/cnelson/bin/OLGenie_simulations_pipeline_POINT_100k.pl ${NBOOTSTRAPS} ${NCPUS} ${FRAME_TYPE} ${DISTANCE} ${R} 2>&1 | tee OLGenie_simulations_${FRAME_TYPE}.out
date

