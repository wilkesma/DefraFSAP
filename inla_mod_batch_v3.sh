#!/bin/bash
REPO=/home/mw22803/Defra/R/github
PROCESSED=$REPO/data/processed

# File: inla_mod_batch_v3.sh
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N inla_mod_final_v3.R
#$ -m n
#$ -t 1-266
#$ -l h_vmem=10G
cd $REPO
R CMD BATCH $REPO/inla_mod_final_v3.R $PROCESSED/inla_mods/logs/inla_mod_final_v3.Rout.$SGE_TASK_ID