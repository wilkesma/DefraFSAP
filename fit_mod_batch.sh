#!/bin/bash
REPO=/home/mw22803/Defra/R/github
PROCESSED=$REPO/data/processed

# File: fit_mod_batch.sh
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N fit_mod.R
#$ -m n
#$ -t 1-56
#$ -tc 56
#$ -l h_vmem=250G
cd $REPO
R CMD BATCH $REPO/fit_mod.R $PROCESSED/mods/fit_mod.Rout.$SGE_TASK_ID