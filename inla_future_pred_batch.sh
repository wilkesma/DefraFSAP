#!/bin/bash

REPO=$(cd "$(dirname "$0")" && pwd)
PROCESSED=$REPO/data/processed

# File: inla_future_pred_batch.sh
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N inla_future_pred.R
#$ -m n
#$ -t 1-266
#$ -l h_vmem=32G,mem_free=200G
#$ -l hostname="!compute-0-36&!compute-1-32&!compute-1-39"
#$ -o /home/mw22803/Defra/R/inla_future_pred/logs/
cd $REPO
R CMD BATCH --no-restore $REPO/inla_future_pred.R $PROCESSED/inla_future_pred/logs/inla_future_pred.Rout.$SGE_TASK_ID