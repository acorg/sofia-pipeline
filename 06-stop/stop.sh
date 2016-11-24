#!/bin/bash -e

#SBATCH -J stop
#SBATCH -A DSMITH-BIOCLOUD
#SBATCH -o slurm-%A.out
#SBATCH -p biocloud-normal
#SBATCH --time=5:00:00

. /home/tcj25/.virtualenvs/35/bin/activate

log=../slurm-pipeline.log

echo "SLURM pipeline finished at `date`" >> $log

touch ../slurm-pipeline.done
