#!/bin/bash -e

#SBATCH -J trim
#SBATCH -A DSMITH-BIOCLOUD
#SBATCH -o slurm-%A.out
#SBATCH -p biocloud-normal
#SBATCH --time=05:00:00

. /home/tcj25/.virtualenvs/35/bin/activate

task=$1
log=../slurm-pipeline.log

# Our data is in a directory with the same basename as ours, but up 4
# directories. The existence of this directory has been checked in
# 00-start/start.sh
dataDir=../../../../$(basename $(dirname $(/bin/pwd)))

fastq1=$dataDir/${task}_R1_001.fastq.gz
fastq2=$dataDir/${task}_R2_001.fastq.gz
out1=${task}_R1_001-trimmed.fastq.gz
out2=${task}_R2_001-trimmed.fastq.gz

function trim()
{
    echo "  AdapterRemoval started at `date`" >> $log
    srun -n 1 AdapterRemoval \
      --basename $task \
      --file1 $fastq1 \
      --file2 $fastq2 \
      --output1 $out1 \
      --output2 $out2 \
      --gzip \
      --minlength 20 \
      --trimns \
      --trimqualities

    echo "  Sleeping 5 seconds to give the filesystem a chance to settle." >> $log
    sleep 5
    echo "  AdapterRemoval stopped at `date`" >> $log
}


echo "01-trim on task $task started at `date`" >> $log
echo "  fastq1 is $fastq1" >> $log
echo "  fastq2 is $fastq2" >> $log

if [ $SP_SIMULATE = "0" ]
then
    echo "  This is not a simulation." >> $log
    if [ -f $out1 -a -f $out2 ]
    then
        if [ $SP_FORCE = "1" ]
        then
            echo "  Pre-existing output files $out1 and $out2 exist, but --force was used. Overwriting." >> $log
            trim
        else
            echo "  Will not overwrite pre-existing output files $out1 and $out2. Use --force to make me." >> $log
        fi
    else
        echo "  Pre-existing output files $out1 and $out2 do not both exist. Trimming." >> $log
        trim
    fi
else
    echo "  This is a simulation." >> $log
fi

echo "01-trim on task $task stopped at `date`" >> $log
echo >> $log
