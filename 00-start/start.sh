#!/bin/bash -e

. /home/tcj25/.virtualenvs/35/bin/activate

log=../slurm-pipeline.log

# Remove the marker file that indicates when a job is fully complete.
rm -f ../slurm-pipeline.done

echo "SLURM pipeline started at `date`" >> $log

echo >> $log
echo "00-start started at `date`" >> $log

# Assume our data is in a directory with the same basename as ours, but up
# 4 directories.
dataDir=../../../../$(basename $(/bin/pwd))

if [ ! -d $dataDir ]
then
    echo "  Data directory '$dataDir' does not exist!" >> $log
    exit 1
fi

#
# Make sure there are only 2 FASTQ files in the data directory, and that
# their names form a pair in the expected way (R1 and R2 must be the only
# difference in the names). E.g.:
#
# 141110-79_S79_L001_R1_001.fastq.gz
# 141110-79_S79_L001_R2_001.fastq.gz
#

IFS=$'\n'
fastq=($(ls $dataDir/*.fastq.gz | sed -e 's:.*/::' | cut -f1 -d.))
unset IFS

case ${#fastq[*]} in
    2) ;;
    *) echo "  Unexpected FASTQ file list: ${fastq[*]}" >> $log; exit 1;;
esac

expectedSecond=`echo ${fastq[0]} | sed -e s/_R1_/_R2_/`

if [ ${fastq[1]} != $expectedSecond ]
then
    echo "  Second FASTQ file name '${fastq[1]}' does not match expected name ($expectedSecond)" >> $log
    exit 1
fi

# The task name is the common prefix of the FASTQ file names, up to the R1.
task=`echo $fastq[0] | sed -e s/_R1_.*//`

echo "  FASTQ files: ${fastq[*]}" >> $log
echo "  task name $task" >> $log

# Emit task name (without job ids as this step does not start any
# SLURM jobs).
echo "TASK: $task"

echo "00-start on task $task stopped at `date`" >> $log
echo >> $log
