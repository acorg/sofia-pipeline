#!/bin/bash -e

. /home/tcj25/.virtualenvs/35/bin/activate

log=../slurm-pipeline.log

# Remove the marker file that indicates when a job is fully complete.
rm -f ../slurm-pipeline.done

echo "SLURM pipeline started at `date`" >> $log

for fastq in "$@"
do
    task=`echo $fastq | cut -f1 -d.`
    echo >> $log
    echo "  FASTQ file $fastq" >> $log
    echo "  task name  $task" >> $log
    # Emit task names (without job ids as this step does not start any
    # SLURM jobs).
    echo "TASK: $task"
done

echo >> $log
