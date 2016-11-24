#!/bin/bash -e

task=$1
log=../$task.log

echo "02-map sbatch.sh running at `date`" >> $log
echo "task is $task" >> $log
echo "dependencies are $SP_DEPENDENCY_ARG" >> $log
echo >> $log

# Request an exclusive machine because we're going to tell bwa and samtools
# to use 24 threads.
jobid=`sbatch -n 1 --exclusive $SP_DEPENDENCY_ARG submit.sh $task | cut -f4 -d' '`
echo "TASK: $task $jobid"
