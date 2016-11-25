#!/bin/bash -e

task=$1

log=../sbatch.log

echo "02-flash sbatch.sh running at `date`" >> $log
echo "task is $task" >> $log
echo "dependencies are $SP_DEPENDENCY_ARG" >> $log
echo >> $log

jobid=`sbatch -n 1 flash.sh $task | cut -f4 -d' '`
echo "TASK: $task $jobid"
