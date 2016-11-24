#!/bin/bash -e

task=$1

log=../$task.log

echo "01-trim sbatch.sh running at `date`" >> $log
echo "task is $task" >> $log
echo "dependencies are $SP_DEPENDENCY_ARG" >> $log
echo >> $log

jobid=`sbatch -n 1 trim.sh $task | cut -f4 -d' '`
echo "TASK: $task $jobid"
