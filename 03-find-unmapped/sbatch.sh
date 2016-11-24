#!/bin/bash -e

task=$1
log=../$task.log

echo "03-find-unmapped sbatch.sh running at `date`" >> $log
echo "task is $task" >> $log
echo "dependencies are $SP_DEPENDENCY_ARG" >> $log
echo >> $log

jobid=`sbatch -n 1 $SP_DEPENDENCY_ARG submit.sh $task | cut -f4 -d' '`
echo "TASK: $task $jobid"
